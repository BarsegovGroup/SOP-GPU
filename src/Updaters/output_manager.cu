/*
 * output_manager.cu
 *
 *  Created on: Apr 7, 2010
 *      Author: zhmurov
 */

#include "../def_param.h"
#include "../gsop.h"
#include "../param_initializer.h"

#include "../IO/configreader.h"
#include "../Util/mystl.h"
#include <string.h>
#include <stdio.h>
#include <cuda.h>
#include "output_manager.h"
#include "output_manager_kernel.cu"

int mode;

void createOutputManager(){
	updaters[updatersCount] = new OutputManager();
	updatersCount++;
}

OutputManager::OutputManager(){
    FILE *dat_file;
	printf("Initializing output manager...\n");
	this->name = "DAT output";
	this->frequency = getIntegerParameter(OUTPUT_FREQUENCY_STRING, DEFAULT_OUTPUT_FREQUENCY, 1);

	// Parameters for std out
	this->outputWidth = getIntegerParameter(OUTPUT_DATA_WIDTH_STRING, DEFAULT_OUTPUT_DATA_WIDTH, 1);
	this->printRuns = getIntegerParameter(OUTPUT_PRINT_RUNS_STRING, DEFAULT_OUTPUT_PRINT_RUNS, 1);

	// If the gyration radius is needed
	this->computeRgFlag = getYesNoParameter(OUTPUT_COMPUTE_RG_STRING, DEFAULT_OUTPUT_COMPUTE_RG, 1);

	// Parameters to define native contacts
	this->R_limit = getFloatParameter(COVALENT_R_LIMIT_STRING, DEFAULT_COVALENT_R_LIMIT, 1);
	this->R_limit_bond = getFloatParameter(NATIVE_R_LIMIT_BOND_STRING, DEFAULT_NATIVE_R_LIMIT_BOND, 1);

	// Deprecated 'mode capsid' (if defined, the gyration radius is not computed due to the performance issues).
	char modeString[100];
	getParameter(modeString, "mode", "", 1);
	if(strcmp(modeString, "capsid") == 0){
		mode = MODE_CAPSID;
		printf("WARNING: 'mode capsid' is deprecated. Use '%s' parameter to turn the Rg computation on/off.\n", OUTPUT_COMPUTE_RG_STRING);
		this->computeRgFlag = false;
	} else {
		mode = 0;
	}

	// Preparing files
	int traj;
    std::string dat_filename = getMaskedParameterAs<std::string>(OUTPUT_FILENAME_STRING, DEFAULT_OUTPUT_FILENAME);
	dat_filenames.resize(gsop.Ntr);
	for(traj = 0; traj < gsop.Ntr; traj++){
        dat_filenames[traj] = string_replace(dat_filename, "<run>", traj+gsop.firstrun);
		dat_file = fopen(dat_filenames[traj].c_str(), "w");
		fclose(dat_file);
	}
	printf("Output will be saved in '%s'.\n", dat_filename.c_str());

	// Allocatng memory for gyration radii
    this->h_rgs = (float*)calloc(gsop.aminoCount, sizeof(float));
    cudaMalloc((void**)&this->d_rgs, gsop.aminoCount*sizeof(float));
    this->rgs = (float*)calloc(gsop.Ntr, sizeof(float));
    this->temperatures = (float*)calloc(gsop.Ntr, sizeof(float));
    this->nativeCounts = (int*)calloc(gsop.Ntr, sizeof(int));

	printf("Done initializing output manager...\n");
}

OutputManager::~OutputManager(){
}

void OutputManager::update(){
	if(step % this->frequency == 0){

		// Compute all values for the output
		int p, id;
		for(p = 0; p < potentialsCount; p++){
			for(id = 0; id < potentials[p]->getEnergiesCount(); id++){
				potentials[p]->computeEnergy(id);
			}
		}

		this->computeTemperatures();
		this->computeNativeCounts();

		if(this->computeRgFlag){
			this->computeRgs();
		}

		// Some info for standard output
		char runstring[20];
		if(gsop.Ntr == 1){
			sprintf(runstring, "Run %d", gsop.firstrun);
		} else {
			sprintf(runstring, "%d runs", gsop.Ntr);
		}
		printf("Writing output at step %ld of %ld. %s on %d [%s].\n", step, numsteps, runstring, gsop.deviceId, deviceProp.name);
		this->printDataTitleToScreen();

		// Output for all trajectories
		int traj;
		for(traj = 0; traj < gsop.Ntr; traj++){

			// Computing values
			this->computeTEAeps(traj);

			// Standard output
			if(traj < this->printRuns){
				this->printDataToScreen(traj);
			}
			if(traj == this->printRuns){
				this->printDataLdotsToScreen();
			}

			// .dat file output
			FILE* dat_file = fopen(dat_filenames[traj].c_str(), "a");
			this->printDataToFile(dat_file, traj);
			fclose(dat_file);
		}
		printf("Done writing output.\n");
    }
}


/*
 * Computation of the average temperature (based on atoms displacement)
 * Averaging is on both particles and time.
 * The temperature values are set to zero on a GPU once the temperature is computed.
 */
void OutputManager::computeTemperatures(){
	copyTemperatureDeviceToHost();
	int i, traj;
	double tempav = 0.0;
	for(traj = 0; traj < gsop.Ntr; traj++){
		for(i = 0; i < sop.aminoCount; i++){
			tempav += gsop.h_T[traj*sop.aminoCount + i];
		}
		tempav /= ((double)(sop.aminoCount*this->frequency));
		this->temperatures[traj] = tempav;
	}
    const int blockSize = gsop.blockSize;
    const int blockNum = gsop.aminoCount/blockSize + 1;
    reset_temperature<<<blockNum, blockSize>>>();
}

/*
 * Number of native contacts that are intact.
 * The contact is considered alive if one of the following conditions is satisfied:
 * 1. The bead-to-bead distance is lower than the R_limit_bond (default is 8A).
 * 2. Absolute value of the difference of the current and native bead-to-bead distances is lower than R_limit (default is 2A).
 */
void OutputManager::computeNativeCounts(){
	copyCoordDeviceToHost();
	int i, j, k, traj;
	for(traj = 0; traj < gsop.Ntr; traj++){
		int natCount = 0;
		for(k = 0; k < sop.nativeCount; k++){
			i = sop.natives[k].i;
			j = sop.natives[k].j;
			float dr;
			float3 r;
			r.x = gsop.h_coord[traj*sop.aminoCount + i].x - gsop.h_coord[traj*sop.aminoCount + j].x;
			r.y = gsop.h_coord[traj*sop.aminoCount + i].y - gsop.h_coord[traj*sop.aminoCount + j].y;
			r.z = gsop.h_coord[traj*sop.aminoCount + i].z - gsop.h_coord[traj*sop.aminoCount + j].z;
			dr = sqrtf(r.x*r.x+r.y*r.y+r.z*r.z);
			if(dr < this->R_limit_bond || fabs(sop.natives[k].r0 - dr) <= this->R_limit){
				natCount ++;
			}
		}
		this->nativeCounts[traj] = natCount;
	}
}

/*
 * Evaluation of the gyration radii. Done on the GPU (see computeRg_kernel).
 */
void OutputManager::computeRgs(){
	const int blockSize = gsop.blockSize;
	const int blockNum = gsop.aminoCount/blockSize + 1;
	checkCUDAError();
	computeRg_kernel<<<blockNum, blockSize>>>(d_rgs);
	checkCUDAError();
	cudaMemcpy(this->h_rgs, this->d_rgs, gsop.aminoCount*sizeof(float), cudaMemcpyDeviceToHost);
	checkCUDAError();
	double rg = 0.0f;
	int i, traj;
	for(traj = 0; traj < gsop.Ntr; traj++){
		for(i = traj*sop.aminoCount; i < (traj+1)*sop.aminoCount; i++){
			rg += h_rgs[i];
		}
		rg = sqrt(rg/(double)(sop.aminoCount*sop.aminoCount));
		this->rgs[traj] = rg;
	}
}

/*
 * Hydrodynamics EPS value
 */
void OutputManager::computeTEAeps(int traj){
	if (integratorTea) {
        TeaIntegrator *tea = (TeaIntegrator*)integrator;
        if (!tea->exact){
    		outputData.tea_eps = tea->h_epsilon[traj];
        }
	}else{
		outputData.tea_eps = 0./0.;
	}
}

/*
 * Print energy titles to the standard output
 */
void OutputManager::printDataTitleToScreen() const{
	int p, id;
	printf("%*s%*s%*s",
			this->outputWidth, "Run",
			this->outputWidth, "TimeStep",
			this->outputWidth, "Temperature");
	for(p = 0; p < potentialsCount; p++){
		if(potentials[p]->getEnergiesCount() == 1){
			printf("%*s", this->outputWidth, potentials[p]->name.c_str());
		} else {
			for(id = 0; id < potentials[p]->getEnergiesCount(); id++){
				printf("%*s%d", this->outputWidth-1, potentials[p]->name.c_str(), id+1);
			}
		}
	}
	printf("%*s%*s%*s",
			this->outputWidth, "Total",
			this->outputWidth, "Native#",
			this->outputWidth, "Rg");
	if(integratorTea){
		printf("%*s",
				this->outputWidth, "TEA-eps");
	}
	printf("\n");
}

/*
 * Print energy values for a given trajectory
 */
void OutputManager::printDataToScreen(int traj) const{
	int p, id;
	printf("%*d%*ld%*f",
			this->outputWidth, traj + gsop.firstrun,
			this->outputWidth, step,
			this->outputWidth, this->temperatures[traj]);

	double total = 0.0;
	for(p = 0; p < potentialsCount; p++){
		for(id = 0; id < potentials[p]->getEnergiesCount(); id++){
			float pot = potentials[p]->getEnergy(traj, id);
			printf("%*f", this->outputWidth, pot);
			total += pot;
		}
	}
	printf("%*f", this->outputWidth, total);
	printf("%*d", this->outputWidth, this->nativeCounts[traj]);
	printf("%*f", this->outputWidth, this->rgs[traj]);

	if(integratorTea){
		printf("%*f", this->outputWidth, outputData.tea_eps);
	}
	printf("\n");
}

/*
 * Print ldots to show that there are more trajectories running
 * (in case the number of trajectories is too large to show them all).
 */
void OutputManager::printDataLdotsToScreen() const{
	int p, id;
	printf("%*s%*s%*s",
			this->outputWidth, "...",
			this->outputWidth, "...",
			this->outputWidth, "...");
	for(p = 0; p < potentialsCount; p++){
		for(id = 0; id < potentials[p]->getEnergiesCount(); id++){
			printf("%*s", this->outputWidth, "...");
		}
	}
	printf("%*s%*s%*s",
			this->outputWidth, "...",
			this->outputWidth, "...",
			this->outputWidth, "...");
	if(integratorTea){
		printf("%*s",
				this->outputWidth, "...");
	}
	printf("\n");
}

/*
 * Save energy values for a given trajectory into the dat file.
 * The legend is:
 * Step
 * Average temperature
 * Energies for potentials (Covalent, Native, Long-Range)
 * Total potential energy
 * Number of native contacts
 * Radius of gyration (zero if not computed)
 * TEA eps (if hydrodynamics is active)
 *
 * Columns in the dat-file correspond to the std output
 */
void OutputManager::printDataToFile(FILE *dat_file, int traj) const{

	fprintf(dat_file, "%12ld\t%8.5f\t", step, this->temperatures[traj]);

	double total = 0.0;
	int p, id;
	for(p = 0; p < potentialsCount; p++){
		for(id = 0; id < potentials[p]->getEnergiesCount(); id++){
			float pot = potentials[p]->getEnergy(traj, id);
			fprintf(dat_file, "%8.5f\t", pot);
			total += pot;
		}
	}
	fprintf(dat_file, "%8.5f\t", total);
	fprintf(dat_file, "%8d\t", this->nativeCounts[traj]);
	fprintf(dat_file, "%8.5f", this->rgs[traj]);

	if(integratorTea){
		fprintf(dat_file, "\t%2.9f", outputData.tea_eps);
	}
	fprintf(dat_file, "\n");
}
