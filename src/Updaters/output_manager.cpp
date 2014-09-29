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
	printf("Initializing output...\n");
	int traj;
    std::string dat_filename = getMaskedParameterAs<std::string>(OUTPUT_FILENAME_STRING, DEFAULT_OUTPUT_FILENAME);
	dat_filenames.resize(gsop.Ntr);
	for(traj = 0; traj < gsop.Ntr; traj++){
        dat_filenames[traj] = string_replace(dat_filename, "<run>", traj+gsop.firstrun);
		dat_file = fopen(dat_filenames[traj].c_str(), "w");
		fclose(dat_file);
	}
	printf("Output will be saved in '%s'.\n", dat_filename.c_str());

    char modeString[100];
	getParameter(modeString, "mode", "", 1);
	if(strcmp(modeString, "capsid") == 0){
		mode = MODE_CAPSID;
	} else {
		mode = 0;
	}

    this->covalent = (CovalentPotential*)potentialByName("Covalent");
    this->native = (NativePotential*)potentialByName("Native");
    this->tea = (TeaIntegrator*)integrator;

	printf("Done initializing output manager...\n");
}

OutputManager::~OutputManager(){
}

void OutputManager::update(){
	if(step % this->frequency == 0){
		copyCoordDeviceToHost();
		int p;
		for(p = 0; p < potentialsCount; p++){
			potentials[p]->computeEnergy();
		}
		int size = gsop.aminoCount*sizeof(float4);
        copyEnergiesDeviceToHost();
		char runstring[20];
		if(gsop.Ntr == 1){
			sprintf(runstring, "Run %d", gsop.firstrun);
		} else {
			sprintf(runstring, "%d runs", gsop.Ntr);
		}
		printf("Writing output at step %ld of %ld. %s on %d [%s].\n", step, numsteps, runstring, gsop.deviceId, deviceProp.name);
		int traj;
		for(traj = 0; traj < gsop.Ntr; traj++){
			this->outputData.step = step;
			this->computeEnergies(traj);
			this->computeNativeNumber(traj);
			if(mode != MODE_CAPSID){
				this->computeRg(traj);
			}
			this->computeTEAeps(traj);

			if(traj == 0){
				this->printDataToScreen();
			}
			FILE* dat_file = fopen(dat_filenames[traj].c_str(), "a");
			this->printDataToFile(dat_file);
			fclose(dat_file);
		}
		printf("Done writing output.\n");
	}
    if (step % updaterByName("Pairlist")->frequency == 0) {
        //TODO: resetting should be done in OutputManager::computeEnergies
        //It is here only to provide complete compatibility with old SOP
        this->resetTemperatureCounter();
    }
}

void OutputManager::computeEnergies(int traj){
	int i;
	outputData.tempav = 0.0; // Average temperature
	outputData.epot_LJ = 0.0; // Total LJ energy
	outputData.epot_fene = 0.0; // Total FENE potential energy
	outputData.epot_tot = 0.0; // Total energy
	outputData.epot_native = 0.0; // Total energy of native interactions
	outputData.epot_longrange = 0.0; // Total energy for pairs interactions

	// Computing a sum of energies and temperature
	for(i = 0; i < sop.aminoCount; i++){
		if(isinf(gsop.h_energies[traj*sop.aminoCount + i].x)){
			printf("WARNING: Bead #%d has infinite FENE energy\n", i);
		}
		outputData.epot_fene += gsop.h_energies[traj*sop.aminoCount + i].x;
		if(isinf(gsop.h_energies[traj*sop.aminoCount + i].y)){
			printf("WARNING: Bead #%d has infinite native interaction energy\n", i);
		}
		outputData.epot_native += gsop.h_energies[traj*sop.aminoCount + i].y;
		if(isinf(gsop.h_energies[traj*sop.aminoCount + i].z)){
			printf("WARNING: Bead #%d has infinite long-range LJ energy\n", i);
		}
		outputData.epot_longrange += gsop.h_energies[traj*sop.aminoCount + i].z;
		outputData.tempav += gsop.h_energies[traj*sop.aminoCount + i].w;
	}
	outputData.epot_fene /= 2.0f;
	outputData.epot_native /= 2.0f;
	outputData.epot_longrange /= 2.0f;

	// Normalizing the temperature
	outputData.tempav /= ((double)(sop.aminoCount*updaterByName("Pairlist")->frequency));
    //TODO: we should use following code:
	//outputData.tempav /= ((double)(sop.aminoCount*this->frequency));
    //this->resetTemperatureCounter();

	// Other energies
	outputData.epot_LJ = outputData.epot_native + outputData.epot_longrange;
	outputData.epot_tot = outputData.epot_LJ + outputData.epot_fene;
}

void OutputManager::computeNativeNumber(int traj){
	int i,j,k;
	outputData.nat_num = 0;
	for(k = 0; k < sop.nativeCount; k++){
		i = sop.natives[k].i;
		j = sop.natives[k].j;
		float dr;
		float3 r;
		r.x = gsop.h_coord[traj*sop.aminoCount + i].x - gsop.h_coord[traj*sop.aminoCount + j].x;
		r.y = gsop.h_coord[traj*sop.aminoCount + i].y - gsop.h_coord[traj*sop.aminoCount + j].y;
		r.z = gsop.h_coord[traj*sop.aminoCount + i].z - gsop.h_coord[traj*sop.aminoCount + j].z;
		dr = sqrtf(r.x*r.x+r.y*r.y+r.z*r.z);
		if(dr < native->R_limit_bond || fabs(sop.natives[k].r0 - dr) <= covalent->R_limit){
			outputData.nat_num ++;
		}
	}
}

void OutputManager::computeRg(int traj){
	int i, j;
	double rgsq;
	rgsq = 0.0;
	int shift = traj*sop.aminoCount;
	for(i = 0 ; i < (sop.aminoCount-1) ; i ++){
		for(j = (i+1) ; j < sop.aminoCount ; j ++){
			rgsq = rgsq + (gsop.h_coord[shift + i].x-gsop.h_coord[shift + j].x)*(gsop.h_coord[shift + i].x-gsop.h_coord[shift + j].x)+
				  (gsop.h_coord[shift + i].y-gsop.h_coord[shift + j].y)*(gsop.h_coord[shift + i].y-gsop.h_coord[shift + j].y)+
				  (gsop.h_coord[shift + i].z-gsop.h_coord[shift + j].z)*(gsop.h_coord[shift + i].z-gsop.h_coord[shift + j].z);
		}
	}

	rgsq = rgsq/(double)(sop.aminoCount*sop.aminoCount);
	outputData.rg = sqrt(rgsq);
}

void OutputManager::computeTEAeps(int traj){
	if (integratorTea && !this->tea->exact){
		outputData.tea_eps = this->tea->h_epsilon[traj];
	}else{
		outputData.tea_eps = 0./0.;
	}
}

void OutputManager::printDataToScreen() const{
	printf("TimeStep\tTemp    \tPotent  \tNative  \tLonRan  \tLJ      \tFENE    \tNative# \tRg");
	if (integratorTea) 
		printf("	\tTEA-eps");
	printf("\n");
	printf("%12ld\t%8.3f\t"
					"%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t"
					"%3d\t%8.3f",
					outputData.step, outputData.tempav,
					outputData.epot_tot, outputData.epot_native,
					outputData.epot_longrange, outputData.epot_LJ, outputData.epot_fene,
					outputData.nat_num, outputData.rg);
	if (integratorTea) 
		printf("\t%2.9f", outputData.tea_eps);
	printf("\n");
}

void OutputManager::printDataToFile(FILE *dat_file) const{
	fprintf(dat_file, "%12ld\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t"
					"%8d\t%8.5f",
					step, outputData.tempav,
					outputData.epot_tot, outputData.epot_native, outputData.epot_longrange,
					outputData.epot_LJ, outputData.epot_fene,
					outputData.nat_num, outputData.rg);
	if (integratorTea) 
		fprintf(dat_file, "\t%2.9f", outputData.tea_eps);
	fprintf(dat_file, "\n");
}

