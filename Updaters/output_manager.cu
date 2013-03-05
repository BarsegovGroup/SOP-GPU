/*
 * output_manager.cu
 *
 *  Created on: Apr 7, 2010
 *      Author: zhmurov
 */

#include "../def_param.h"
#include "../gsop.cuh"
#include <string.h>
#include <stdio.h>
#include "output_manager.cuh"

FILE* dat_file;
char** dat_filenames;
char tempFilename[100];

inline void computeEnergies(int traj);
inline void computeNativeNumber(int traj);
inline void computeRg(int traj);

inline void printDataToScreen();
inline void printDataToFile(int traj);

extern void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace);

void createOutputManager(){
	printf("Initializing output manager...\n");
	sprintf(outputManager.name, "DAT output");
	outputManager.update = &printStep;
	outputManager.destroy = &closeDAT;
	outputManager.frequency = getIntegerParameter(OUTPUT_FREQUENCY_STRING, DEFAULT_OUTPUT_FREQUENCY, 1);
	updaters[updatersCount] = &outputManager;
	updatersCount++;
	initOutputManager();
	printf("Done initializing output manager...\n");
}

void initOutputManager(){
	printf("Initializing output...\n");
	int traj;
	/*dat_file = (FILE**)calloc(Ntr, sizeof(FILE*));
	for(traj = 0; traj < Ntr; traj++){
		char trajnum[10];
		sprintf(trajnum, "%d\0", traj+firstrun);
		replaceString(tempFilename, dat_filename, trajnum, "<run>");
		dat_file[traj] = fopen(tempFilename, "w");
	}*/
	dat_filenames = (char**)calloc(Ntr, sizeof(char*));
	char dat_filename[100];
	getMaskedParameter(dat_filename, OUTPUT_FILENAME_STRING, DEFAULT_OUTPUT_FILENAME, 1);;
	for(traj = 0; traj < Ntr; traj++){
		dat_filenames[traj] = (char*)calloc(100, sizeof(char));
		char trajnum[10];
		sprintf(trajnum, "%d\0", traj+firstrun);
		replaceString(dat_filenames[traj], dat_filename, trajnum, "<run>");
		dat_file = fopen(dat_filenames[traj], "w");
		fclose(dat_file);
	}
	printf("Output will be saved in '%s'.\n", dat_filename);
}

void closeDAT(){
	/*int traj;
	for(traj = 0; traj < Ntr; traj++){
		fclose(dat_file[traj]);
	}*/
}

inline void printStep(){
	if(step % outputManager.frequency == 0){
		copyCoordDeviceToHost();
		int p;
		for(p = 0; p < potentialsCount; p++){
			potentials[p]->computeEnergy();
		}
		int size = gsop.aminoCount*sizeof(float4);
		cudaMemcpy(gsop.h_energies, gsop.d_energies, size, cudaMemcpyDeviceToHost);
		//cudaMemcpy(gsop.d_coordToSave, gsop.d_coord, size, cudaMemcpyDeviceToDevice);
		//cudaMemcpy(gsop.d_energiesToSave, gsop.d_energies, size, cudaMemcpyDeviceToDevice);
		char runstring[20];
		if(Ntr == 1){
			sprintf(runstring, "Run %d", run);
		} else {
			sprintf(runstring, "%d runs", Ntr);
		}
		printf("Writing output at step %ld of %ld. %s on %d.\n", step, numsteps, runstring, device);
		int traj;
		for(traj = 0; traj < Ntr; traj++){
			dat_file = fopen(dat_filenames[traj], "a");
			outputData.step = step;
			//outputData.temp = temp;
			computeEnergies(traj);
			computeNativeNumber(traj);
			if(mode != MODE_CAPSID){
				computeRg(traj);
			}
			if(traj == 0){
				printDataToScreen();
			}
			printDataToFile(traj);
			fflush(dat_file);
			fclose(dat_file);
		}
		//cudaMemcpyAsync(gsop.h_coord, gsop.d_coordToSave, size, cudaMemcpyDeviceToHost, 0);
		//cudaMemcpyAsync(gsop.h_energies, gsop.d_energiesToSave, size, cudaMemcpyDeviceToHost, 0);
		printf("Done writing output.\n");
	}
}

inline void computeEnergies(int traj){
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
	outputData.tempav /= ((double)(sop.aminoCount*pairlistMaker.frequency));

	// Other energies
	outputData.epot_LJ = outputData.epot_native + outputData.epot_longrange;
	outputData.epot_tot = outputData.epot_LJ + outputData.epot_fene;
}

inline void computeNativeNumber(int traj){
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
		if(dr < native.R_limit_bond || fabs(sop.natives[k].r0 - dr) <= covalent.R_limit){
			outputData.nat_num ++;
		}
	}
}

inline void computeRg(int traj){
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

/*inline void computeEndToEnd(int traj){
	int shift = traj*sop.aminoCount;
	outputData.endToEnd_x = (gsop.h_coord[shift + pulled_end].x - gsop.h_coord[shift + fixed_end].x)*pull_vector[traj].x +
			(gsop.h_coord[shift + pulled_end].y - gsop.h_coord[shift + fixed_end].y)*pull_vector[traj].y +
			(gsop.h_coord[shift + pulled_end].z - gsop.h_coord[shift + fixed_end].z)*pull_vector[traj].z;
	outputData.endToEnd = sqrtf((gsop.h_coord[shift + pulled_end].x-gsop.h_coord[shift + fixed_end].x)*
			(gsop.h_coord[shift + pulled_end].x-gsop.h_coord[shift + fixed_end].x) +
			(gsop.h_coord[shift + pulled_end].y-gsop.h_coord[shift + fixed_end].y)*
			(gsop.h_coord[shift + pulled_end].y-gsop.h_coord[shift + fixed_end].y) +
			(gsop.h_coord[shift + pulled_end].z-gsop.h_coord[shift + fixed_end].z)*
			(gsop.h_coord[shift + pulled_end].z-gsop.h_coord[shift + fixed_end].z));
}*/

/*inline void computePullForce(int traj){
	outputData.f = ext_force[traj].x*pull_vector[traj].x + ext_force[traj].y*pull_vector[traj].y + ext_force[traj].z*pull_vector[traj].z;
	outputData.fx = ext_force[traj].x;
	outputData.fy = ext_force[traj].y;
	outputData.fz = ext_force[traj].z;
	outputData.xt = xt;
}*/


/*inline void computeCapsidR(int traj){
	int i, shift;
	shift = traj*sop.aminoCount;
	for(i = 0; i < sop.aminoCount; i++){
		outputData.R += sqrtf(gsop.h_coord[shift + i].x*gsop.h_coord[shift + i].x +
				gsop.h_coord[shift + i].y*gsop.h_coord[shift + i].y +
				gsop.h_coord[shift + i].z*gsop.h_coord[shift + i].z);
	}
	outputData.R /= sop.aminoCount;

	float surf = 4.0*M_PI*outputData.R*outputData.R;
	outputData.V = surf*outputData.R/3.0;
}*/

inline void computeIndentationForce(int traj){

}

inline void printDataToScreen(){
	/*switch (stage) {
		case pull_stage:
			printf("\tTimeStep\tNative#\tR\tRx\tFext\tfx\tfy\tfz\n");
			printf("%s: %12ld\t%3d\t"
					"%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t\n",
						stageString,
						outputData.step, outputData.nat_num,
						outputData.endToEnd);//, outputData.endToEnd_x,
						//outputData.f, outputData.fx, outputData.fy, outputData.fz);
			break;
		case heat_stage:
			printf("\tTimeStep\tNative#\tR\tTemp\n");
			printf("%s: %12ld\t%3d\t"
					"%5.3f\t%5.3f\n",
						stageString,
						outputData.step, outputData.nat_num,
						outputData.endToEnd, outputData.temp);
			break;
		case inden_stage:
			printf("\tTimeStep\tNative#\txt\tFin\tfx\tfy\tfz\n");
			printf("%s: %12ld\t%3d\t"
					"%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t\n",
						stageString,
						outputData.step, outputData.nat_num,
						outputData.xt, outputData.f, outputData.fx, outputData.fy, outputData.fz);
			break;
		default:
			printf("\tTimeStep\tNative#\tR\n");
			printf("%s: %12ld\t%3d\t"
					"%5.3f\n",
						stageString,
						outputData.step, outputData.nat_num,
						outputData.endToEnd);
			break;
	}
	if(mode == MODE_CAPSID){
		printf("\tTimeStep\tNative#\tV\tR\tTemp\n");
		printf("%s: %12ld\t%3d\t"
					"%8.1f\t%5.3f\t%5.3f\n",
						stageString,
						outputData.step, outputData.nat_num,
						outputData.V, outputData.R, outputData.temp);
	}*/
	printf("TimeStep\tTemp    \tPotent  \tNative  \tLonRan  \tLJ      \tFENE    \tNative# \tRg\n");
	printf("%12ld\t%8.3f\t"
					"%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t"
					"%3d\t%8.3f\n",
					outputData.step, outputData.tempav,
					outputData.epot_tot, outputData.epot_native,
					outputData.epot_longrange, outputData.epot_LJ, outputData.epot_fene,
					outputData.nat_num, outputData.rg);
}

inline void printDataToFile(int traj){
	if(mode != MODE_CAPSID){
		fprintf(dat_file, "%12ld\t%8.5f\t"
						"%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t"
						"%8d\t%8.5f\n",
						step, outputData.tempav,
						outputData.epot_tot, outputData.epot_native, outputData.epot_longrange,
						outputData.epot_LJ, outputData.epot_fene,
						outputData.nat_num, outputData.rg);
	} else {
		fprintf(dat_file, "%12ld\t%8.5f\t"
						"%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t"
						"%8d\t%8.5f\n",
						step, outputData.tempav,
						outputData.epot_tot, outputData.epot_native, outputData.epot_longrange,
						outputData.epot_LJ, outputData.epot_fene,
						outputData.nat_num, outputData.rg);
	}
	//fflush(dat_file);
}
