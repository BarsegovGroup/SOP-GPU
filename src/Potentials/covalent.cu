/*
 * covalent.cu
 *
 *  Created on: Jun 5, 2009
 *      Author: zhmurov
 */

#include "../gsop.cuh"
#include "covalent.cuh"
#include "covalent_kernel.cu"
#include "../IO/configreader.h"

/*
 * Covalent potential. Adds itself into the list of potentials and creates timer.
 */
void createCovalentPotential(){
	sprintf(covalentPotential.name, "Covalent");
	covalentPotential.compute = &computeCovalent;
	covalentPotential.computeEnergy = &computeCovalentEnergy;
	potentials[potentialsCount] = &covalentPotential;
	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(covalent_kernel, cudaFuncCachePreferL1);
		cudaFuncSetCacheConfig(covalentEnergy_kernel, cudaFuncCachePreferL1);
	}
	potentialsCount++;
	initCovalent();
}

/*
 * Construct map of covalent bonds, copy it to the device memory
 */
void initCovalent(){
	printf("Building a map of covalent bounds...\n");
	// Reading parameters from config file
	covalent.kspring_cov = getFloatParameter(COVALENT_KS_STRING, DEFAULT_COVALENT_KS, 1);
	covalent.R_limit = getFloatParameter(COVALENT_R_LIMIT_STRING, DEFAULT_COVALENT_R_LIMIT, 1);
	covalent.R_limit_sq = covalent.R_limit*covalent.R_limit;
	covalent.max_covalent = getIntegerParameter(MAX_COVALENT_STRING, DEFAULT_MAX_COVALENT, 1);
	covalent.blockSize = getIntegerParameter(COVALENT_BLOCK_SIZE_STRING, BLOCK_SIZE, 1);
	covalent.blockNum = gsop.aminoCount/covalent.blockSize + 1;
	// Allocating memory
	covalent.h_bonds = (GCovalentBond*)calloc(gsop.aminoCount*covalent.max_covalent, sizeof(GCovalentBond));
	covalent.h_covalentCount = (int*)calloc(gsop.aminoCount, sizeof(int));
	cudaMalloc((void**)&covalent.d_bonds, gsop.aminoCount*covalent.max_covalent*sizeof(GCovalentBond));
	cudaMalloc((void**)&covalent.d_covalentCount, gsop.aminoCount*sizeof(int));
	checkCUDAError();

	// Constructing map from sop-model
	// TODO: Construct directly from topology file to make it usefull for all-atom
	int totalCovalent = 0;
	int i, j, k;
	printf("Constructing map of covalent bonds.\n");
	for(k = 0; k < sop.bondCount; k++){
		i = sop.bonds[k].i;
		j = sop.bonds[k].j;
		covalent.h_bonds[covalent.h_covalentCount[i]*gsop.aminoCount + i].i2 = j;
		covalent.h_bonds[covalent.h_covalentCount[j]*gsop.aminoCount + j].i2 = i;
		covalent.h_bonds[covalent.h_covalentCount[i]*gsop.aminoCount + i].r0 = sop.bonds[k].r0;
		covalent.h_bonds[covalent.h_covalentCount[j]*gsop.aminoCount + j].r0 = sop.bonds[k].r0;
		covalent.h_covalentCount[i] ++;
		covalent.h_covalentCount[j] ++;
		if(covalent.h_covalentCount[i] > covalent.max_covalent || covalent.h_covalentCount[j] > covalent.max_covalent){
			printf("ERROR: Maximum number of covalent bonds exceeded the limit of %d.\n", covalent.max_covalent);
			exit(-1);
		}
		totalCovalent++;
	}

	// Multiplying data for many trajectories (Many-runs-per-GPU)
	int traj;
	for(traj = 1; traj < gsop.Ntr; traj++){
		for(i = 0; i < sop.aminoCount; i++){
			for(k = 0; k < covalent.max_covalent; k++){
				covalent.h_bonds[traj*sop.aminoCount + i + k*gsop.aminoCount].i2 =
							covalent.h_bonds[i + k*gsop.aminoCount].i2 + traj*sop.aminoCount;
				covalent.h_bonds[traj*sop.aminoCount + i + k*gsop.aminoCount].r0 =
							covalent.h_bonds[i + k*gsop.aminoCount].r0;
			}
			covalent.h_covalentCount[traj*sop.aminoCount + i] = covalent.h_covalentCount[i];
		}
	}

	/*cudaMemcpyToArray(gsop.d_covalent, 0, 0, gsop.h_covalent, size, cudaMemcpyHostToDevice);
	cudaMemcpyFromArray((void*)gsop.h_covalent, gsop.d_covalent, 0, 0, size, cudaMemcpyHostToDevice);*/
	#ifdef DEBUG
	printf("Covalent bounds (number of bounds, #s of beads, zeros):\n");
	for(int i = 0; i < gsop.aminoCount; i++){
		printf("%d(%d): ", i, gsop.h_covalent[i*width + 0]);
		for(int j = 0; j < max_covalent; j++){
			printf("%d(%3.2f) ", gsop.h_covalent[i*width + j + 1], gsop.h_covalentR0[i*max_covalent + j]);
		}
		printf("\n");
	}
	#endif
	//exit(-1);

	// Copying data to device and binding texture
	cudaMemcpy(covalent.d_bonds, covalent.h_bonds, covalent.max_covalent*gsop.aminoCount*sizeof(GCovalentBond), cudaMemcpyHostToDevice);
	cudaMemcpy(covalent.d_covalentCount, covalent.h_covalentCount, gsop.aminoCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_covalent, &covalent, sizeof(Covalent), 0, cudaMemcpyHostToDevice);
	checkCUDAError();

	//printf("error = %s\n", cudaGetErrorString(error));
	printf("Total number of covalent bonds: %d\n", totalCovalent);

}

/*
 * Start/stop timer + execute kernel
 */
inline void computeCovalent(){
	covalent_kernel<<<covalent.blockNum, covalent.blockSize>>>();
	checkCUDAError();
}

/*
 * Compute energy for .dat output
 */
inline void computeCovalentEnergy(){
	covalentEnergy_kernel<<<covalent.blockNum, covalent.blockSize>>>();
	checkCUDAError();
}

