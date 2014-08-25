/*
 * possiblepairlist.cu
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "possiblepairlist.h"

PossiblepairList possiblepairList;
__device__ __constant__ PossiblepairList c_possiblepairList;
SOPUpdater possiblepairlistMaker;

#include "possiblepairlist_kernel.cu"

void createPossiblepairlistUpdater(){
	sprintf(possiblepairlistMaker.name, "Possiblepairlist");
	possiblepairlistMaker.update = &generatePossiblepairlist;
	possiblepairlistMaker.destroy = &deletePossiblepairlist;
	possiblepairlistMaker.frequency = getIntegerParameter("possiblepairs_freq", 100000, 1);
	updaters[updatersCount] = &possiblepairlistMaker;
	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(generate_possiblepairs, cudaFuncCachePreferL1);
	}
	updatersCount++;
	initPossiblepairlist();
}

/*
 * Building map of possible long-range contacts
 * (i.e. all pairs excluding native and covalent)
 * (for map description see initCovalent())
 */
void initPossiblepairlist(){
	printf("Initializing possible pairs list generator...\n");
	possiblepairList.blockSize = getIntegerParameter(POSSIBLEPAIRS_BLOCK_SIZE_STRING, gsop.blockSize, 1);
	possiblepairList.blockNum = gsop.aminoCount/possiblepairList.blockSize + 1;
	possiblepairList.pairsThreshold = getFloatParameter(POSSIBLEPAIRS_CUTOFF_STRING, DEFAULT_POSSIBLEPAIRS_CUTOFF, 1);
	cudaMemcpyToSymbol(c_possiblepairList, &possiblepairList, sizeof(PossiblepairList), 0, cudaMemcpyHostToDevice);
	checkCUDAError();
}

void deletePossiblepairlist(){

}

inline void generatePossiblepairlist(){
	if(step % possiblepairlistMaker.frequency == 0){
		printf("Regenerating the list of possible pairs...");
		generate_possiblepairs<<<possiblepairList.blockNum, possiblepairList.blockSize>>>();
		checkCUDAError();
		printf("done.\n");
	}
}

