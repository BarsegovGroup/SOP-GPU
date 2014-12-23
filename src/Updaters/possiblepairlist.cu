/*
 * possiblepairlist.cu
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "possiblepairlist.h"

#include "possiblepairlist_kernel.cu"

void createPossiblepairlistUpdater(){
	updaters[updatersCount] = new PossiblepairList();
	updatersCount++;
}

/*
 * Building map of possible long-range contacts
 * (i.e. all pairs excluding native and covalent)
 * (for map description see initCovalent())
 */
PossiblepairList::PossiblepairList(){
	this->name = "Possiblepairlist";
	this->frequency = getIntegerParameter("possiblepairs_freq", 100000, 1);
	printf("Initializing possible pairs list generator...\n");
	this->blockSize = getIntegerParameter(POSSIBLEPAIRS_BLOCK_SIZE_STRING, gsop.blockSize, 1);
	this->blockNum = gsop.aminoCount/this->blockSize + 1;
	this->pairsThreshold = getFloatParameter(POSSIBLEPAIRS_CUTOFF_STRING, DEFAULT_POSSIBLEPAIRS_CUTOFF, 1);

    this->updateParametersOnGPU();

	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(generate_possiblepairs, cudaFuncCachePreferL1);
	}
}

PossiblepairList::~PossiblepairList(){

}

void PossiblepairList::update(){
	if(step % this->frequency == 0){
		//printf("Regenerating the list of possible pairs...");
		generate_possiblepairs<<<this->blockNum, this->blockSize>>>();
		checkCUDAError();
		//printf("done.\n");
	}
}

void PossiblepairList::updateParametersOnGPU(){
    hc_possiblepairList.pairsThreshold = this->pairsThreshold;
	cudaMemcpyToSymbol(c_possiblepairList, &hc_possiblepairList, sizeof(PossiblepairListConstant), 0, cudaMemcpyHostToDevice);
	checkCUDAError();
}

