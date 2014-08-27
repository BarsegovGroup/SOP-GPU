/*
 * pairs.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "pairs.h"


#include "pairs_kernel.cu"

void createPairsPotential(){
	potentials[potentialsCount] = new PairsPotential();
	potentialsCount++;
}
/*
 * Constructing initial pairlist map
 * (see description of map for covalent bonds)
 */
PairsPotential::PairsPotential(){
    this->name = "Long Range";
	printf("Initializing pairlist...");
	this->blockSize = getIntegerParameter(PAIRS_BLOCK_SIZE_STRING, gsop.blockSize, 1);
	this->blockNum = gsop.aminoCount/this->blockSize + 1;
	this->a = getFloatParameter(PAIRS_A_STRING, DEFAULT_PAIRS_A, 1);
	this->a2 = this->a*this->a;
	this->el = getFloatParameter(PAIRS_EL_STRING, DEFAULT_PAIRS_EL, 1);
	this->minus6elovera2 = -6.0f*this->el/this->a2;
	this->max_pairs = getIntegerParameter(MAX_PAIRS_STRING, DEFAULT_MAX_PAIRS, 1);
	this->pairsCutoff = getFloatParameter(PAIRS_CUTOFF_STRING, DEFAULT_PAIRS_CUTOFF, 1);
	this->pairsCutoff2 = this->pairsCutoff*this->pairsCutoff;
	// Allocating memory
	unsigned int width = gsop.aminoCount;
	unsigned int height = this->max_pairs;
	unsigned int size = width*height*sizeof(int);
	this->h_pairs = (int*)calloc(width*height, sizeof(int));
	this->h_pairsCount = (int*)calloc(gsop.aminoCount, sizeof(int));
	// Allocating memory on the Device
	cudaMalloc((void**)&this->d_pairs, size);
	checkCUDAError();
	cudaMalloc((void**)&this->d_pairsCount, gsop.aminoCount*sizeof(int));
    checkCUDAError();

    this->updateParametersOnGPU();
	printf("done.\n");

	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(pairs_kernel, cudaFuncCachePreferL1);
		cudaFuncSetCacheConfig(pairsEnergy_kernel, cudaFuncCachePreferL1);
	}
}

void PairsPotential::updateParametersOnGPU(){
    hc_pairs.pairsCutoff2 = this->pairsCutoff2;
    hc_pairs.a2 = this->a2;
    hc_pairs.el = this->el;
    hc_pairs.minus6elovera2 = this->minus6elovera2;
    hc_pairs.d_pairs = this->d_pairs;
    hc_pairs.d_pairsCount = this->d_pairsCount;
    hc_pairs.max_pairs = this->max_pairs;
	cudaMemcpyToSymbol(c_pairs, &hc_pairs, sizeof(PairsConstant), 0, cudaMemcpyHostToDevice);
	checkCUDAError();
}

void PairsPotential::compute(){
	pairs_kernel<<<this->blockNum, this->blockSize>>>();
	checkCUDAError();
}

void PairsPotential::computeEnergy(){
	pairsEnergy_kernel<<<this->blockNum, this->blockSize>>>();
	checkCUDAError();
}

