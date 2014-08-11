/*
 * pairs.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "pairs.cuh"
#include "pairs_kernel.cu"

void createPairsPotential(){
	sprintf(pairsPotential.name, "Long Range");
	pairsPotential.compute = &computePairs;
	pairsPotential.computeEnergy = &computePairsEnergy;
	potentials[potentialsCount] = &pairsPotential;
	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(pairs_kernel, cudaFuncCachePreferL1);
		cudaFuncSetCacheConfig(pairsEnergy_kernel, cudaFuncCachePreferL1);
	}
	potentialsCount++;
	initPairs();
}
/*
 * Constructing initial pairlist map
 * (see description of map for covalent bonds)
 */
void initPairs(){
	printf("Initializing pairlist...");
	pairs.blockSize = getIntegerParameter(PAIRS_BLOCK_SIZE_STRING, BLOCK_SIZE, 1);
	pairs.blockNum = gsop.aminoCount/pairs.blockSize + 1;
	pairs.a = getFloatParameter(PAIRS_A_STRING, DEFAULT_PAIRS_A, 1);
	pairs.a2 = pairs.a*pairs.a;
	pairs.el = getFloatParameter(PAIRS_EL_STRING, DEFAULT_PAIRS_EL, 1);
	pairs.minus6elovera2 = -6.0f*pairs.el/pairs.a2;
	pairs.max_pairs = getIntegerParameter(MAX_PAIRS_STRING, DEFAULT_MAX_PAIRS, 1);
	pairs.pairsCutoff = getFloatParameter(PAIRS_CUTOFF_STRING, DEFAULT_PAIRS_CUTOFF, 1);
	pairs.pairsCutoff2 = pairs.pairsCutoff*pairs.pairsCutoff;
	// Allocating memory
	unsigned int width = gsop.aminoCount;
	unsigned int height = pairs.max_pairs;
	unsigned int size = width*height*sizeof(int);
	pairs.h_pairs = (int*)calloc(width*height, sizeof(int));
	pairs.h_pairsCount = (int*)calloc(gsop.aminoCount, sizeof(int));
	// Allocating memory on the Device
	cudaMalloc((void**)&pairs.d_pairs, size);
	checkCUDAError();
	cudaMalloc((void**)&pairs.d_pairsCount, gsop.aminoCount*sizeof(int));
	checkCUDAError();
	cudaMemcpyToSymbol(c_pairs, &pairs, sizeof(Pairs), 0, cudaMemcpyHostToDevice);
	checkCUDAError();
	printf("done.\n");
}

inline void computePairs(){
	pairs_kernel<<<pairs.blockNum, pairs.blockSize>>>();
	checkCUDAError();
}

inline void computePairsEnergy(){
	pairsEnergy_kernel<<<pairs.blockNum, pairs.blockSize>>>();
	checkCUDAError();
}
