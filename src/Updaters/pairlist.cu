/*
 * pairlist.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "pairlist.h"
//#define DEBUGPAIRLIST

#include "pairlist_kernel.cu"

void createPairlistUpdater(){
	updaters[updatersCount] = new PairList();
	updatersCount++;
}

/*
 * Building map of possible long-range contacts
 * (i.e. all pairs excluding native and covalent)
 * (for map description see initCovalent())
 */
PairList::PairList(){
	this->name = "Pairlist";
	this->frequency = getIntegerParameter("pairs_freq", 1000, 1);
	printf("Searching for all possible pairs for repulsive LJ...\n");
	this->blockSize = getIntegerParameter(PAIRLIST_BLOCK_SIZE_STRING, gsop.blockSize, 1);
	this->blockNum = gsop.aminoCount/this->blockSize + 1;
	this->max_possiblePairs = getIntegerParameter(MAX_POSSIBLEPAIRS_STRING, DEFAULT_MAX_POSSIBLEPAIRS, 1);
	this->pairlistCutoff = getFloatParameter(PAIRLIST_CUTOFF_STRING, DEFAULT_PAIRLIST_CUTOFF, 1);
	// Allocating memory
	this->h_possiblePairs = (int*)calloc(gsop.aminoCount*this->max_possiblePairs, sizeof(int));
	this->h_possiblePairsCount = (int*)calloc(gsop.aminoCount, sizeof(int));
	cudaMalloc((void**)&this->d_possiblePairs, gsop.aminoCount*this->max_possiblePairs*sizeof(int));
	cudaMalloc((void**)&this->d_possiblePairsCount, gsop.aminoCount*sizeof(int));
	checkCUDAError();
	// Building a map
	int totalPairs = 0;
	int i, j, k;
	for(k = 0; k < sop.pairCount; k++){
		i = sop.pairs[k].i;
		j = sop.pairs[k].j;

		this->h_possiblePairs[this->h_possiblePairsCount[i]*gsop.aminoCount + i] = j;
		this->h_possiblePairs[this->h_possiblePairsCount[j]*gsop.aminoCount + j] = i;
		totalPairs++;

		this->h_possiblePairsCount[i] ++;
		this->h_possiblePairsCount[j] ++;
		if(this->h_possiblePairsCount[i] > this->max_possiblePairs ||
				this->h_possiblePairsCount[j] > this->max_possiblePairs){
			printf("ERROR: Maximum number of possible pairs exceeded the limit of %d.\n", this->max_possiblePairs);
			exit(-1);
		}
	}

	//Duplicating data for mass-production (if gsop.Ntr > 1)
	int traj;
	for(traj = 1; traj < gsop.Ntr; traj++){
		for(i = 0; i < sop.aminoCount; i++){
			for(k = 0; k < this->max_possiblePairs; k++){
				this->h_possiblePairs[traj*sop.aminoCount + i + k*gsop.aminoCount] =
						this->h_possiblePairs[i + k*gsop.aminoCount] + traj*sop.aminoCount;
			}
			this->h_possiblePairsCount[traj*sop.aminoCount + i] =
					this->h_possiblePairsCount[i];
		}
	}

	#ifdef DEBUGPAIRLIST
	printf("Possible pairs:\n");
	for(i = 0; i < gsop.aminoCount; i++){
		printf("%d (%d): ", i, this->h_possiblePairsCount[i]);
		for(j = 0; j < this->h_possiblePairsCount[i]; j++){
			printf("%d ", this->h_possiblePairs[j*gsop.aminoCount + i]);
		}
		printf("\n");
	}
	#endif
	// Copying to the Device
	cudaMemcpy(this->d_possiblePairs, this->h_possiblePairs, gsop.aminoCount*this->max_possiblePairs*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_possiblePairsCount, this->h_possiblePairsCount, gsop.aminoCount*sizeof(int), cudaMemcpyHostToDevice);

    this->updateParametersOnGPU();

	printf("Total number of LJ pairs possible: %d\n", totalPairs);

	if(deviceProp.major == 2){ // TODO: >=  2
		cudaFuncSetCacheConfig(generate_pairs, cudaFuncCachePreferL1);
	}
}

PairList::~PairList(){
	free(this->h_possiblePairs);
	free(this->h_possiblePairsCount);
	cudaFree(this->d_possiblePairs);
	cudaFree(this->d_possiblePairsCount);
}

void PairList::update(){
	if(step % this->frequency == 0){
		//printf("Generating new pairlist...");
		generate_pairs<<<this->blockNum, this->blockSize>>>();
		cudaDeviceSynchronize();
		checkCUDAError();
		//printf("done.\n");

#ifdef DEBUGPAIRLIST
		checkCUDAError();
		cudaMemcpy(pairs.h_pairs, pairs.d_pairs, max_pairs*gsop.aminoCount*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(pairs.h_pairsCount, pairs.d_pairsCount, gsop.aminoCount*sizeof(int), cudaMemcpyDeviceToHost);
		printf("LJ pairs (number, #s of beads, zeros):\n");
		int counter = 0;
		for(int i = 0; i < gsop.aminoCount; i++){
		//int i = 66;
			if(i==66 || i == 13)printf("%d: (%d)", i, pairs.h_pairsCount[i]);
			for(int j = 0; j < max_pairs; j++){
				if(pairs.h_pairs[j*gsop.aminoCount + i] != -1){
					counter++;
				}
				if(i==66 || i == 13)printf("%d ", pairs.h_pairs[j*gsop.aminoCount + i]);
			}
			if(i==66 || i == 13)printf("\n");
		}
		printf("Total number of pairs: %d\n", counter);
#endif
	}
}

void PairList::updateParametersOnGPU(){
    hc_pairList.pairlistCutoff = this->pairlistCutoff;
    hc_pairList.d_possiblePairs = this->d_possiblePairs;
    hc_pairList.d_possiblePairsCount = this->d_possiblePairsCount;
	cudaMemcpyToSymbol(c_pairList, &hc_pairList, sizeof(PairListConstant), 0, cudaMemcpyHostToDevice);
    checkCUDAError();
}

