/*
 * pairlist.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "pairlist.cuh"
#include "pairlist_kernel.cu"
//#define DEBUGPAIRLIST

void createPairlistUpdater(){
	sprintf(pairlistMaker.name, "Pairlist");
	pairlistMaker.update = &generatePairlist;
	pairlistMaker.destroy = &deletePairlist;
	pairlistMaker.frequency = getIntegerParameter("pairs_freq", 1000, 1);
	updaters[updatersCount] = &pairlistMaker;
	if(deviceProp.major == 2){ // TODO: >=  2
		cudaFuncSetCacheConfig(generate_pairs, cudaFuncCachePreferL1);
	}
	updatersCount++;
	initPairlist();
}

/*
 * Building map of possible long-range contacts
 * (i.e. all pairs excluding native and covalent)
 * (for map description see initCovalent())
 */
void initPairlist(){
	printf("Searching for all possible pairs for repulsive LJ...\n");
	pairList.blockSize = getIntegerParameter(PAIRLIST_BLOCK_SIZE_STRING, BLOCK_SIZE, 1);
	pairList.blockNum = gsop.aminoCount/pairList.blockSize + 1;
	pairList.max_possiblePairs = getIntegerParameter(MAX_POSSIBLEPAIRS_STRING, DEFAULT_MAX_POSSIBLEPAIRS, 1);
	pairList.pairlistCutoff = getFloatParameter(PAIRLIST_CUTOFF_STRING, DEFAULT_PAIRLIST_CUTOFF, 1);
	// Allocating memory
	pairList.h_possiblePairs = (int*)calloc(gsop.aminoCount*pairList.max_possiblePairs, sizeof(int));
	pairList.h_possiblePairsCount = (int*)calloc(gsop.aminoCount, sizeof(int));
	cudaMalloc((void**)&pairList.d_possiblePairs, gsop.aminoCount*pairList.max_possiblePairs*sizeof(int));
	cudaMalloc((void**)&pairList.d_possiblePairsCount, gsop.aminoCount*sizeof(int));
	checkCUDAError();
	// Building a map
	int totalPairs = 0;
	int i, j, k;
	for(k = 0; k < sop.pairCount; k++){
		i = sop.pairs[k].i;
		j = sop.pairs[k].j;

		pairList.h_possiblePairs[pairList.h_possiblePairsCount[i]*gsop.aminoCount + i] = j;
		pairList.h_possiblePairs[pairList.h_possiblePairsCount[j]*gsop.aminoCount + j] = i;
		totalPairs++;

		pairList.h_possiblePairsCount[i] ++;
		pairList.h_possiblePairsCount[j] ++;
		if(pairList.h_possiblePairsCount[i] > pairList.max_possiblePairs ||
				pairList.h_possiblePairsCount[j] > pairList.max_possiblePairs){
			printf("ERROR: Maximum number of possible pairs exceeded the limit of %d.\n", pairList.max_possiblePairs);
			exit(-1);
		}
	}

	//Duplicating data for mass-production (if gsop.Ntr > 1)
	int traj;
	for(traj = 1; traj < gsop.Ntr; traj++){
		for(i = 0; i < sop.aminoCount; i++){
			for(k = 0; k < pairList.max_possiblePairs; k++){
				pairList.h_possiblePairs[traj*sop.aminoCount + i + k*gsop.aminoCount] =
						pairList.h_possiblePairs[i + k*gsop.aminoCount] + traj*sop.aminoCount;
			}
			pairList.h_possiblePairsCount[traj*sop.aminoCount + i] =
					pairList.h_possiblePairsCount[i];
		}
	}

	#ifdef DEBUGPAIRLIST
	printf("Possible pairs:\n");
	for(i = 0; i < gsop.aminoCount; i++){
		printf("%d (%d): ", i, pairList.h_possiblePairsCount[i]);
		for(j = 0; j < pairList.h_possiblePairsCount[i]; j++){
			printf("%d ", pairList.h_possiblePairs[j*gsop.aminoCount + i]);
		}
		printf("\n");
	}
	#endif
	// Copying to the Device
	cudaMemcpy(pairList.d_possiblePairs, pairList.h_possiblePairs, gsop.aminoCount*pairList.max_possiblePairs*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(pairList.d_possiblePairsCount, pairList.h_possiblePairsCount, gsop.aminoCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_pairList, &pairList, sizeof(PairList), 0, cudaMemcpyHostToDevice);

	checkCUDAError();
	printf("Total number of LJ pairs possible: %d\n", totalPairs);
}

void deletePairlist(){
	free(pairList.h_possiblePairs);
	free(pairList.h_possiblePairsCount);
	cudaFree(pairList.d_possiblePairs);
	cudaFree(pairList.d_possiblePairsCount);

}

inline void generatePairlist(){
	if(step % pairlistMaker.frequency == 0){
		//printf("Generating new pairlist...");
		generate_pairs<<<pairList.blockNum, pairList.blockSize>>>();
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
