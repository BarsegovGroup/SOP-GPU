/*
 * pairlist.cuh
 *
 *  Created on: Mar 11, 2010
 *      Author: zhmurov
 */

#ifndef PAIRLIST_CUH_
#define PAIRLIST_CUH_

#define PAIRLIST_CUTOFF_STRING			"pairlist_cutoff"
#define PAIRLIST_BLOCK_SIZE_STRING		"block_size_pairlist"
#define MAX_POSSIBLEPAIRS_STRING		"max_possiblePairs"

#define DEFAULT_PAIRLIST_CUTOFF			20.0f
#define DEFAULT_MAX_POSSIBLEPAIRS		4096

struct PairList{

	float pairlistCutoff;
	int max_possiblePairs;

	int* h_possiblePairs;
	int* h_possiblePairsCount;

	int* d_possiblePairs;
	int* d_possiblePairsCount;

	int blockSize;
	int blockNum;
};

PairList pairList;
__device__ __constant__ PairList c_pairList;

SOPUpdater pairlistMaker;


void createPairlist();
void initPairlist();
void deletePairlist();
inline void generatePairlist();

#endif /* PAIRLIST_CUH_ */
