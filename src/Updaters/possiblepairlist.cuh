/*
 * possiblepairlist.cuh
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#ifndef POSSIBLEPAIRSLIST_CUH_
#define POSSIBLEPAIRSLIST_CUH_

#define POSSIBLEPAIRS_CUTOFF_STRING			"pairs_threshold"
#define POSSIBLEPAIRS_BLOCK_SIZE_STRING		"block_size_possiblepairs"

#define DEFAULT_POSSIBLEPAIRS_CUTOFF			200.0f

#include "../gsop.cuh"

struct PossiblepairList{
	float pairsThreshold;
	int blockSize;
	int blockNum;
};

PossiblepairList possiblepairList;
__device__ __constant__ PossiblepairList c_possiblepairList;

SOPUpdater possiblepairlistMaker;

void createPossiblepairlistUpdater();
void initPossiblepairlist();
void deletePossiblepairlist();
inline void generatePossiblepairlist();

#endif /* POSSIBLEPAIRSLIST_CUH_ */
