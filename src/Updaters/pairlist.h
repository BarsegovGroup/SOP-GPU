/*
 * pairlist.cuh
 *
 *  Created on: Mar 11, 2010
 *      Author: zhmurov
 */

#pragma once

#define PAIRLIST_CUTOFF_STRING			"pairlist_cutoff"
#define PAIRLIST_BLOCK_SIZE_STRING		"block_size_pairlist"
#define MAX_POSSIBLEPAIRS_STRING		"max_possiblePairs"

#define DEFAULT_PAIRLIST_CUTOFF			20.0f
#define DEFAULT_MAX_POSSIBLEPAIRS		4096

class PairList : public SOPUpdater{
public:
    PairList();
    virtual ~PairList();
    virtual void update();
private:
	float pairlistCutoff;
	int max_possiblePairs;

	int* h_possiblePairs;
	int* h_possiblePairsCount;

	int* d_possiblePairs;
	int* d_possiblePairsCount;

	int blockSize;
	int blockNum;
};

void createPairlist();

