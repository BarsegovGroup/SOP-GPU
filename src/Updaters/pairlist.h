/*
 * pairlist.cuh
 *
 *  Created on: Mar 11, 2010
 *      Author: zhmurov
 */

#pragma once

#include "../Util/parameters.h"

PARAMETER(pairlist_cutoff, float, 20.0f, "?", "...")
PARAMETER(pairs_freq, int, 1000, "steps", "...")
PARAMETER(max_possiblePairs, int, 4096, "?", "...")

class PairList : public SOPUpdater{
public:
    PairList();
    virtual ~PairList();
    virtual void update();
private:
    void updateParametersOnGPU();
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

