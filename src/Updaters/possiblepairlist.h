/*
 * possiblepairlist.cuh
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#pragma once

#define POSSIBLEPAIRS_CUTOFF_STRING			"pairs_threshold"
#define POSSIBLEPAIRS_BLOCK_SIZE_STRING		"block_size_possiblepairs"

#define DEFAULT_POSSIBLEPAIRS_CUTOFF			200.0f

class PossiblepairList : public SOPUpdater{
public:
    PossiblepairList();
    virtual ~PossiblepairList();
    virtual void update();
private:
	float pairsThreshold;
	int blockSize;
	int blockNum;
};

void createPossiblepairlistUpdater();
