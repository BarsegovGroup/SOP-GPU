/*
 * possiblepairlist.cuh
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 *
 *      TODO: replace headers everywhere to reference SOP-GPU
 */

#pragma once

#include "../Util/parameters.h"

PARAMETER(pairs_threshold, float, 200.0f, "?", "...")
PARAMETER(possiblepairs_freq, int, 100000, "?", "...")

class PossiblepairList : public SOPUpdater{
public:
    PossiblepairList();
    virtual ~PossiblepairList();
    virtual void update();
private:
    void updateParametersOnGPU();
	float pairsThreshold;
	int blockSize;
	int blockNum;
};

void createPossiblepairlistUpdater();

