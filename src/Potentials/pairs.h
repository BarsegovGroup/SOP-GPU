/*
 * pairs.cuh
 *
 *  Created on: Mar 11, 2010
 *      Author: zhmurov
 */

#pragma once

#include "../Util/parameters.h"

PARAMETER(pairs_cutoff, float, 20.0f, "?", "...")
PARAMETER(el, float, 1.0f, "?", "...")
PARAMETER(a, float, 3.8f, "?", "...")
PARAMETER(max_pairs, int, 512, "?", "...")

class PairsPotential : public SOPPotential {
public:
    PairsPotential();
    virtual ~PairsPotential() { }
    virtual void compute();
    virtual int getEnergiesCount() const;
    virtual float* computeEnergy(int id);
    virtual float getEnergy(int traj, int id);

private:
    void updateParametersOnGPU();
	float pairsCutoff;
	float pairsCutoff2;
	float el;
	float a;
	float a2;
	float minus6elovera2;
	int max_pairs;

	int* h_pairs;
	int* h_pairsCount;

	int* d_pairs;
	int* d_pairsCount;

	float* h_energies;
	float* d_energies;
	float* energies;

	int blockSize;
	int blockNum;
};

void createPairsPotential();

