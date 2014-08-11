/*
 * pairs.cuh
 *
 *  Created on: Mar 11, 2010
 *      Author: zhmurov
 */

#ifndef PAIRS_CUH_
#define PAIRS_CUH_

#define PAIRS_CUTOFF_STRING			"pairs_cutoff"
#define PAIRS_EL_STRING				"el"
#define PAIRS_A_STRING				"a"
#define MAX_PAIRS_STRING			"max_pairs"
#define PAIRS_BLOCK_SIZE_STRING		"block_size_pairs"

#define DEFAULT_PAIRS_CUTOFF			20.0f
#define DEFAULT_PAIRS_EL				1.0f
#define DEFAULT_PAIRS_A					3.8f
#define DEFAULT_MAX_PAIRS				512

struct Pairs{

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

	int blockSize;
	int blockNum;
};

Pairs pairs;
__device__ __constant__ Pairs c_pairs;
SOPPotential pairsPotential;

void createPairsPotential();
void initPairs();
inline void computePairs();
inline void computePairsEnergy();

#endif /* PAIRS_CUH_ */