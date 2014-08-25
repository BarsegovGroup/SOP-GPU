/*
 * possiblepairlist_kernel.cu
 *
 *  Created on: Apr 9, 2010
 *      Author: zhmurov
 */

#include "../gsop.cuh"

struct PossiblepairListConstant{
    float pairsThreshold;
};

PossiblepairListConstant hc_possiblepairList;
__device__ __constant__ PossiblepairListConstant c_possiblepairList;

__global__ void generate_possiblepairs(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	int counter;
	if(d_i < c_gsop.aminoCount){
		// Get coordinate of current amino acid
#ifdef NOTEXTURE
		float4 r1 = c_gsop.d_coord[d_i];
#else
		float4 r1 = tex1Dfetch(t_coord, d_i);
#endif
		// Number of pairs counter
		counter = 0;
		// Iterate over all possible pairs from t_possiblepairs map
		for(int i = 0; i < c_gsop.aminoCount; i++){
			if(i != d_i){
#ifdef NOTEXTURE
				float4 r2 = c_gsop.d_coord[i];
#else
				float4 r2 = tex1Dfetch(t_coord, i);
#endif
				// Compute the distance
				r2.x -= r1.x;
				r2.y -= r1.y;
				r2.z -= r1.z;
				r2.w = sqrtf(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z);
				// If distance < threshold - the pair is qualified and added to the map
				if(r2.w < c_possiblepairList.pairsThreshold){
					int add = 1;
					for(int j = 0; j < c_covalent.d_covalentCount[d_i]; j++){
						if(i == c_covalent.d_bonds[j*c_gsop.aminoCount + d_i].i2){
							add = 0;
						}
					}
					if(add == 1){
						for(int j = 0; j < c_native.d_nativeCount[d_i]; j++){
							if(i == c_native.d_native[j*c_gsop.aminoCount + d_i]){
								add = 0;
							}
						}
					}
					//counter ++;
					if(add == 1){
						c_pairList.d_possiblePairs[counter*c_gsop.aminoCount + d_i] = i;
						counter ++;
					}
				}
			}
		}
		// Write number of pairs added
		c_pairList.d_possiblePairsCount[d_i] = counter;
	}
}

