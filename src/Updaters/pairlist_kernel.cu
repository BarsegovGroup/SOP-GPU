/*
 * pairlist.cu
 *
 *  Created on: Oct 27, 2008
 *      Author: zhmurov
 */
#include "../gsop.cuh"

/*
 *  Build a pairlist for long range interactions
 */


/*
 * CUDA Kernel.
 * d_pairs - place for writing output
 * N = sop.aminoCount
 */
__global__ void generate_pairs(){
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
		for(int i = 0; i < c_pairList.d_possiblePairsCount[d_i]; i++){
			// Get index of second amino acid
			int i2 = c_pairList.d_possiblePairs[i*c_gsop.aminoCount + d_i];
			// Get coordinate of second amino acid
#ifdef NOTEXTURE
		float4 r2 = c_gsop.d_coord[i2];
#else
		float4 r2 = tex1Dfetch(t_coord, i2);
#endif
			// Compute the distance
  			r2.x -= r1.x;
			r2.y -= r1.y;
			r2.z -= r1.z;
			r2.w = sqrtf(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z);
			// If distance < cutoff - the pair is qualified and added to the map
			if(r2.w < c_pairList.pairlistCutoff){
				//counter ++;
				c_pairs.d_pairs[counter*c_gsop.aminoCount + d_i] = i2;
				counter ++;
			}

		}
		for(int i = counter; i < c_pairs.max_pairs; i++){
			c_pairs.d_pairs[i*c_gsop.aminoCount + d_i] = -1;
		}
		// Write number of pairs added
		c_pairs.d_pairsCount[d_i] = counter;
		c_gsop.d_energies[d_i].w = 0.0f;
	}
}
