/*
 * externalForce_kernel.cu
 *
 *  Created on: Jun 20, 2009
 *      Author: zhmurov
 */
#include "gsop.cuh"

__device__ int do_pulling(int d_i){
	// Fixed beads
	int mask = tex1Dfetch(t_masks, d_i);
	if(mask == 1){
		f.x = 0.0f;
		f.y = 0.0f;
		f.z = 0.0f;
	}
	// Pulled beads
	if(mask == 2){
		f.x += c_ext_force.x;
		f.y += c_ext_force.y;
		f.z += c_ext_force.z;
	}
}
