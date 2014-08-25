/*
 * pulling_kernel.cu
 *
 *  Created on: May 25, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"

struct PullingConstant{
	float4* d_extForces;
};

PullingConstant hc_pulling;
__device__ __constant__ PullingConstant c_pulling;

__global__ void pulling_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		float4 f = c_gsop.d_forces[d_i];
		float4 extF = c_pulling.d_extForces[d_i];
		if(extF.w == 1.0f){
			f.x = 0.0f;
			f.y = 0.0f;
			f.z = 0.0f;
		}
		// Pulled beads
		if(extF.w == 2.0f){
			f.x += extF.x;
			f.y += extF.y;
			f.z += extF.z;
		}
		c_gsop.d_forces[d_i] = f;
	}
}

