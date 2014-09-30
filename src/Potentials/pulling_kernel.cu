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

