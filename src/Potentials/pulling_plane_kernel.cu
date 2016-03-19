/*
 * pulling_plane_kernel.cu
 *
 *  Created on: Jan 17, 2014
 *      Author: alekseenko
 */
#include "../gsop.cuh"

struct PullingPlaneConstant{
	float* d_extForces;
	int* d_masks;
	float* d_r0;
	float3 pullVector;
	float3 planeCoord0;
	float* d_planeDispl;
	float Ks;
};

PullingPlaneConstant hc_pullingPlane;
__device__ __constant__ PullingPlaneConstant c_pullingPlane;
// TODO: add function to potential to copy these parameters to constant memory

__global__ void pullingPlane_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		int mask = c_pullingPlane.d_masks[d_i];
		if(mask != 0){
			float4 f = c_gsop.d_forces[d_i];
			float4 coord = c_gsop.d_coord[d_i];
			float3 planeCoord = c_pullingPlane.planeCoord0;
			if(mask == 2){
				int aminosPerTraj = c_gsop.aminoCount/c_gsop.Ntr;
				int traj = d_i/aminosPerTraj;
				float d = c_pullingPlane.d_planeDispl[traj];
				planeCoord.x += d*c_pullingPlane.pullVector.x;
				planeCoord.y += d*c_pullingPlane.pullVector.y;
				planeCoord.z += d*c_pullingPlane.pullVector.z;
			}

			float dx = coord.x - planeCoord.x;
			float dy = coord.y - planeCoord.y;
			float dz = coord.z - planeCoord.z;
			float dr = dx*c_pullingPlane.pullVector.x + dy*c_pullingPlane.pullVector.y + dz*c_pullingPlane.pullVector.z;
			float df = -c_pullingPlane.Ks*(dr - c_pullingPlane.d_r0[d_i]);
			f.x += df*c_pullingPlane.pullVector.x;
			f.y += df*c_pullingPlane.pullVector.y;
			f.z += df*c_pullingPlane.pullVector.z;
			c_pullingPlane.d_extForces[d_i] = df;
			c_gsop.d_forces[d_i] = f;
		}
	}
}

