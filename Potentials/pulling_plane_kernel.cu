/*
 * pulling_plane_kernel.cu
 *
 *  Created on: Jan 17, 2014
 *      Author: alekseenko
 */
#include "../gsop.cuh"


__global__ void pullingPlane_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		float4 f = c_gsop.d_forces[d_i];
		float4 extF = c_pullingPlane.d_extForces[d_i];
		if(extF.w == 1.0f){
			f.x = 0.0f;
			f.y = 0.0f;
			f.z = 0.0f;
		}
		// Pulled beads
		if(extF.w == 2.0f){
#ifdef NOTEXTURE
    		float4 coord = c_gsop.d_coord[d_i];
#else
	    	float4 coord = tex1Dfetch(t_coord, d_i);
#endif
            float3 norm = c_pullingPlane.pullVector;
            float dis = norm.x*coord.x + norm.y*coord.y + norm.z*coord.z + c_pullingPlane.d;
            dis = dis*c_pullingPlane.Ks;
            extF.x = dis*norm.x;
            extF.y = dis*norm.y;
            extF.z = dis*norm.z;
			f.x += extF.x;
			f.y += extF.y;
			f.z += extF.z;
		}
		c_gsop.d_forces[d_i] = f;
		c_pullingPlane.d_extForces[d_i] = extF;
	}
}

