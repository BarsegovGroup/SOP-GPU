/*
 * langevin_kernel.cu
 *
 *  Created on: Jun 4, 2009
 *      Author: zhmurov
 */

struct LangevinConstant{
    float var;
    float hOverZeta;
    float tempNorm;
};

LangevinConstant hc_langevin;
__device__ __constant__ LangevinConstant c_langevin;

__global__ void integrateLangevin_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){

		float4 coord = c_gsop.d_coord[d_i];
		float4 f = c_gsop.d_forces[d_i];
		c_gsop.d_forces[d_i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);

		// Random force
		if(c_gsop.minimizationOn != 1){
			float var = c_langevin.var;
			float4 df = rforce(d_i);
			f.x += df.x*var;
			f.y += df.y*var;
			f.z += df.z*var;
		}

		if(c_gsop.pullingOn == 1){
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
		}

		if(c_gsop.pullingPlaneOn == 1){
			float4 extF = c_pullingPlane.d_extForces[d_i];
			if(extF.w == 1.0f){
				f.x = 0.0f;
				f.y = 0.0f;
				f.z = 0.0f;
            }
        }

		// Integration step
		float mult = c_langevin.hOverZeta;
		float3 dr = make_float3(mult*f.x, mult*f.y, mult*f.z);
		coord.x += dr.x;
		coord.y += dr.y;
		coord.z += dr.z;
		c_gsop.d_T[d_i] += c_langevin.tempNorm*(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		c_gsop.d_coord[d_i] = coord;
	}
}


