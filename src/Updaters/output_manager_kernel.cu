#include "../gsop.h"
#include "output_manager.h"

__global__ void reset_temperature(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		c_gsop.d_T[d_i] = 0.0f;
    }
}



__global__ void computeRg_kernel(float* d_rgs){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		float4 r1 = c_gsop.d_coord[d_i];
		float rg = 0.0f;
		int i;
		int aminoCount = c_gsop.aminoCount/c_gsop.Ntr;
		int traj = d_i/aminoCount;
		for(i = traj*aminoCount; i < d_i; i++){
			float4 r2 = c_gsop.d_coord[i];
			r2.x -= r1.x;
			r2.y -= r1.y;
			r2.z -= r1.z;
			rg += r2.x*r2.x + r2.y*r2.y + r2.z*r2.z;
		}
		d_rgs[d_i] = rg;
	}
}

