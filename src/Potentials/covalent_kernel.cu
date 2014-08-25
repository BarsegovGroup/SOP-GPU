/*
 * covalent_kernel.cu
 *
 *  Created on: Jun 5, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"

struct CovalentConstant {
	GCovalentBond* d_bonds; // Same on device
	int* d_covalentCount;
    float kspring_cov, R_limit_sq;
};

CovalentConstant hc_covalent;
__device__ __constant__ CovalentConstant c_covalent;

/*
 * Kernel to compute FENE potential for covalent bonds
 */
__global__ void covalent_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
#ifdef NOTEXTURE
		float4 coord = c_gsop.d_coord[d_i];
#else
		float4 coord = tex1Dfetch(t_coord, d_i);
#endif
		float4 f = c_gsop.d_forces[d_i];
		GCovalentBond bond = c_covalent.d_bonds[d_i];
#ifdef NOTEXTURE
		float4 r2 = c_gsop.d_coord[bond.i2];
#else
		float4 r2 = tex1Dfetch(t_coord, bond.i2);
#endif
		int i;
		for(i = 1; i <= c_covalent.d_covalentCount[d_i]; i++){
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			r2.w = sqrtf(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z);
			r2.w -= bond.r0;
			r2.w = c_covalent.kspring_cov*r2.w/((r2.w + bond.r0)*(1.0f-r2.w*r2.w/(c_covalent.R_limit_sq)));
			bond = c_covalent.d_bonds[i*c_gsop.aminoCount + d_i];
			f.x += r2.w*r2.x;
			f.y += r2.w*r2.y;
			f.z += r2.w*r2.z;
#ifdef NOTEXTURE
			r2 = c_gsop.d_coord[bond.i2];
#else
			r2 = tex1Dfetch(t_coord, bond.i2);
#endif
		}
		c_gsop.d_forces[d_i] = f;
	}
}

/*
 * Kernel to compute FENE potential
 */
__global__ void covalentEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
#ifdef NOTEXTURE
		float4 coord = c_gsop.d_coord[d_i];
#else
		float4 coord = tex1Dfetch(t_coord, d_i);
#endif
		float pot = 0.0f;
		GCovalentBond bond = c_covalent.d_bonds[d_i];
#ifdef NOTEXTURE
		float4 r2 = c_gsop.d_coord[bond.i2];
#else
		float4 r2 = tex1Dfetch(t_coord, bond.i2);
#endif
		int i;
		for(i = 1; i <= c_covalent.d_covalentCount[d_i]; i++){
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			r2.w = sqrtf(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z);
			r2.w -= bond.r0;
			pot -= (c_covalent.kspring_cov/2.0)*(c_covalent.R_limit_sq)*__logf(1.0f - r2.w*r2.w/(c_covalent.R_limit_sq));
#ifdef NOTEXTURE
			r2 = c_gsop.d_coord[bond.i2];
#else
			r2 = tex1Dfetch(t_coord, bond.i2);
#endif
		}
		c_gsop.d_energies[d_i].x = pot;
	}
}
