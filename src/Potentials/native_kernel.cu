/*
 * native_kernel.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"

struct NativeConstant {
	int* d_native;
	int* d_nativeCount;
	GNativeParameters* d_nativeParameters;
	float* d_energies;
};

NativeConstant hc_native;
__device__ __constant__ NativeConstant c_native;

/*
 * CUDA Kernel for computation of full LJ (Native contacts). The general formula is:
 * U = eh*((r0/r)^12-2*(r0/r)^6)
 * eh - strength of the interaction (in general, depend on a particular contact)
 * r0 - equilibrium distance (in general, depend on a particular contact)
 * Computational procedure is done in N threads, where N is the total number of particles.
 * d_i - thread index - is an index of a particle, thread computes force for.
 */
__global__ void native_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x; // Thread id, or particle id in this case (1 register)
	if(d_i < c_gsop.aminoCount){
		float4 coord = c_gsop.d_coord[d_i];
		float4 f = c_gsop.d_forces[d_i]; // Force, acting on a particle (4 registers)
		int i2 = c_native.d_native[d_i]; // First interacting particle (1 register, reading ASAP)
		GNativeParameters par = c_native.d_nativeParameters[d_i]; // par.x = r0*r0, par.y = -12.0*eh/(r0*r0) (2 registers)
		int i; // Iterating index (1 register)
		float4 r2; // Coordinates of the second particle (4 registers)
		for(i = 1; i <= c_native.d_nativeCount[d_i]; i++){
			r2 = c_gsop.d_coord[i2];
			i2 = c_native.d_native[i*c_gsop.aminoCount + d_i]; // Index of the next interacting particle (reading ASAP)
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			r2.w = par.r02/(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z); // par.x = r0*r0
			float p6 = r2.w*r2.w*r2.w; // (1 register)
			r2.w = r2.w*(p6*p6 - p6)*par.minus12ehOverR02; // par.y = -12.0*eh/(r0*r0)
			par = c_native.d_nativeParameters[i*c_gsop.aminoCount + d_i]; // Parameters for the next interaction (reading ASAP)
			f.x += r2.x * r2.w;
			f.y += r2.y * r2.w;
			f.z += r2.z * r2.w;
		}
		c_gsop.d_forces[d_i] = f;
	}
}

/*
 * Kernel to compute energy native interactions for a given frame. Not optimized.
 */
__global__ void nativeEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		float4 coord = c_gsop.d_coord[d_i];
		int i2 = c_native.d_native[d_i];
		GNativeParameters par = c_native.d_nativeParameters[d_i]; // par.x = r0*r0, par.y = -12.0*eh/(r0*r0)
		float pot = 0.0f;
		int i;
		float4 r2;
		for(i = 1; i <= c_native.d_nativeCount[d_i]; i++){
			r2 = c_gsop.d_coord[i2];
			i2 = c_native.d_native[i*c_gsop.aminoCount + d_i];
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			r2.w = par.r02/(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z);
			float p6 = r2.w*r2.w*r2.w;
			float p12 = p6*p6;
			p6 = -2.0f*p6;
			par.minus12ehOverR02 = - par.minus12ehOverR02*par.r02/12.0; // eh (par.y = -12.0*eh/(r0*r0))
			pot +=  par.minus12ehOverR02*(p6 + p12);
			par = c_native.d_nativeParameters[i*c_gsop.aminoCount + d_i];
		}
		c_native.d_energies[d_i] = pot;
	}
}
