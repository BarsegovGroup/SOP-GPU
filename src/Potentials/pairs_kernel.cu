/*
 * pairs_kernel.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"

struct PairsConstant {
	float pairsCutoff2;
	float a2;
	float el;
	float minus6elovera2;
	int* d_pairs;
	int* d_pairsCount;
	float* d_energies;
	int max_pairs;
};

PairsConstant hc_pairs;
__device__ __constant__ PairsConstant c_pairs;

/*
 * CUDA Kernel for computation of repulsive LJ. The general formula is:
 * U = el*(a/r)^6
 * fx1=-(dU)/(x1)=-6*el*(a/r)^6*(x2-x1)/r12^2
 * el - strength of the interaction
 * a - repel distance
 * Computational procedure is done in N threads, where N is the total number of particles.
 * d_i - thread index - is an index of a particle, thread computes force for.
 * This kernel uses 15 registers, leaving enough resources for all 1024 possible threads on a multiprocessor.
 */
__global__ void pairs_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x; // Thread id, or particle id in this case (1 register)
	if(d_i < c_gsop.aminoCount){
		float4 coord = c_gsop.d_coord[d_i]; // Coordinates of the particle associated with a thread
		float4 f = c_gsop.d_forces[d_i]; // Total force, acting on a particle (4 registers)
		int i2 = c_pairs.d_pairs[d_i];  // Reading the first interacting particle index from a map (1 register, reading ASAP)
		int i; // Loop counter (1 register)
		float4 r2; // Coordinates of second interacting particle (4 registers)
		for(i = 1; i <= c_pairs.d_pairsCount[d_i]; i++){ // Iterating over all interacting pairs
			r2 = c_gsop.d_coord[i2];
			i2 = c_pairs.d_pairs[i*c_gsop.aminoCount + d_i]; // Index of the next interacting particle (reading ASAP)
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			r2.w = r2.x*r2.x+r2.y*r2.y+r2.z*r2.z; // r2.w = (r2-r1)^2
			if(r2.w < c_pairs.pairsCutoff2){ // pairsCutoff2 = pairsCutoff*pairsCutoff
				r2.w = c_pairs.a2/r2.w; // a2 = a*a; r2.w = (a/r12)^2
				r2.w = r2.w*r2.w; // r2.w = (a/r12)^4
				r2.w = c_pairs.minus6elovera2*r2.w*r2.w; // minus6elovera2 = -6*el/(a*a); r2.w = -6*el*a^6/r12^8
				f.x += r2.x*r2.w; // fx = -6*el*(a/r12)^6*(x2-x1)/r12^2
				f.y += r2.y*r2.w; // fy = -6*el*(a/r12)^6*(y2-y1)/r12^2
				f.z += r2.z*r2.w; // fz = -6*el*(a/r12)^6*(z2-z1)/r12^2
			}
		}
		c_gsop.d_forces[d_i] = f; // Saving the resulting force
	}
}

/*
 * Kernel to compute energy of LJ-repulsive interaction at a given time-frame.
 * Since it is only called once in a while, it is not optimized.
 * This kernel source follows the source of force-computing kernel apart from
 * the actual value it computes and saves.
 * Note, that saving energies of all potentials one float4 variable for each particle
 * is an old gsop feature. Repulsive LJ is saved into 'z' component. This will be changed soon.
 */
__global__ void pairsEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		float4 coord = c_gsop.d_coord[d_i];
		float pot = 0.0f;
		int i2 = c_pairs.d_pairs[d_i];
		int i;
		float4 r2;
		for(i = 1; i <= c_pairs.d_pairsCount[d_i]; i++){
			r2 = c_gsop.d_coord[i2];
			i2 = c_pairs.d_pairs[i*c_gsop.aminoCount + d_i];
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			r2.w = r2.x*r2.x+r2.y*r2.y+r2.z*r2.z;
			if(r2.w < c_pairs.pairsCutoff2){
				r2.w = c_pairs.a2/r2.w; //Computing energy instead of force
				pot += c_pairs.el*r2.w*r2.w*r2.w;
			}
		}
		c_pairs.d_energies[d_i] = pot;
	}
}
