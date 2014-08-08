/*
 * pairs_kernel.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"

/*
 * CUDA Kernel for computation of repulsive LJ. The general formula is:
 * U = el*(a/r)^6
 * el - strength of the interaction
 * a - repel distance
 * Computational procedure is done in N threads, where N is the total number of particles.
 * d_i - thread index - is an index of a particle, thread computes force for.
 * This kernel uses 15 registers, leaving enough resources for all 1024 possible threads on a multiprocessor.
 */
__global__ void pairs_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x; // Thread id, or particle id in this case (1 register)
	if(d_i < c_gsop.aminoCount){
#ifdef NOTEXTURE
		float4 coord = c_gsop.d_coord[d_i]; // Coordinates of the particle associated with a thread
#else
		float4 coord = tex1Dfetch(t_coord, d_i); // Coordinates of the particle associated with a thread read using texture reference (4 registers)
#endif
		float4 f = c_gsop.d_forces[d_i]; // Total force, acting on a particle (4 registers)
		int i2 = c_pairs.d_pairs[d_i];  // Reading the first interacting particle index from a map (1 register, reading ASAP)
		int i; // Loop counter (1 register)
		float4 r2; // Coordinates of second interacting particle (4 registers)
		for(i = 1; i <= c_pairs.d_pairsCount[d_i]; i++){ // Iterating over all interacting pairs
#ifdef NOTEXTURE
			r2 = c_gsop.d_coord[i2];
#else
			r2 = tex1Dfetch(t_coord, i2); // Reading the coordinates of the current interacting particle
#endif
			i2 = c_pairs.d_pairs[i*c_gsop.aminoCount + d_i]; // Index of the next interacting particle (reading ASAP)
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			r2.w = r2.x*r2.x+r2.y*r2.y+r2.z*r2.z;
			if(r2.w < c_pairs.pairsCutoff2){ // pairsCutoff2 = pairsCutoff*pairsCutoff
				r2.w = c_pairs.a2/r2.w; // a2 = a*a
				r2.w = r2.w*r2.w;
				r2.w = c_pairs.minus6elovera2*r2.w*r2.w; // minus6elovera2 = -6.0*el/(a*a)
				f.x += r2.x*r2.w;
				f.y += r2.y*r2.w;
				f.z += r2.z*r2.w;
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
#ifdef NOTEXTURE
		float4 coord = c_gsop.d_coord[d_i];
#else
		float4 coord = tex1Dfetch(t_coord, d_i);
#endif
		float pot = 0.0f;
		int i2 = c_pairs.d_pairs[d_i];
		int i;
		float4 r2;
		for(i = 1; i <= c_pairs.d_pairsCount[d_i]; i++){
#ifdef NOTEXTURE
			r2 = c_gsop.d_coord[i2];
#else
			r2 = tex1Dfetch(t_coord, i2);
#endif
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
		c_gsop.d_energies[d_i].z = pot;
	}
}
