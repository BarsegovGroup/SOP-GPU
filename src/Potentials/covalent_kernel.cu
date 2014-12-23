/*
 * covalent_kernel.cu
 *
 *  Created on: Jun 5, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"

struct CovalentConstant {
	GCovalentBond* d_bonds; // Information about the bonds (second particle index, equilibrium distance)
	float* d_energies; // Energies per particle (on GPU)
	int* d_covalentCount; // Number of covalent bonds per particle (on GPU)
    float kspring_cov, R_limit_sq; // FENE potential constants Ks and Rc^2
};

CovalentConstant hc_covalent;
__device__ __constant__ CovalentConstant c_covalent;

/*
 * Kernel to compute FENE potential for covalent bonds.
 * U_FENE=(Ks/2)*(Rc^2)*log(1-((r-r0)/Rc)^2);
 * Ks - spring constant (in kcal/mol)
 * Rc - bond extesibility (in A)
 * r - distance between bonded atoms
 * r0 - equilibrium (native) distance between bonded atoms (in A)
 * Each thread in the kernel computes value of FENE force for one particle (with index d_i).
 */
__global__ void covalent_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		float4 coord = c_gsop.d_coord[d_i]; // Coordinates of the first particle
		float4 f = c_gsop.d_forces[d_i]; // Force acting on the particle
		GCovalentBond bond = c_covalent.d_bonds[d_i]; // Bond parameters
		float4 r2 = c_gsop.d_coord[bond.i2]; // Coordinates of the second particle (reading ASAP)
		int i;
		for(i = 1; i <= c_covalent.d_covalentCount[d_i]; i++){
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			r2.w = sqrtf(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z); // r12=r2-r1
			r2.w -= bond.r0; // r12-r0
			r2.w = c_covalent.kspring_cov*r2.w/((r2.w + bond.r0)*(1.0f-r2.w*r2.w/(c_covalent.R_limit_sq))); //Ks/(r12*(1-((r12-r0)/Rc)^2))
			bond = c_covalent.d_bonds[i*c_gsop.aminoCount + d_i]; // Reading the information about the second particle ASAP
			f.x += r2.w*r2.x; // fx = (x2-x1)*Ks/(r12*(1-((r12-r0)/Rc)^2))
			f.y += r2.w*r2.y; // fy = (y2-y1)*Ks/(r12*(1-((r12-r0)/Rc)^2))
			f.z += r2.w*r2.z; // fz = (z2-z1)*Ks/(r12*(1-((r12-r0)/Rc)^2))
			r2 = c_gsop.d_coord[bond.i2]; // Reading the coordinates of the second particle ASAP
		}
		c_gsop.d_forces[d_i] = f; // Saving the resulting force
	}
}

/*
 * Kernel to compute FENE potential energy
 * Each thread in the kernel computes value of FENE energy for one particle (with index d_i).
 */
__global__ void covalentEnergy_kernel(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		float4 coord = c_gsop.d_coord[d_i]; // Coordinates of the first particle
		float pot = 0.0f; // Place to accumulate potential energy
		GCovalentBond bond = c_covalent.d_bonds[d_i]; // Information about the bond (reading ASAP)
		float4 r2 = c_gsop.d_coord[bond.i2]; // Coordinates of the second particle (reading ASAP)
		int i;
		for(i = 1; i <= c_covalent.d_covalentCount[d_i]; i++){
			r2.x -= coord.x;
			r2.y -= coord.y;
			r2.z -= coord.z;
			r2.w = sqrtf(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z); // r12 = r2-r1
			r2.w -= bond.r0; // r12-r0
			pot -= (c_covalent.kspring_cov/2.0)*(c_covalent.R_limit_sq)*__logf(1.0f - r2.w*r2.w/(c_covalent.R_limit_sq)); // u=(Ks/2)*(Rc^2)*log(1-((r12-r0)/Rc)^2);
			bond = c_covalent.d_bonds[i*c_gsop.aminoCount + d_i]; // Information about the second bond (reading ASAP)
			r2 = c_gsop.d_coord[bond.i2]; // Coordinates of the second particle (Reading ASAP)
		}
		c_covalent.d_energies[d_i] = pot; // Saving the resulting energy per particle (summation is done on CPU)
	}
}
