/*
 * bdhitea_kernel.cu
 *
 *  Created on: Mar 4, 2013
 *	  Author: alekseenko
 * 
 * All equations referenced here are from Geyer & Winter, 2009 [doi:10.1063/1.3089668]
 */



#include "../gsop.cuh"
#include "bdhitea.cuh"

struct float6 {
	float3 xx; // xx, xy, xz
	float3 yz; // yy, yz, zz
};

// Macroses to facilitate access to float6 members
#define _XX xx.x
#define _XY xx.y
#define _XZ xx.z
#define _YX xx.y
#define _YY yz.x
#define _YZ yz.y
#define _ZX xx.z
#define _ZY yz.y
#define _ZZ yz.z

__device__ inline void operator+=(float4 &a, const float4 &b){
	a.x+=b.x;
	a.y+=b.y;
	a.z+=b.z;
	a.w+=b.w;
}

__global__ void integrateTea_prepare(){
	// Precalculate random forces and apply pulling forces
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		// Random force
		const float var = c_langevin.var;
		float4 df = rforce(d_i);
		df.x *= var;
		df.y *= var;
		df.z *= var;
		c_tea.rforce[d_i] = df;

		// We take pulling forces into account here...
		if(c_gsop.pullingOn == 1){
			float4 f = c_gsop.d_forces[d_i];
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
			c_gsop.d_forces[d_i] = f;
		}
	}
}

__device__ inline float6 integrateTea_RPY(const float4 dr){
	// This functions requires `dr` to be normalized and have its original length in its 4th component
	// Returns Dij(dr)/Dii, where Dij is Rotne-Prager-Yamakawa tensor submatrix, eq. (3-5)
	float6 ret;
	const float ra = dr.w / c_tea.a;
	float coeffrr, coeffii;
	if (ra >= 2.f){
		coeffrr = 0.75f/ra * (1.f - 2.f/ra/ra);
		coeffii = 0.75f/ra * (1.f + 2.f/3.f/ra/ra);
	}else{
		coeffrr = 3.f*ra/32.f;
		coeffii = 1.f - 9.f*ra/32.f;
	}
	coeffrr /= 1.0f;
	coeffii /= 1.0f;
	ret._XX = dr.x*dr.x*coeffrr + coeffii;
	ret._XY = dr.x*dr.y*coeffrr;
	ret._XZ = dr.x*dr.z*coeffrr;
	ret._YY = dr.y*dr.y*coeffrr + coeffii;
	ret._YZ = dr.y*dr.z*coeffrr;
	ret._ZZ = dr.z*dr.z*coeffrr + coeffii;
	return ret;
}

__device__ inline float4 integrateTea_epsilon_local(const float4& coord1, const int idx2){
	// Functions calculates various statistics for sub-tensor responsible for interactions between two given particles
	// .x, .y, .z --- sum of squares of squares of normalized tensor components, 3 rows, used for C_i in eq. (19)
	// .w --- sum of normalized tensor components (\sum Dij/Dii), used for \epsilon in eq. (22)
	float4 dr = c_gsop.d_coord[idx2];
	dr.x -= coord1.x;
	dr.y -= coord1.y;
	dr.z -= coord1.z;
	dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
	dr.x /= dr.w;
	dr.y /= dr.w;
	dr.z /= dr.w;
	float6 d = integrateTea_RPY(dr);
	dr.w = d._XX + 2*d._XY + 2*d._XZ + d._YY + 2*d._YZ + d._ZZ; // Sum of all 3*3 tensor components
	dr.x = d._XX*d._XX + d._XY*d._XY + d._XZ*d._XZ;
	dr.y = d._YX*d._YX + d._YY*d._YY + d._YZ*d._YZ;
	dr.z = d._ZX*d._ZX + d._ZY*d._ZY + d._ZZ*d._ZZ;
	return dr;
}

__global__ void integrateTea_epsilon(){
	// Compute average coupling \epsilon and data for C_i
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		int i;
		float4 coord = c_gsop.d_coord[d_i];
		float4 sum = make_float4(0.f, 0.f, 0.f, 0.f);
		for(i = 0; i < c_covalent.d_covalentCount[d_i]; i++){ // HI interactions with covalent neighbors
			GCovalentBond bond = c_covalent.d_bonds[i*c_gsop.aminoCount + d_i];
			sum += integrateTea_epsilon_local(coord, bond.i2);
		}
		for(i = 0; i < c_native.d_nativeCount[d_i]; i++){ // HI interactions with native neighbors
			sum += integrateTea_epsilon_local(coord, c_native.d_native[i*c_gsop.aminoCount + d_i]);
		}
		for(i = 0; i < c_pairs.d_pairsCount[d_i]; i++){ // HI interactions with neighbors in pairlist
			sum += integrateTea_epsilon_local(coord, c_pairs.d_pairs[i*c_gsop.aminoCount + d_i]);
		}
		c_tea.d_ci[d_i] = make_float3(sum.x, sum.y, sum.z);
		c_tea.d_epsilon[d_i] = sum.w; // Should later be divided by number of degrees-of-freedom in single trajectory
	}
}

__device__ inline float4 integrateTea_force(const float4& coord1, const int idx2, float beta_ij, const float3& ci, const int idx1){
	// Calculate the effective force acting on particle with coordinates `coord1` from particle with index `idx2`
	// eq. (13,14,19)
	float4 dr = c_gsop.d_coord[idx2];
	dr.x -= coord1.x;
	dr.y -= coord1.y;
	dr.z -= coord1.z;
	dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
	dr.x /= dr.w;
	dr.y /= dr.w;
	dr.z /= dr.w;
	float4 f = c_gsop.d_forces[idx2];
	//f = make_float4(0.f, 0.f, 0.f, 0.f);
	float6 D = integrateTea_RPY(dr);
	/*
	printf("d[%d, %d]=%f\n", 3*idx1+0, 3*idx2+0, D._XX);
	printf("d[%d, %d]=%f\n", 3*idx1+0, 3*idx2+1, D._XY);
	printf("d[%d, %d]=%f\n", 3*idx1+0, 3*idx2+2, D._XZ);
	printf("d[%d, %d]=%f\n", 3*idx1+1, 3*idx2+0, D._YX);
	printf("d[%d, %d]=%f\n", 3*idx1+1, 3*idx2+1, D._YY);
	printf("d[%d, %d]=%f\n", 3*idx1+1, 3*idx2+2, D._YZ);
	printf("d[%d, %d]=%f\n", 3*idx1+2, 3*idx2+0, D._ZX);
	printf("d[%d, %d]=%f\n", 3*idx1+2, 3*idx2+1, D._ZY);
	printf("d[%d, %d]=%f\n", 3*idx1+2, 3*idx2+2, D._ZZ);
	*/
	float4 r = c_tea.rforce[idx2];
//	printf("%d: %.3f %.3f %.3f <-> %.3f %.3f %.3f\n", idx2, f.x, f.y, f.z, r.x, r.y, r.z);
//	printf("%f %f %f %f %f %f\n", D.xx.x, D.xx.y, D.xx.z, D.yz.x, D.yz.y, D.yz.z);
	f.x += beta_ij * r.x / sqrtf(1 + beta_ij*beta_ij*ci.x);
	f.y += beta_ij * r.y / sqrtf(1 + beta_ij*beta_ij*ci.y);
	f.z += beta_ij * r.z / sqrtf(1 + beta_ij*beta_ij*ci.z);
	return make_float4(D._XX*f.x + D._XY*f.y + D._XZ*f.z, D._YX*f.x + D._YY*f.y + D._YZ*f.z, D._ZX*f.x + D._ZY*f.y + D._ZZ*f.z, 0.f);
}

__global__ void integrateTea_kernel(){
	// Based on integrateLangevin_kernel from langevin_kernel.cu
	// Basically, the random force generation is delegated to another kernel, since all random forces must be already known by now
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		int i;
		float4 coord = c_gsop.d_coord[d_i];
		float4 f = c_gsop.d_forces[d_i]; // Previously calculated "Mechanical" force
//		printf("> %d : %f %f %f\n", d_i, f.x, f.y, f.z);
		float4 df;
		const int tr = d_i / c_tea.ntr;
		const float beta_ij = c_tea.d_beta_ij[tr];
		const float3 ci = c_tea.d_ci[d_i];

		// Calculate effective force from neighbors
//		printf("%d : %d %d %d\n",d_i, c_covalent.d_covalentCount[d_i], c_native.d_nativeCount[d_i], c_pairs.d_pairsCount[d_i]);
		for(i = 0; i < c_covalent.d_covalentCount[d_i]; i++){ // HI interactions with covalent neighbors
			GCovalentBond bond = c_covalent.d_bonds[i*c_gsop.aminoCount + d_i];
//			printf(">C %d-%d\n", d_i, bond.i2);
			df = integrateTea_force(coord, bond.i2, beta_ij, ci, d_i);
//			printf("HI %d %d : %f %f %f vs. %f %f %f\n", d_i, bond.i2, df.x, df.y, df.z, f.x, f.y, f.z);
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		for(i = 0; i < c_native.d_nativeCount[d_i]; i++){ // HI interactions with native neighbors
//			printf(">N %d-%d\n", d_i, c_native.d_native[i*c_gsop.aminoCount + d_i]);
			df = integrateTea_force(coord, c_native.d_native[i*c_gsop.aminoCount + d_i], beta_ij, ci, d_i);
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		for(i = 0; i < c_pairs.d_pairsCount[d_i]; i++){ // HI interactions with neighbors in pairlist
//			printf(">P %d-%d\n", d_i, c_pairs.d_pairs[i*c_gsop.aminoCount + d_i]);
			df = integrateTea_force(coord, c_pairs.d_pairs[i*c_gsop.aminoCount + d_i], beta_ij, ci, d_i);
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		df = c_tea.rforce[d_i];
		f.x += df.x;
		f.y += df.y;
		f.z += df.z;

		if(c_gsop.pullingOn == 1){
			float4 extF = c_pulling.d_extForces[d_i];
			if(extF.w == 1.0f){
				f.x = 0.0f;
				f.y = 0.0f;
				f.z = 0.0f;
			}
			// We have already added pulling force to pulled beads during preparation
		}

		// Integration step
		// We've replaced all forces with their `effective` counterparts, so this part of integration process stays unaltered
		const float mult = c_langevin.hOverZeta;
		const float3 dr = make_float3(mult*f.x, mult*f.y, mult*f.z);
		coord.x += dr.x;
		coord.y += dr.y;
		coord.z += dr.z;
		c_gsop.d_energies[d_i].w += c_langevin.tempNorm*(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		c_gsop.d_coord[d_i] = coord;
	}
}

__global__ void integrateTea_cleanup(){
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		// TODO: do force copying and zeroing in _prepare
		c_gsop.d_forces[d_i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
	}
}

#undef _XX
#undef _XY
#undef _XZ
#undef _YX
#undef _YY
#undef _YZ
#undef _ZX
#undef _ZY
#undef _ZZ

