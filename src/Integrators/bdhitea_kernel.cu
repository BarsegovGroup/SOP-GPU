/*
 * bdhitea_kernel.cu
 *
 *  Created on: Mar 4, 2013
 *	  Author: alekseenko
 * 
 * All equations referenced here are from Geyer & Winter, 2009 [doi:10.1063/1.3089668]
 */

//#define PRINT_HI_TENSORS // Enable to print tensor values. Should never be used in production. Calling `cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 900000000);` is recommended if you do not want to lose any values!

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
	a.x+=b.x; a.y+=b.y; a.z+=b.z; a.w+=b.w;
}
__device__ inline void operator+=(double4 &a, const float4 &b){
	a.x+=b.x; a.y+=b.y; a.z+=b.z; a.w+=b.w;
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
		float4 f = c_gsop.d_forces[d_i];
		if(c_gsop.pullingOn){
			float4 extF = c_pulling.d_extForces[d_i];
			if(extF.w == 1.0f){
				f.x = 0.0f;
				f.y = 0.0f;
				f.z = 0.0f;
			}
			if(extF.w == 2.0f){
				f.x += extF.x;
				f.y += extF.y;
				f.z += extF.z;
			}
		}
		// Copy forces and coordinates to auxillary arrays to avoid races during integration phase
		c_tea.mforce[d_i] = f;
		c_gsop.d_forces[d_i] = make_float4(0.f, 0.f, 0.f, 0.f);
		c_tea.coords[d_i] = c_gsop.d_coord[d_i];
	}
}

__device__ inline float6 integrateTea_RPY(const float4& dr){
	// This functions requires `dr` to be normalized and have its original length in its W component
	// Returns Dij(dr)/Dii, where Dij is Rotne-Prager-Yamakawa tensor submatrix, eq. (3-5)
	float6 ret;
	const float ra = dr.w / c_tea.a;
	float coeffrr, coeffii;
	if (ra > 2.f){
		coeffrr = 0.75f/ra * (1.f - 2.f/ra/ra);
		coeffii = 0.75f/ra * (1.f + 2.f/3.f/ra/ra);
	}else{
		coeffrr = 3.f*ra/32.f;
		coeffii = 1.f - 9.f*ra/32.f;
	}
	ret._XX = dr.x*dr.x*coeffrr + coeffii;
	ret._XY = dr.x*dr.y*coeffrr;
	ret._XZ = dr.x*dr.z*coeffrr;
	ret._YY = dr.y*dr.y*coeffrr + coeffii;
	ret._YZ = dr.y*dr.z*coeffrr;
	ret._ZZ = dr.z*dr.z*coeffrr + coeffii;
	return ret;
}

__global__ void integrateCholesky_D(float* D, int stride){
	// First step of Cholesky-based integration: 
    //   compute D matrices and calculate D*F
    // assert (stride >= c_tea.namino * 3)
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		const float4 coord = c_gsop.d_coord[d_i];
	    float4 f = c_tea.mforce[d_i];
		const int traj = (int)(d_i / c_tea.namino);
        const int i0 = traj * c_tea.namino;
        // Pointer to the start of matrix for current trajectory
        int didx = traj * stride * 3 * c_tea.namino;
        // Pointer to the start of 3-column for current bead
        didx += (d_i-i0)*3;
		int i;
		for(i = i0; i < i0 + c_tea.namino; i++){
            float4 dr = c_gsop.d_coord[i];
        	dr.x -= coord.x;
            dr.y -= coord.y;
            dr.z -= coord.z;
            dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
            if (i == d_i) {
                D[didx+0] = 1.; // Probably should cast to float3, but I believe in NVCC
                D[didx+1] = 0.;
                D[didx+2] = 0.;
                didx += stride; // Next row
                D[didx+0] = 0.;
                D[didx+1] = 1.;
                D[didx+2] = 0.;
                didx += stride; // Next row
                D[didx+0] = 0.;
                D[didx+1] = 0.;
                D[didx+2] = 1.;
                didx += stride; // Next row              
            }else{
                dr.x /= dr.w;
                dr.y /= dr.w;
                dr.z /= dr.w;
                float6 d = integrateTea_RPY(dr);
                D[didx+0] = d._XX; // Probably should cast to float3, but I believe in NVCC
                D[didx+1] = d._XY;
                D[didx+2] = d._XZ;
                didx += stride; // Next row
                D[didx+0] = d._YX;
                D[didx+1] = d._YY;
                D[didx+2] = d._YZ;
                didx += stride; // Next row
                D[didx+0] = d._ZX;
                D[didx+1] = d._ZY;
                D[didx+2] = d._ZZ;
                didx += stride; // Next row
#ifdef TEA_TEXTURE
            	float4 df = tex1Dfetch(t_mforce, i);
#else
            	float4 df = c_tea.mforce[i];
#endif
                f.x += d._XX * df.x + d._XY * df.y + d._XZ * df.z;
                f.y += d._YX * df.x + d._YY * df.y + d._YZ * df.z;
                f.z += d._ZX * df.x + d._ZY * df.y + d._ZZ * df.z;
            }
		}
        c_tea.mforce[d_i] = f;
    }
}

__global__ void integrateCholesky_decompose(float* D, int stride){
    // Second step of Cholesky-based integration
    //   compute Cholesky decompositions of several matrices
    // One block for one submatrix, 3*c_tea.namino threads in each block, 
    //  up to CHOLSIZE threads per block (due to shared mem)
    // This version is not ueber-optimal, but seems to perform on par with CuBCholS, 
    //  and its authors were proud enough to make a poster talk about it

    // assert(stride >= 3*c_tea.namino); assert(stride <= CHOLSIZE)
    __shared__ float diag;
    __shared__ float row[CHOLSIZE];
    float val;
    const int tid = threadIdx.x;
    const int m = c_tea.namino * 3;
    float* M = D + (blockIdx.y*gridDim.x + blockIdx.x)*m*stride;
    int j,k;
    for(j = 0; j < m; j++){
        row[tid] = M[j*stride + tid];
        __syncthreads();
        val = row[tid];
        if (tid >= j){
            for (k = 0; k < j; ++k){
                val -= M[k*stride + tid] * row[k];
            }
        }
        if (tid == j) {
           diag = sqrtf(val);
        }
        __syncthreads();
        if (tid >= j){
            val = val / diag;
            M[tid*stride + j] = val;
            M[j*stride + tid] = val;
        }
        __syncthreads(); // CUDA THREADS SYNCHRONIZE KAMUI SENKETSU!!!
    }
}

__device__ inline float4 integrateTea_epsilon_local(const float4& coord1, const int idx2){
	// Function calculates various statistics for sub-tensor responsible for interactions between two given particles
	// .x, .y, .z --- sum of squares of normalized tensor components, 3 rows, used for C_i in eq. (19)
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
		c_tea.d_ci[d_i] = sum;//make_float3(sum.x, sum.y, sum.z);
		c_tea.d_epsilon[d_i] = sum.w; // Should later be divided by number of non-diagonal degrees-of-freedom in single trajectory
	}
}

__global__ void integrateTea_epsilon_unlisted(){
	// Like integrateTea_epsilon, but calculate all-vs-all
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		const int i0 = ((int)(d_i / c_tea.namino))*c_tea.namino;
		int i;
		float4 coord = c_gsop.d_coord[d_i];
		float4 sum = make_float4(0.f, 0.f, 0.f, 0.f);
		for(i = i0; i < i0 + c_tea.namino; i++){
			if (i != d_i)
				sum += integrateTea_epsilon_local(coord, i);
		}
		c_tea.d_ci[d_i] = sum;//make_float3(sum.x, sum.y, sum.z);
		c_tea.d_epsilon[d_i] = sum.w; // Should later be divided by number of non-diagonal degrees-of-freedom in single trajectory
	}
}


__device__ inline float4 integrateTea_force(const float4& coord1, const int idx2, const float3& ci, const int idx1){
	// Calculate the effective force acting on particle with coordinates `coord1` from particle with index `idx2`
	// eq. (13,14,19)
#ifndef NOTEXTURE
	float4 dr = tex1Dfetch(t_coord, idx2);
#else
	float4 dr = c_tea.coords[idx2];
#endif
	dr.x -= coord1.x;
	dr.y -= coord1.y;
	dr.z -= coord1.z;
	dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
	dr.x /= dr.w;
	dr.y /= dr.w;
	dr.z /= dr.w;
#ifdef TEA_TEXTURE
	float4 f = tex1Dfetch(t_mforce, idx2);
	float4 r = tex1Dfetch(t_rforce, idx2);
#else
	float4 f = c_tea.mforce[idx2];
	float4 r = c_tea.rforce[idx2];
#endif
	f.x += r.x * ci.x;
	f.y += r.y * ci.y;
	f.z += r.z * ci.z;
	float6 D = integrateTea_RPY(dr);
#ifdef PRINT_HI_TENSORS
	printf("d[%d, %d]=%f\n", 3*idx1+0, 3*idx2+0, D._XX);
	printf("d[%d, %d]=%f\n", 3*idx1+0, 3*idx2+1, D._XY);
	printf("d[%d, %d]=%f\n", 3*idx1+0, 3*idx2+2, D._XZ);
	printf("d[%d, %d]=%f\n", 3*idx1+1, 3*idx2+0, D._YX);
	printf("d[%d, %d]=%f\n", 3*idx1+1, 3*idx2+1, D._YY);
	printf("d[%d, %d]=%f\n", 3*idx1+1, 3*idx2+2, D._YZ);
	printf("d[%d, %d]=%f\n", 3*idx1+2, 3*idx2+0, D._ZX);
	printf("d[%d, %d]=%f\n", 3*idx1+2, 3*idx2+1, D._ZY);
	printf("d[%d, %d]=%f\n", 3*idx1+2, 3*idx2+2, D._ZZ);
	printf("b[%d, %d]=%f\n", 3*idx1+0, 3*idx2+0, D._XX*ci.x);
	printf("b[%d, %d]=%f\n", 3*idx1+0, 3*idx2+1, D._XY*ci.x);
	printf("b[%d, %d]=%f\n", 3*idx1+0, 3*idx2+2, D._XZ*ci.x);
	printf("b[%d, %d]=%f\n", 3*idx1+1, 3*idx2+0, D._YX*ci.y);
	printf("b[%d, %d]=%f\n", 3*idx1+1, 3*idx2+1, D._YY*ci.y);
	printf("b[%d, %d]=%f\n", 3*idx1+1, 3*idx2+2, D._YZ*ci.y);
	printf("b[%d, %d]=%f\n", 3*idx1+2, 3*idx2+0, D._ZX*ci.z);
	printf("b[%d, %d]=%f\n", 3*idx1+2, 3*idx2+1, D._ZY*ci.z);
	printf("b[%d, %d]=%f\n", 3*idx1+2, 3*idx2+2, D._ZZ*ci.z);
#endif
	return make_float4( D._XX*f.x + D._XY*f.y + D._XZ*f.z, 
						D._YX*f.x + D._YY*f.y + D._YZ*f.z, 
						D._ZX*f.x + D._ZY*f.y + D._ZZ*f.z, 0.f);
}

__global__ void integrateCholesky_L(const float* L, int stride){
    // Third and final step of Cholesky-based integration:
    //   calculate B*R and final displacement
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		int i;
		float4 f = c_tea.mforce[d_i];
		const int traj = (int)(d_i / c_tea.namino);
        const int i0 = traj * c_tea.namino;
        // Pointer to the start of matrix for current trajectory
        const float* L0 = L + traj * stride * 3 * c_tea.namino;
        // Pointer to the start of 3-column for current bead
        L0 += (d_i-i0)*3;
		for(i = i0; i <= d_i; i++){
#ifdef TEA_TEXTURE
        	float4 df = tex1Dfetch(t_rforce, i);
#else
        	float4 df = c_tea.rforce[i];
#endif
            float6 d;
            d._XX = L0[0];
            d._XY = L0[1];
            d._XZ = L0[2];
            L0 += stride;
            d._YY = L0[1];
            d._YZ = L0[2];
            L0 += stride;
            d._ZZ = L0[2];
            L0 += stride;
            f.x += d._XX * df.x + d._XY * df.y + d._XZ * df.z;
            f.y += d._YX * df.x + d._YY * df.y + d._YZ * df.z;
            f.z += d._ZX * df.x + d._ZY * df.y + d._ZZ * df.z;
        //    printf("%d %d - %.3f %.3f %.3f - %.3f %.3f %.3f %.3f %.3f %.3f\n", d_i, i, f.x, f.y, f.z, d._XX, d._XY, d._XZ, d._YY, d._YZ, d._ZZ);
		}
        if(c_gsop.pullingOn){
			float4 extF = c_pulling.d_extForces[d_i];
			if(extF.w == 1.0f){
				f.x = 0.0f;
				f.y = 0.0f;
				f.z = 0.0f;
			}
			// We have already added pulling force to pulled beads during preparation
		}

		// Integration step
		// We've replaced all forces with their `effective` counterparts, so this part of integration process stays the same as in simple langevin integrator
		const float mult = c_langevin.hOverZeta;
		const float3 dr = make_float3(mult*f.x, mult*f.y, mult*f.z);
		float4 coord = c_tea.coords[d_i];
		coord.x += dr.x;
		coord.y += dr.y;
		coord.z += dr.z;
		c_gsop.d_coord[d_i] = coord;
		// Update energies
		coord = c_gsop.d_energies[d_i];
		coord.w += c_langevin.tempNorm*(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		c_gsop.d_energies[d_i] = coord;
	}
}


__global__ void integrateTea_kernel(){
	// Based on integrateLangevin_kernel from langevin_kernel.cu
	// The random force generation is delegated to another kernel, since all random forces must be already known by now
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		int i;
		float4 coord = c_tea.d_ci[d_i]; // Not coord yet!
		float4 f = c_tea.mforce[d_i];
		float4 df = c_tea.rforce[d_i];
		const int tr = d_i / c_tea.namino;
		const float beta_ij = c_tea.d_beta_ij[tr];
		float3 ci;
		// Make ci to be actual C_i
		coord.w = beta_ij*beta_ij;
		ci.x = 1.f/sqrtf(1.f + coord.w * coord.x);
		ci.y = 1.f/sqrtf(1.f + coord.w * coord.y);
		ci.z = 1.f/sqrtf(1.f + coord.w * coord.z);
		coord = c_tea.coords[d_i]; // And now it's actually bead coordinates
		f.x += df.x * ci.x;
		f.y += df.y * ci.y;
		f.z += df.z * ci.z;
		ci.x *= beta_ij;
		ci.y *= beta_ij;
		ci.z *= beta_ij;

		// Calculate effective force from neighbors
		const int nc = c_covalent.d_covalentCount[d_i];
		for(i = 0; i < nc; i++){ // HI interactions with covalent neighbors
			GCovalentBond bond = c_covalent.d_bonds[i*c_gsop.aminoCount + d_i];
			df = integrateTea_force(coord, bond.i2, ci, d_i);
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		const int nn = c_native.d_nativeCount[d_i];
		for(i = 0; i < nn; i++){ // HI interactions with native neighbors
			df = integrateTea_force(coord, c_native.d_native[i*c_gsop.aminoCount + d_i], ci, d_i);
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		const int np = c_pairs.d_pairsCount[d_i];
		for(i = 0; i < np; i++){ // HI interactions with neighbors in pairlist
			df = integrateTea_force(coord, c_pairs.d_pairs[i*c_gsop.aminoCount + d_i], ci, d_i);
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
#ifdef PRINT_HI_TENSORS
		printf("d[%d, %d]=%f\n", 3*d_i+0, 3*d_i+0, 1.0);
		printf("d[%d, %d]=%f\n", 3*d_i+0, 3*d_i+1, 0.0);
		printf("d[%d, %d]=%f\n", 3*d_i+0, 3*d_i+2, 0.0);
		printf("d[%d, %d]=%f\n", 3*d_i+1, 3*d_i+0, 0.0);
		printf("d[%d, %d]=%f\n", 3*d_i+1, 3*d_i+1, 1.0);
		printf("d[%d, %d]=%f\n", 3*d_i+1, 3*d_i+2, 0.0);
		printf("d[%d, %d]=%f\n", 3*d_i+2, 3*d_i+0, 0.0);
		printf("d[%d, %d]=%f\n", 3*d_i+2, 3*d_i+1, 0.0);
		printf("d[%d, %d]=%f\n", 3*d_i+2, 3*d_i+2, 1.0);
		printf("b[%d, %d]=%f\n", 3*d_i+0, 3*d_i+0, ci.x/beta_ij);
		printf("b[%d, %d]=%f\n", 3*d_i+0, 3*d_i+1, 0.0);
		printf("b[%d, %d]=%f\n", 3*d_i+0, 3*d_i+2, 0.0);
		printf("b[%d, %d]=%f\n", 3*d_i+1, 3*d_i+0, 0.0);
		printf("b[%d, %d]=%f\n", 3*d_i+1, 3*d_i+1, ci.y/beta_ij);
		printf("b[%d, %d]=%f\n", 3*d_i+1, 3*d_i+2, 0.0);
		printf("b[%d, %d]=%f\n", 3*d_i+2, 3*d_i+0, 0.0);
		printf("b[%d, %d]=%f\n", 3*d_i+2, 3*d_i+1, 0.0);
		printf("b[%d, %d]=%f\n", 3*d_i+2, 3*d_i+2, ci.z/beta_ij);
#endif

		if(c_gsop.pullingOn){
			float4 extF = c_pulling.d_extForces[d_i];
			if(extF.w == 1.0f){
				f.x = 0.0f;
				f.y = 0.0f;
				f.z = 0.0f;
			}
			// We have already added pulling force to pulled beads during preparation
		}

		// Integration step
		// We've replaced all forces with their `effective` counterparts, so this part of integration process stays the same as in simple langevin integrator
		const float mult = c_langevin.hOverZeta;
		const float3 dr = make_float3(mult*f.x, mult*f.y, mult*f.z);
		coord.x += dr.x;
		coord.y += dr.y;
		coord.z += dr.z;
		c_gsop.d_coord[d_i] = coord;
		// Update energies
		coord = c_gsop.d_energies[d_i];
		coord.w += c_langevin.tempNorm*(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		c_gsop.d_energies[d_i] = coord;
	}
}


__global__ void integrateTea_kernel_unlisted(){
	// Pairist-free version of  integrateTea_kernel
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		int i;
		float4 coord = c_tea.d_ci[d_i]; // Not coord yet!
		float4 f = c_tea.mforce[d_i];
		float4 df = c_tea.rforce[d_i];
		const int tr = d_i / c_tea.namino;
		const float beta_ij = c_tea.d_beta_ij[tr];
		float3 ci;
		// Make ci to be actual C_i
		coord.w = beta_ij*beta_ij;
		ci.x = 1.f/sqrtf(1.f + coord.w * coord.x);
		ci.y = 1.f/sqrtf(1.f + coord.w * coord.y);
		ci.z = 1.f/sqrtf(1.f + coord.w * coord.z);
		coord = c_tea.coords[d_i]; // And now it's actually bead coordinates
		f.x += df.x * ci.x;
		f.y += df.y * ci.y;
		f.z += df.z * ci.z;
		ci.x *= beta_ij;
		ci.y *= beta_ij;
		ci.z *= beta_ij;

		// Calculate effective force
		const int i0 = ((int)(d_i / c_tea.namino))*c_tea.namino;
		for(i = i0; i < i0 + c_tea.namino; i++){
			if (i == d_i) continue;
			df = integrateTea_force(coord, i, ci, d_i);
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		if(c_gsop.pullingOn){
			float4 extF = c_pulling.d_extForces[d_i];
			if(extF.w == 1.0f){
				f.x = 0.0f;
				f.y = 0.0f;
				f.z = 0.0f;
			}
			// We have already added pulling force to pulled beads during preparation
		}

		// Integration step
		// We've replaced all forces with their `effective` counterparts, so this part of integration process stays the same as in simple langevin integrator
		const float mult = c_langevin.hOverZeta;
		const float3 dr = make_float3(mult*f.x, mult*f.y, mult*f.z);
		coord.x += dr.x;
		coord.y += dr.y;
		coord.z += dr.z;
		c_gsop.d_coord[d_i] = coord;
		// Update energies
		coord = c_gsop.d_energies[d_i];
		coord.w += c_langevin.tempNorm*(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		c_gsop.d_energies[d_i] = coord;
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

