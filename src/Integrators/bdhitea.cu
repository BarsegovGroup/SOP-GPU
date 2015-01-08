/*
 * bdhitea.cu
 *
 *  Created on: Mar 01, 2013
 *	  Author: alekseenko
 * 
 * The description of algorithm is available in Geyer & Winter, 2009 [doi:10.1063/1.3089668]
 * All the equations referenced here are from the very same paper
 */


#include "../gsop.cuh"
#include "ht.cu"
#include "bdhitea.h"
#include "langevin.h"

#ifdef TEA_TEXTURE
texture<float4, 1, cudaReadModeElementType> t_rforce;
texture<float4, 1, cudaReadModeElementType> t_mforce;
#endif

#include "bdhitea_kernel.cu"

// TODO: the code is nonintuitive and variables are scattered over structures

void createTeaIntegrator(){
	integrator = new TeaIntegrator();
}

TeaIntegrator::TeaIntegrator() : LangevinIntegrator() {
	this->name = "BDHI-TEA";
	if (gsop.minimizationOn){
		DIE("Using minimization with BDHI integrator is not supported and will never be.");
	}

	hc_tea.namino = sop.aminoCount;
	hc_tea.a = parameters::tea_a.get();
    this->exact = parameters::tea_exact.get();
    if (!this->exact){
    	this->capricious = parameters::tea_capricious.get();
	    this->epsilon_freq = parameters::tea_epsilon_freq.get();
    	this->epsmax = parameters::tea_epsmax.get();
    	this->unlisted = parameters::tea_unlisted.get();
    }

	cudaMalloc(&hc_tea.rforce, gsop.aminoCount * sizeof(float4));
	cudaMalloc(&hc_tea.mforce, gsop.aminoCount * sizeof(float4));
	cudaMalloc(&hc_tea.coords, gsop.aminoCount * sizeof(float4));
    if (!this->exact){
    	cudaMalloc(&hc_tea.d_epsilon, gsop.aminoCount * sizeof(float));
    	cudaMalloc(&hc_tea.d_ci, gsop.aminoCount * sizeof(float4));
    	cudaMalloc(&hc_tea.d_beta_ij, gsop.Ntr * sizeof(float));
	    this->h_epsilon = (float*) malloc(gsop.aminoCount * sizeof(float));
    	this->h_beta_ij = (float*) malloc(gsop.Ntr * sizeof(float));
    }else{
    	cudaMalloc(&hc_tea.d_tensor, gsop.Ntr * hc_tea.namino * 3 * hc_tea.namino * 3 * sizeof(float));
    }
	checkCUDAError();
	cudaMemcpyToSymbol(c_tea, &hc_tea, sizeof(TeaConstant), 0, cudaMemcpyHostToDevice);
	checkCUDAError();

    if(!this->exact){
	    updaters[updatersCount] = new TeaUpdater(this);
        updatersCount++;
	    printf("TEA integrator initialized, a = %f A, freq = %d steps; capricious mode: %s; using pairlist: %s; block size: %d\n", hc_tea.a, this->epsilon_freq, (this->capricious ? "on" : "off"), (this->unlisted ? "no" : "yes"), gsop.blockSize);
    }else{
    	printf("Cholesky integrator initialized, a = %f A\n", hc_tea.a);
    }

#ifndef TEA_TEXTURE
	cudaFuncSetCacheConfig(integrateTea_kernel, cudaFuncCachePreferL1);
	checkCUDAError();
#endif
}

void TeaIntegrator::integrate(){
	// Pregenerate random forces
	integrateTea_prepare<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
	// Integrate
#ifdef TEA_TEXTURE
	cudaBindTexture(0, t_rforce, hc_tea.rforce, gsop.aminoCount * sizeof(float4));
	cudaBindTexture(0, t_mforce, hc_tea.mforce, gsop.aminoCount * sizeof(float4));
#endif
    if (!this->exact){
    	if (this->unlisted)
    		integrateTea_kernel_unlisted<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
    	else
    		integrateTea_kernel<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
    }else{
        integrateCholesky_D<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>(hc_tea.d_tensor, 3*hc_tea.namino);
        integrateCholesky_decompose<<<gsop.Ntr, 3*hc_tea.namino>>>(hc_tea.d_tensor, 3*hc_tea.namino);
    /*
    float *d = (float*) malloc(gsop.Ntr * hc_tea.namino * 3 * hc_tea.namino * 3 * sizeof(flonat));
    cudaMemcpy(d, hc_tea.d_tensor, gsop.Ntr * hc_tea.namino * 3 * hc_tea.namino * 3 * sizeof(float), cudaMemcpyDeviceToHost);
    for (int i = 0; i < 3*5; ++i) {
        for (int j = 0; j < 3*5; ++j) {
            printf("%.3f ", d[i*3*hc_tea.namino + j]);
        }
        printf("\n");
    }
    */
        integrateCholesky_L<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>(hc_tea.d_tensor, 3*hc_tea.namino);
    }
#ifdef TEA_TEXTURE
	cudaUnbindTexture(t_rforce);
	cudaUnbindTexture(t_mforce);
#endif
	checkCUDAError();
}

TeaIntegrator::~TeaIntegrator(){
	cudaFree(hc_tea.rforce);
	cudaFree(hc_tea.mforce);
	cudaFree(hc_tea.coords);
    if (!this->exact){
    	cudaFree(hc_tea.d_epsilon);
    	cudaFree(hc_tea.d_beta_ij);
    	cudaFree(hc_tea.d_ci);
	    free(this->h_epsilon);
    	free(this->h_beta_ij);
    }else{
        cudaFree(hc_tea.d_tensor);
    }
}


TeaUpdater::TeaUpdater(TeaIntegrator *tea_integrator){
	this->name = "BDHI-TEA updater";
    this->tea = tea_integrator;
	this->frequency = tea_integrator->epsilon_freq;
}

void TeaUpdater::update(){
	const int update_epsilon = (gsop.step % this->frequency) == 0;
	const int N = sop.aminoCount;
	if (update_epsilon){
		// Calculate relative coupling
		if (tea->unlisted)
			integrateTea_epsilon_unlisted<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
		else
			integrateTea_epsilon<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
		// Dowload epsilon`s
		cudaMemcpy(this->tea->h_epsilon, hc_tea.d_epsilon, gsop.aminoCount * sizeof(float), cudaMemcpyDeviceToHost);
		checkCUDAError();
		for (int t = 0; t < gsop.Ntr; ++t){
			double epsilon = 0.0;
			for (int i = 0; i < N; ++i){
				epsilon += this->tea->h_epsilon[t*N + i];
			}
			epsilon /= 3.*N*(3.*N - 3.); // Averaging, off-diagonal elements only
			if (epsilon > 1.0){
				if (tea->capricious){
					DIE("HI tensor is not diagonal enough for trajectory %d: epsilon = %lf -> 1.0!\n", t, epsilon);
				}
				epsilon = 1.0;
			}
			if (epsilon > tea->epsmax){
				DIE("HI tensor is not diagonal enough for trajectory %d: epsilon = %lf > %f = tea_epsmax!\n", t, epsilon, tea->epsmax);
			}
			double a = (3.*N-1.)*epsilon*epsilon - (3.*N-2.)*epsilon;
			if (fabs(a) < 1e-7){ // To avoid 0/0 division in eq. (26) we explicitly handle small a's
				this->tea->h_beta_ij[t] = .5f;
				if(tea->capricious && hc_tea.a > 0.0f){
					DIE("HI tensor is too diagonal for trajectory %d: a = %lf, beta: %lf -> 0.5!\n", t, a, (1. - sqrt(1. - a)) / a);
				}
			} else {
				this->tea->h_beta_ij[t] = (1. - sqrt(1. - a)) / a; // eq. (26)
			}
			this->tea->h_epsilon[t] = epsilon; // We slowly overwrite the beginning of h_epsilon with per-trajectory epsilons to later output them
		}
		cudaMemcpy(hc_tea.d_beta_ij, this->tea->h_beta_ij, gsop.Ntr * sizeof(float), cudaMemcpyHostToDevice);
		checkCUDAError();
	}
}

TeaUpdater::~TeaUpdater(){
}

