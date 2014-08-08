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
#include "bdhitea.cuh"
#include "langevin.cuh"
#include "bdhitea_kernel.cu"

void createTeaIntegrator(){
#ifndef TEA_TEXTURE
	cudaFuncSetCacheConfig(integrateTea_kernel, cudaFuncCachePreferL1);
#endif
	checkCUDAError();
	sprintf(teaIntegrator.name, "BDHI-TEA");
	teaIntegrator.integrate = &integrateTea;
	teaIntegrator.destroy = &deleteTeaIntegrator;
	integrator = &teaIntegrator;
	initTeaIntegrator();
}

void initTeaIntegrator(){
	if (gsop.minimizationOn){
		printf("Using minimization with BDHI integrator is not supported and will never be.");
		exit(-1);
	}
	initLangevinIntegrator();
	teaIntegrator.h = langevinIntegrator.h;
	tea.namino = sop.aminoCount;
	tea.a = getFloatParameter(BDHITEA_A_STRING, 1.8f, 1); // in Angstroms
    tea.exact = getYesNoParameter(BDHITEA_EXACT_STRING, 0, 1);
    if (!tea.exact){
    	tea.capricious = getYesNoParameter(BDHITEA_CAPRICIOUS_STRING, 1, 1); // Be paranoid about tensor values?
	    tea.epsilon_freq = getIntegerParameter(BDHITEA_EPSILONFREQ_STRING, 0, 0); // How often recalculate epsilon?
    	tea.epsmax = getFloatParameter(BDHITEA_EPSMAX_STRING, 999.f, 1); // Epsilon will never exceed 1, so epsmax=999 will never trigger halt by itself; used in capricious mode
    	if (getYesNoParameter(BDHITEA_UNLISTED_STRING, 1, 0)) { // use all-to-all or pairlists?
    		tea.unlisted = 1;
    	} else tea.unlisted = 0;
    }
    /*else{
        int arch = getCudaArchitecture();
        printf("Code compiled for CC %.1f\n", arch*.01f);
        if (getCudaArchitecture() >= 200){
            // Seems like I'm hitting sometihng related to https://devtalk.nvidia.com/default/topic/527307/-39-cicc-39-compilation-error-and-debug-flag/
            // Anyway, I'm too lazy to file a bug report
            printf("! Please compile sop-gpu for CC 1.3 or less\n");
            printf("! Compiling for CCs 2.0, 2.1 and 3.0 with nvcc 5.0, 5.5 and 6.0rc lead to incorrect behavior of Cholesky-based integrator\n");
            printf("!  (particularly, integrateCholesky_D kernel)");
            printf("! It will probably work with nvcc 6.0 or even 5.x if you are lucky\n");
            printf("! You can try to disable this safeguard by removing exit() call around %s:%d\n",__FILE__, __LINE__);
            printf("! Or you can try compiling with -G nvcc switch. It seems to work.\n");
            printf("! If something is wrong, trajectory NaN's immediately, so you will easily whether it works\n");
            printf("! It was nice talking to you. Goodbye and good luck.\n");
            exit(-1);
        }
    }*/

	cudaMalloc(&tea.rforce, gsop.aminoCount * sizeof(float4));
	cudaMalloc(&tea.mforce, gsop.aminoCount * sizeof(float4));
	cudaMalloc(&tea.coords, gsop.aminoCount * sizeof(float4));
    if (!tea.exact){
    	cudaMalloc(&tea.d_epsilon, gsop.aminoCount * sizeof(float));
    	cudaMalloc(&tea.d_ci, gsop.aminoCount * sizeof(float4));
    	cudaMalloc(&tea.d_beta_ij, Ntr * sizeof(float));
	    tea.h_epsilon = (float*) malloc(gsop.aminoCount * sizeof(float));
    	tea.h_beta_ij = (float*) malloc(Ntr * sizeof(float));
    }else{
    	cudaMalloc(&tea.d_tensor, Ntr * tea.namino * 3 * tea.namino * 3 * sizeof(float));
    }
	checkCUDAError();
	cudaMemcpyToSymbol(c_tea, &tea, sizeof(Tea), 0, cudaMemcpyHostToDevice);
	checkCUDAError();

    if(!tea.exact){
	    createTeaUpdater();
	    printf("TEA integrator initialized, a = %f A, freq = %d steps; capricious mode: %s; using pairlist: %s; block size: %d\n", tea.a, tea.epsilon_freq, (tea.capricious ? "on" : "off"), (tea.unlisted ? "no" : "yes"), BLOCK_SIZE);
    }else{
    	printf("Cholesky integrator initialized, a = %f A\n", tea.a);
    }
}

void integrateTea(){
	// Pregenerate random forces
	integrateTea_prepare<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
	// Integrate
#ifdef TEA_TEXTURE
	cudaBindTexture(0, t_rforce, tea.rforce, gsop.aminoCount * sizeof(float4));
	cudaBindTexture(0, t_mforce, tea.mforce, gsop.aminoCount * sizeof(float4));
#endif
    if (!tea.exact){
    	if (tea.unlisted)
    		integrateTea_kernel_unlisted<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
    	else
    		integrateTea_kernel<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
    }else{
        integrateCholesky_D<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>(tea.d_tensor, 3*tea.namino);
        integrateCholesky_decompose<<<Ntr, 3*tea.namino>>>(tea.d_tensor, 3*tea.namino);
    /*
    float *d = (float*) malloc(Ntr * tea.namino * 3 * tea.namino * 3 * sizeof(float));
    cudaMemcpy(d, tea.d_tensor, Ntr * tea.namino * 3 * tea.namino * 3 * sizeof(float), cudaMemcpyDeviceToHost);
    for (int i = 0; i < 3*5; ++i) {
        for (int j = 0; j < 3*5; ++j) {
            printf("%.3f ", d[i*3*tea.namino + j]);
        }
        printf("\n");
    }
    */
        integrateCholesky_L<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>(tea.d_tensor, 3*tea.namino);
    }
#ifdef TEA_TEXTURE
	cudaUnbindTexture(t_rforce);
	cudaUnbindTexture(t_mforce);
#endif
	checkCUDAError();
}

void deleteTeaIntegrator(){
	deleteLangevinIntegrator();
	cudaFree(tea.rforce);
	cudaFree(tea.mforce);
	cudaFree(tea.coords);
    if (!tea.exact){
    	cudaFree(tea.d_epsilon);
    	cudaFree(tea.d_beta_ij);
    	cudaFree(tea.d_ci);
	    free(tea.h_epsilon);
    	free(tea.h_beta_ij);
    }else{
        cudaFree(tea.d_tensor);
    }
}


void createTeaUpdater(){
	sprintf(teaUpdater.name, "BDHI-TEA updater");
	teaUpdater.update = &updateTea;
	teaUpdater.destroy = &deleteTeaUpdater;
	teaUpdater.frequency = tea.epsilon_freq;
	updaters[updatersCount] = &teaUpdater;
	updatersCount++;
	initTeaUpdater();
}

void initTeaUpdater(){
}

void updateTea(){
	const int update_epsilon = (step % tea.epsilon_freq) == 0;
	const int N = sop.aminoCount;
	if (update_epsilon){
		// Calculate relative coupling
		if (tea.unlisted)
			integrateTea_epsilon_unlisted<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
		else
			integrateTea_epsilon<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
		// Dowload epsilon`s
		cudaMemcpy(tea.h_epsilon, tea.d_epsilon, gsop.aminoCount * sizeof(float), cudaMemcpyDeviceToHost);
		checkCUDAError();
//		printf("epsilon: [ ");
		for (int t = 0; t < Ntr; ++t){
			double epsilon = 0.0;
			for (int i = 0; i < N; ++i){
				epsilon += tea.h_epsilon[t*N + i];
			}
			epsilon /= 3.*N*(3.*N - 3.); // Averaging, off-diagonal elements only
			if (epsilon > 1.0){
				if (tea.capricious){
					printf("HI tensor is not diagonal enough for trajectory %d: epsilon = %lf -> 1.0!\n", t, epsilon);
					exit(-1);
				}
				epsilon = 1.0;
			}
			if (epsilon > tea.epsmax){
				printf("HI tensor is not diagonal enough for trajectory %d: epsilon = %lf > %f = tea_epsmax!\n", t, epsilon, tea.epsmax);
				exit(-1);
			}
			double a = (3.*N-1.)*epsilon*epsilon - (3.*N-2.)*epsilon;
			if (fabs(a) < 1e-7){ // To avoid 0/0 division in eq. (26) we explicitly handle small a's
				tea.h_beta_ij[t] = .5f;
				if(tea.capricious && tea.a > 0.0f){
					printf("HI tensor is too diagonal for trajectory %d: a = %lf, beta: %lf -> 0.5!\n", t, a, (1. - sqrt(1. - a)) / a);
					exit(-1);
				}
			} else {
				tea.h_beta_ij[t] = (1. - sqrt(1. - a)) / a; // eq. (26)
			}
			tea.h_epsilon[t] = epsilon; // We slowly overwrite the beginning of h_epsilon with per-trajectory epsilons to later output them
//			printf("%lf ", epsilon);
		}
//		printf("]\n");
		cudaMemcpy(tea.d_beta_ij, tea.h_beta_ij, Ntr * sizeof(float), cudaMemcpyHostToDevice);
		checkCUDAError();
	}
}

void deleteTeaUpdater(){
}

