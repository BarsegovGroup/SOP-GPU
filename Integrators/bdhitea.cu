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

	tea.namino = sop.aminoCount;
	tea.epsilon_freq = getIntegerParameter(BDHITEA_EPSILONFREQ_STRING, 0, 0);
	tea.a = getFloatParameter(BDHITEA_A_STRING, 1.2f, 1);
	tea.capricious = getYesNoParameter(BDHITEA_CAPRICIOUS_STRING, 1, 1);
	tea.epsmax = getFloatParameter(BDHITEA_EPSMAX_STRING, 413.f, 1); // Epsilon will never exceed 1, so epsmax=413 will never trigger halt by itself

	cudaMalloc(&tea.rforce, gsop.aminoCount * sizeof(float4));
	cudaMalloc(&tea.mforce, gsop.aminoCount * sizeof(float4));
	cudaMalloc(&tea.coords, gsop.aminoCount * sizeof(float4));
	cudaMalloc(&tea.d_epsilon, gsop.aminoCount * sizeof(float));
	cudaMalloc(&tea.d_ci, gsop.aminoCount * sizeof(float4));
	cudaMalloc(&tea.d_beta_ij, Ntr * sizeof(float));
	checkCUDAError();
	tea.h_epsilon = (float*) malloc(gsop.aminoCount * sizeof(float));
	tea.h_beta_ij = (float*) malloc(Ntr * sizeof(float));
	cudaMemcpyToSymbol(c_tea, &tea, sizeof(Tea), 0, cudaMemcpyHostToDevice);
	checkCUDAError();

	createTeaUpdater();
	printf("TEA integrator initialized, a = %f A, freq = %d steps; capricious mode %s; block size %d\n", tea.a, tea.epsilon_freq, (tea.capricious ? "on" : "off"), BLOCK_SIZE);
}

void integrateTea(){
	// Pregenerate random forces
	integrateTea_prepare<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
	// Integrate
#ifdef TEA_TEXTURE
	cudaBindTexture(0, t_rforce, tea.rforce, gsop.aminoCount * sizeof(float4));
	cudaBindTexture(0, t_mforce, tea.mforce, gsop.aminoCount * sizeof(float4));
#endif
	integrateTea_kernel<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
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
	cudaFree(tea.d_epsilon);
	cudaFree(tea.d_beta_ij);
	cudaFree(tea.d_ci);
	free(tea.h_epsilon);
	free(tea.h_beta_ij);
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
				epsilon = 1.0;
				if (tea.capricious){
					printf("HI tensor is not diagonal enough for trajectory %d: epsilon = %lf -> 1.0!\n", t, epsilon);
					exit(-1);
				}
			}
			if (epsilon > tea.epsmax){
				printf("HI tensor is not diagonal enough for trajectory %d: epsilon = %lf > %f = tea_epsmax!\n", t, epsilon, tea.epsmax);
				exit(-1);
			}
			double a = (3.*N-1.)*epsilon*epsilon - (3.*N-2.)*epsilon;
			if (fabs(a) < 1e-5){ // To avoid 0/0 division in eq. (26) we explicitly handle small a's
				tea.h_beta_ij[t] = .5f;
				if(tea.capricious){
					printf("HI tensor is too diagonal for trajectory %d: a = %lf, beta: %lf -> 0.5!\n", t, a, (1. - sqrt(1. - a)) / a);
					exit(-1);
				}
			} else {
				tea.h_beta_ij[t] = (1. - sqrt(1. - a)) / a; // eq. (26)
			}
//			printf("%.3f ", epsilon);
		}
//		printf("]\n");
		cudaMemcpy(tea.d_beta_ij, tea.h_beta_ij, Ntr * sizeof(float), cudaMemcpyHostToDevice);
		checkCUDAError();
	}
}

void deleteTeaUpdater(){
}

