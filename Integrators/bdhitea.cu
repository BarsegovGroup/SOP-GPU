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
	printf("TEA integrator is untested and should not be used!\n");
	exit(-1);
	cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 655350000);
	checkCUDAError();
	sprintf(teaIntegrator.name, "BDHI-TEA");
	teaIntegrator.integrate = &integrateTea;
	teaIntegrator.destroy = &deleteTeaIntegrator;
	integrator = &teaIntegrator;
	initTeaIntegrator();
}

void initTeaIntegrator(){
	if (gsop.minimizationOn) {
		printf("Using minimization with BDHI integrator is not supported and will never be. I'm so sorry.");
		exit(-1);
	}
	initLangevinIntegrator();

	tea.ntr = gsop.aminoCount / sop.aminoCount;
	tea.epsilon_freq = getIntegerParameter(BDHITEA_EPSILONFREQ_STRING, 0, 0);
	tea.a = getFloatParameter(BDHITEA_A, 3.8f, 1);

	cudaMalloc(&tea.rforce, gsop.aminoCount * sizeof(float4));
	cudaMalloc(&tea.d_epsilon, gsop.aminoCount * sizeof(float));
	cudaMalloc(&tea.d_ci, gsop.aminoCount * sizeof(float3));
	cudaMalloc(&tea.d_beta_ij, Ntr * sizeof(float));
	checkCUDAError();
	tea.h_epsilon = (float*) malloc(gsop.aminoCount * sizeof(float));
	tea.h_beta_ij = (float*) malloc(tea.ntr * sizeof(float));
	cudaMemcpyToSymbol(c_tea, &tea, sizeof(Tea), 0, cudaMemcpyHostToDevice);
	checkCUDAError();

	createTeaUpdater();
	printf("TEA integrator initialized, a = %.2f A, freq = %d steps; total %d trajectories\n", tea.a, tea.epsilon_freq, tea.ntr);
}

void integrateTea(){
	// Pregenerate random forces
	integrateTea_prepare<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
	cudaDeviceSynchronize();
	checkCUDAError();
	// Integrate
	integrateTea_kernel<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
	cudaDeviceSynchronize();
	checkCUDAError();
	// Clean up forces
	integrateTea_cleanup<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
	cudaDeviceSynchronize();
	checkCUDAError();
}

void deleteTeaIntegrator(){
	deleteLangevinIntegrator();
	cudaFree(tea.rforce);
	cudaFree(tea.d_epsilon);
	cudaFree(tea.d_beta_ij);
	free(tea.h_epsilon);
	free(tea.h_beta_ij);
}


void createTeaUpdater(){
	sprintf(teaUpdater.name, "BDHI-TEA updater");
	teaUpdater.update = &updateTea;
	teaUpdater.destroy = &deleteTeaUpdater;
	teaUpdater.frequency = getIntegerParameter(BDHITEA_EPSILONFREQ_STRING, 0, 0);
	updaters[updatersCount] = &teaUpdater;
	updatersCount++;
	initTeaUpdater();
}

void initTeaUpdater(){

}

void updateTea(){
	int update_epsilon = (step % tea.epsilon_freq) == 0;
	if (update_epsilon){
		//printf("updating epsilon @ step %ld\n",step);
		// Calculate final relative coupling
		integrateTea_epsilon<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
		checkCUDAError();
		// Do host-based many-runs-per-gpu-compatible reduction
		cudaMemcpy(tea.h_epsilon, tea.d_epsilon, gsop.aminoCount * sizeof(float), cudaMemcpyDeviceToHost);
		checkCUDAError();
		const int N = sop.aminoCount;
		printf("Beta_ij: [ ");
		for (int t = 0; t < tea.ntr; ++t){
			float epsilon = 0.0;
			for (int i = 0; i < N; ++i){
				epsilon += tea.h_epsilon[t*sop.aminoCount + i];
			}
			epsilon += N*3; // Diagonal elements;
			epsilon /= 3*N*3*N; // Averaging
			tea.h_beta_ij[t] = (1.f - sqrtf(1.f - (3*N-1.0f)*epsilon*epsilon + (3*N-2.f)*epsilon)) / ((3*N-1.f)*epsilon*epsilon - (3*N-2.f)*epsilon); // eq. (26)
			printf("%.3f ", tea.h_beta_ij[t]);
		}
		printf("]\n");
		cudaMemcpy(tea.d_beta_ij, tea.h_beta_ij, Ntr * sizeof(float), cudaMemcpyHostToDevice);
		checkCUDAError();
	}
}

void deleteTeaUpdater(){

}


