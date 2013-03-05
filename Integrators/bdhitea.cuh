/*
 * bdhitea.cuh
 *
 *  Created on: Mar 4, 2013
 *	  Author: alekseenko
 */

#ifndef BDHITEA_CUH_
#define BDHITEA_CUH_
#include "../gsop.cuh"

#define BDHITEA_EPSILONFREQ_STRING "tea_epsilon_freq"
#define BDHITEA_ON_STRING "tea_on"
#define BDHITEA_A "tea_a"


struct Tea {
	float4 *rforce; // Precomputed random forces
//	float epsilon; // Relative coupling, $\epsilon$ // TODO: probably, unnecessary: epsilon_tmp[0] could be used as well; or we might cache its value on host...
	float *d_beta_ij; // $\beta_{ij}$ from eq. (14) in Geyer&Winter, 2009; its value is not copied to device through this structure
	float *h_beta_ij;
	float *d_epsilon; // Array for reduction
	float *h_epsilon; // Array for reduction
	float3 *d_ci; // Does not actually store `C_i`, only \sum (Dij/Dii)^2
	int epsilon_freq;
	int ntr; // Number of trajectories
	float a; // Bead hydrodynamic radius
};

Tea tea;
__device__ __constant__ Tea c_tea;

SOPIntegrator teaIntegrator;
SOPUpdater teaUpdater;

void createTeaIntegrator();
void initTeaIntegrator();
void integrateTea();
void deleteTeaIntegrator();

void createTeaUpdater();
void initTeaUpdater();
void updateTea();
void deleteTeaUpdater();

#endif /* BDHITEA_CUH_ */
