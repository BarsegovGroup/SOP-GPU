/*
 * bdhitea.cuh
 *
 *  Created on: Mar 4, 2013
 *	  Author: alekseenko
 */

#ifndef BDHITEA_CUH_
#define BDHITEA_CUH_
#include "../gsop.cuh"

// TODO: Rename this module, since it not inly implements TEA-HI, but also exact (Cholesky-based) HI

#define BDHITEA_ON_STRING "tea_on" // enable/disable TEA
#define BDHITEA_EXACT_STRING "tea_exact" // use Cholesky-based HI treatment; it's not TEA anymore
#define BDHITEA_A_STRING "tea_a" // set hydrodynamic radii of single bead, Angstroms
#define BDHITEA_EPSILONFREQ_STRING "tea_epsilon_freq" // frequncy of updating ersatz coefficients, steps
#define BDHITEA_CAPRICIOUS_STRING "tea_capricious" // whether to abort execution on weird values of HI tensor
#define BDHITEA_UNLISTED_STRING "tea_unlisted" // whether to calculate all-to-all interactions, or just the ones in pairlist
#define BDHITEA_EPSMAX_STRING "tea_epsmax" // Abort simulation if epsilon reaches this value, unitless

// max. size of matrix for single system when using cholesky integrator
#define CHOLSIZE ((3*128))

// undefine to disable use of textures in TEA
#define TEA_TEXTURE

struct Tea {
	float4 *rforce; // Precomputed random forces
	float4 *mforce; // Copy of molecular forces for each bead
	float4 *coords; // Copy of coordinates for each bead
	float *d_beta_ij; // $\beta_{ij}$ from eq. (14) in Geyer&Winter, 2009;
	float *h_beta_ij;
	float *d_epsilon; // Array of epsilon values for individual beads, evaluated on device and used to fill d_beta_ij
	float *h_epsilon;
	float4 *d_ci; // Does not actually store $C_i$, only $\sum (Dij/Dii)^2$
    float *d_tensor; // Memory for tensor when running in `exact` mode
	int epsilon_freq; // How often to update epsilon, beta_ij and c_i
	int namino; // Number of aminos per trajectory, nothing special
	float a; // Bead hydrodynamic radius, in A
	int capricious; // If != 0, then the simulation will stop if HI tensor has abnormal values. If zero, the simulation will continue anyway (and it is probably perfectly fine).
	int unlisted; // If ==0, then beads will interact hydrodynamically with their friends in covalent, native and pairs lists
    int exact; // If > 0, use Cholesky-based treatment
	float epsmax; // If epsilon exceeds this value, then abort simulation. Default: never [in capricious mode, epsilon > 1.0 will trigger stop anyway]
};

Tea tea;
__device__ __constant__ Tea c_tea;

#ifdef TEA_TEXTURE
texture<float4, 1, cudaReadModeElementType> t_rforce;
texture<float4, 1, cudaReadModeElementType> t_mforce;
#endif


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

// SOP-GPU with TEA --- Standard Operating Procedures for Generator Protection Unit with Technical and Economic Analysis

