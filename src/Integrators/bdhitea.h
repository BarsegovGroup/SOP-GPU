/*
 * bdhitea.cuh
 *
 *  Created on: Mar 4, 2013
 *	  Author: alekseenko
 */

#pragma once

#include "langevin.h"
#include "../Util/parameters.h"

// TODO: Rename this module, since it not inly implements TEA-HI, but also exact (Cholesky-based) HI

PARAMETER(tea_on, bool, false, "true/false", "enable/disable TEA")
PARAMETER(tea_exact, bool, false, "true/false", "use Cholesky-based HI treatment; it's not TEA anymore")
PARAMETER(tea_a, float, 1.8f, "A", "set hydrodynamic radii of single bead")
PARAMETER_MANDATORY(tea_epsilon_freq, int, "steps", "frequncy of updating ersatz coefficients")
PARAMETER(tea_capricious, bool, true, "true/false", "whether to abort execution on weird values of HI tensor")
PARAMETER(tea_unlisted, bool, true, "true/false", "whether to calculate all-to-all interactions, or just the ones in pairlist")
PARAMETER(tea_epsmax, float, 999.f, "", "Abort simulation if epsilon reaches this value and tea_capricious=on; because epsilon will never exceed 1, by default will never trigger")

// max. size of matrix for single system when using cholesky integrator
#define CHOLSIZE ((3*128))

// undefine to disable use of textures in TEA
#define TEA_TEXTURE

struct TeaConstant {
	float4 *rforce; // Precomputed random forces
	float4 *mforce; // Copy of molecular forces for each bead
	float4 *coords; // Copy of coordinates for each bead
	float *d_beta_ij; // $\beta_{ij}$ from eq. (14) in Geyer&Winter, 2009;
	float *d_epsilon; // Array of epsilon values for individual beads, evaluated on device and used to fill d_beta_ij
	float *h_epsilon;
	float4 *d_ci; // Does not actually store $C_i$, only $\sum (Dij/Dii)^2$
    float *d_tensor; // Memory for tensor when running in `exact` mode
	int namino; // Number of aminos per trajectory, nothing special
	float a; // Bead hydrodynamic radius, in A
};

class TeaIntegrator : public LangevinIntegrator{
public:
    TeaIntegrator();
    virtual ~TeaIntegrator();
    virtual void integrate();
	int capricious; // If != 0, then the simulation will stop if HI tensor has abnormal values. If zero, the simulation will continue anyway (and it is probably perfectly fine).
	int unlisted; // If ==0, then beads will interact hydrodynamically with their friends in covalent, native and pairs lists
    int exact; // If > 0, use Cholesky-based treatment
	float epsmax; // If epsilon exceeds this value, then abort simulation. Default: never [in capricious mode, epsilon > 1.0 will trigger stop anyway]
	int epsilon_freq; // How often to update epsilon, beta_ij and c_i

	float *h_epsilon;
	float *h_beta_ij;
};

class TeaUpdater : public SOPUpdater{
public:
    TeaUpdater(TeaIntegrator *tea);
    virtual ~TeaUpdater();
    virtual void update();
private:
    TeaIntegrator *tea;
};

void createTeaIntegrator();

// SOP-GPU with TEA --- Standard Operating Procedures for Generator Protection Unit with Technical and Economic Analysis

