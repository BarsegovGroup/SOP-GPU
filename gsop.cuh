#ifndef GSOP_CUH
#define GSOP_CUH
/*
 * gsop.cuh
 *
 *  Created on: Jun 4, 2009
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "def_param.h"

#define max_potentials 10
#define max_updaters 10
#define max_timers 20

//#define NOTEXTURE


//#define capsid

struct GSOP{

	cudaDeviceProp deviceProp;

	int pullingOn;
	int minimizationOn;
	int aminoCount;
	int width;

	float4* h_coord;
	float4* d_coord;
	//float4* d_coordToSave;

	float4* h_energies;   // Energies for output (see printStep())
	float4* d_energies;
	//float4* d_energiesToSave;

	float4* h_forces;
	float4* d_forces;
};

struct SOPPotential{
	char name[100];
	unsigned int timer;
	void (*compute)();
	void (*computeEnergy)();
	void (*destroy)();
};

struct SOPUpdater{
	char name[100];
	unsigned int timer;
	int frequency;
	void (*update)();
	void (*destroy)();
};

struct SOPIntegrator{
	char name[100];
	unsigned int timer;
	float h;
	void (*integrate)();
	void (*destroy)();
};

struct SOPTimer{
	char name[100];
	unsigned int timer;
};

extern float3* ext_force;

extern SOPPotential** potentials;
extern SOPUpdater** updaters;
extern SOPIntegrator* integrator;
extern SOPTimer** timers;

extern int potentialsCount;
extern int updatersCount;
extern int timersCount;

extern GSOP gsop;
__device__ __constant__ GSOP c_gsop;

#ifndef NOTEXTURE
texture<float4, 1, cudaReadModeElementType> t_coord; // Coordinates
#endif

extern void checkCUDAError();
extern void copyCoordDeviceToHost();

#endif
