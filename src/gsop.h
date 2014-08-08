#ifndef GSOP_H
#define GSOP_H
/*
 * gsop.h
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
    int deviceId;
	cudaDeviceProp deviceProp;

	int pullingOn; // TODO: bitfield?
	int pullingPlaneOn;
	int minimizationOn;
	int heatingOn;
	int indentationOn;

	int aminoCount; // total aminoCount in all trajectories
	int width; // aligned aminoCount
    int Ntr;
    int firstrun;

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

extern SOPPotential** potentials;
extern SOPUpdater** updaters;
extern SOPIntegrator* integrator;

extern int potentialsCount;
extern int updatersCount;
extern int integratorTea;

extern GSOP gsop;

#endif

