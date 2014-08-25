#pragma once
/*
 * gsop.h
 *
 *  Created on: Jun 4, 2009
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector_types.h>
#include <cuda.h>
#include "IO/topio.h"

#define max_potentials 10
#define max_updaters 10
#define max_timers 20

//#define NOTEXTURE
const int DEFAULT_BLOCK_SIZE = 256;


//#define capsid

struct GSOP{
    int deviceId;
    int blockSize;

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

	float4* h_energies;   // Energies for output (see printStep())
	float4* d_energies;

	float4* h_forces;
	float4* d_forces;
};

class SOPPotential{
public:
	virtual ~SOPPotential() { }
    std::string name;

	virtual void compute() = 0;
	virtual void computeEnergy() = 0;
};

class SOPUpdater{
public:
	virtual ~SOPUpdater() { }
    std::string name;
	int frequency;

	virtual void update() = 0;
};

class SOPIntegrator{
public:
	virtual ~SOPIntegrator() { }
    std::string name;
	float h;

	virtual void integrate() = 0;
};

extern SOPPotential** potentials;
extern SOPUpdater** updaters;
extern SOPIntegrator* integrator;

extern int potentialsCount;
extern int updatersCount;
extern int integratorTea;

extern GSOP gsop;
extern SOP sop;

SOPPotential* potentialByName(const char* name);
SOPUpdater* updaterByName(const char* name);

extern cudaDeviceProp deviceProp;

void initGPU();
void initFF();
void runGPU();
void bindTextures();
void __checkCUDAError(const char *file, int line);

void initCoordinates();
void copyCoordinates();
void copyCoordinatesTrajectory(int traj);
void initForces();
void initEnergies();
void copyEnergiesDeviceToHost();

void __checkCUDAError(const char *file, int line);
void copyCoordDeviceToHost();
#define checkCUDAError() __checkCUDAError(__FILE__, __LINE__)

