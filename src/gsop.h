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
#include "Util/wrapper.h"
#include "Util/parameters.h"

#define max_potentials 10 //TODO : use vector
#define max_updaters 10
#define max_timers 20

#define VERSION "trunk"

//#define NOTEXTURE

PARAMETER_MANDATORY(topology, std::string, "path", "...")
PARAMETER_MANDATORY(coordinates, std::string, "path", "...")
PARAMETER(reffilename, std::string, "<name>.ref.pdb", "path", "...")
PARAMETER(name, std::string, "unnamed", "", "...")
PARAMETER(run, int, -1, "?", "...")
PARAMETER_MANDATORY(runnum, int, "?", "...")
PARAMETER_MANDATORY(firstrun, int, "?", "...")
PARAMETER_MANDATORY(numsteps, long, "?", "...")
PARAMETER(device, int, 0, "", "...")
PARAMETER(nav, int, 1000, "steps", "...")


//#define capsid

struct GSOP{
    int deviceId;
    int blockSize;

	int pullingOn; // TODO: bitfield?
	int pullingPlaneOn;
	int minimizationOn;
	int heatingOn;
	int indentationOn;

	/*
	 * These do not belong here.
	 */
	long int step;
	long int numsteps;
	int nav;

	int aminoCount; // total aminoCount in all trajectories
	int width; // aligned aminoCount
    int Ntr;
    int firstrun;

	float4* h_coord;
	float4* d_coord;

	float* h_T;   // Energies for output (see printStep())
	float* d_T;

	float4* h_forces;
	float4* d_forces;
};

class SOPPotential{
public:
	virtual ~SOPPotential() { }
    std::string name;

	virtual void compute() = 0;

	virtual int getEnergiesCount() const {
		return 0;
	}
	virtual float* computeEnergy(int id){ // TODO: make const
		DIE("No energy computation is implemented for this potential");
	}
	virtual float getEnergy(int traj, int id){
		DIE("No energy computation is implemented for this potential");
	}
protected:
	void sumEnergies(const float *h_energies, float *energies);
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
void initTemperature();
void copyTemperatureDeviceToHost();

void __checkCUDAError(const char *file, int line);
void copyCoordDeviceToHost();
#define checkCUDAError() __checkCUDAError(__FILE__, __LINE__)

