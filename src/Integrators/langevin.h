/*
 * langevin.cuh
 *
 *  Created on: Apr 27, 2010
 *      Author: zhmurov
 */

#pragma once
#include "../gsop.h"
#include "../Util/parameters.h"

PARAMETER_MANDATORY(timestep, float, "?; truly ?", "...")
PARAMETER(zeta, float, 50.0f, "?", "...")
PARAMETER(temperature, float, 0.6, "?", "...")

PARAMETER(minimization, bool, false, "true/false", "...")

PARAMETER(heating, bool, false, "true/false", "...")
PARAMETER_MANDATORY(tempFreq, int, "steps", "...")
PARAMETER_MANDATORY(initialT, float, "?", "...")
PARAMETER_MANDATORY(deltaT, float, "?", "...")

PARAMETER_LAZY(seed, int, time(NULL), "", "...")

class LangevinIntegrator : public SOPIntegrator {
public:
    LangevinIntegrator();
    virtual ~LangevinIntegrator();
    virtual void integrate();
    void setTemperature(float temp);
	
    float hOverZeta; // Accessed in Potentials/pulling_plane.cu
private:
    void updateParametersOnGPU();
	float zeta;
	float var;
	float tempNorm;
	float temp;
};

class TemperatureUpdater : public SOPUpdater {
public:
    TemperatureUpdater(LangevinIntegrator *langevin);
    virtual ~TemperatureUpdater();
    virtual void update();
private:
    LangevinIntegrator *langevin;
	float initialT;
	float deltaT;
};

void createLangevinIntegrator();

