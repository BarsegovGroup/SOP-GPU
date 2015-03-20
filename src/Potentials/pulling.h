/*
 * pulling.cuh
 *
 *  Created on: May 25, 2010
 *      Author: zhmurov
 */

#pragma once

#include "../Util/parameters.h"

#include <vector>
#include <string>

PARAMETER(pulling, bool, false, "true/false", "...")
PARAMETER(deltax, float, 0.0, "?", "...")
PARAMETER(k_trans, float, 0.05f, "?", "...")
PARAMETER(fconst, float, 0.0, "?", "...")
PARAMETER_MANDATORY(fixed, std::vector<int>, "", "...")
PARAMETER_MANDATORY(pulled, std::vector<int>, "", "...")
PARAMETER(pullDirection, std::string, "endToEnd", "endToEnd/vector", "...")
PARAMETER_MANDATORY(pullVector, float3, "?", "...")
PARAMETER_MANDATORY(fixedEnd, int, "bead ID", "...")
PARAMETER_MANDATORY(pulledEnd, int, "bead ID", "...")
PARAMETER(pullOutput, std::string, "pull.<name>_<run>.dat", "path", "...")
PARAMETER(pullFreq, int, parameters::nav, "steps", "...")

class PullingPotential : public SOPPotential{
public:
    PullingPotential();
    virtual ~PullingPotential() { }
    virtual void compute();
    void updateForces(float xt);
    void savePullingData();
	float deltax; // Used by updater
private:
    void updateParametersOnGPU();
    float Ks;
	float fconst;
	int fixedEnd;
	int pulledEnd;
    std::vector<int> fixed;
    std::vector<int> pulled;
	float3* pullVector;
	float3* extForce;
	float4* h_extForces;
	float4* d_extForces;
	float3* cantCoord0;
	float3* cantCoord;

    int blockSize, blockNum;
    std::vector<std::string> pullFilenames;

    float3 computeForce(const float4 &coordN, int traj) const;
};

class PullingUpdater : public SOPUpdater{
public:
    PullingUpdater(PullingPotential *pulling);
    virtual ~PullingUpdater() { }
    virtual void update();
private:
    PullingPotential *pulling;
};

void createPullingPotential();

