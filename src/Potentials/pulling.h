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
#define PULLING_FIXED_COUNT_STRING	"fixed_beads"
#define PULLING_FIXED_STRING		"fixed"
#define PULLING_PULLED_COUNT_STRING	"pulled_beads"
#define PULLING_PULLED_STRING		"pulled"
PARAMETER(pullDirection, std::string, "endToEnd", "endToEnd/vector", "...")
PARAMETER_MANDATORY(pullVector, float3, "?", "...")
PARAMETER_MANDATORY(fixedEnd, int, "bead ID", "...")
PARAMETER_MANDATORY(pulledEnd, int, "bead ID", "...")
PARAMETER(pullOutput, std::string, "pull.<name>_<author><run>_<stage>.dat", "path", "...")
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
	int fixedCount;
	int pulledCount;
	int* fixed;
	int* pulled;
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

