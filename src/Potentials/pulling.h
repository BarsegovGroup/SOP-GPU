/*
 * pulling.cuh
 *
 *  Created on: May 25, 2010
 *      Author: zhmurov
 */

#pragma once

#include <vector>
#include <string>

#define PULLING_ON_STRING			"pulling"

#define PULLING_DELTAX_STRING		"deltax"
#define PULLING_KS_STRING			"k_trans"
#define PULLING_FCONST_STRING		"fconst"
#define PULLING_FIXED_COUNT_STRING	"fixed_beads"
#define PULLING_FIXED_STRING		"fixed"
#define PULLING_PULLED_COUNT_STRING	"pulled_beads"
#define PULLING_PULLED_STRING		"pulled"
#define PULLING_DIRECTION_STRING	"pullDirection"
#define PULLING_VECTOR_STRING		"pullVector"
#define PULLING_FIXED_END_STRING	"fixedEnd"
#define PULLING_PULLED_END_STRING	"pulledEnd"
#define PULLING_FILENAME			"pullOutput"
#define PULLING_FREQ				"pullFreq"

#define PULLING_DIRECTION_ENDTOEND_STRING	"endToEnd"
#define PULLING_DIRECTION_VECTOR_STRING		"vector"

#define DEFAULT_PULLING_KS			0.05f
#define DEFAULT_PULLING_DIRECTION	PULLING_DIRECTION_ENDTOEND_STRING
#define DEFAULT_PULLING_FILENAME	"pull.<name>_<author><run>_<stage>.dat"

class PullingPotential : public SOPPotential{
public:
    PullingPotential();
    virtual ~PullingPotential() { }
    virtual void compute();
    virtual void computeEnergy();
    void updateForces(float xt);
    void savePullingData();
	float deltax; // Used by updater
private:
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

    float3 computeForce(float4 coordN, int traj) const;
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

