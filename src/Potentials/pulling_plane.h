/*
 * pulling_plane.cuh
 *
 *  Created on: Jan 17, 2014
 *      Author: alekseenko
 */

#pragma once

#define PULLINGPLANE_ON_STRING			"pullingPlane"

#define PULLINGPLANE_DELTAX_STRING		"pullingPlaneDeltax"
#define PULLINGPLANE_KS_STRING			"pullingPlaneKs"

#define PULLINGPLANE_FIXED_COUNT_STRING	"plane_fixed_beads"
#define PULLINGPLANE_FIXED_STRING		"plane_fixed"
#define PULLINGPLANE_PULLED_COUNT_STRING	"plane_pulled_beads"
#define PULLINGPLANE_PULLED_STRING		"plane_pulled"

#define PULLINGPLANE_PULLVECTOR_STRING		"pullingPlaneDir"
#define PULLINGPLANE_ZEROVECTOR_STRING		"pullingPlanePos"

#define PULLINGPLANE_FILENAME			"pullingPlaneOutput"
#define PULLINGPLANE_FREQ				"pullingPlaneFreq"

#define DEFAULT_PULLINGPLANE_KS			0.05f
#define DEFAULT_PULLINGPLANE_FILENAME	"pullplane.<name>_<author><run>_<stage>.dat"

class PullingPlanePotential : public SOPPotential {
public:
    PullingPlanePotential();
    virtual ~PullingPlanePotential() { }
    virtual void compute();
    virtual void computeEnergy();

	float Ks;
	float deltax;
    float mass;

	int fixedCount;
	int pulledCount;
	int* fixed;
	int* pulled;

	float3 extForce; // Force acting on plane from particles
    float3 cantForce; // Force acting on cantilever from plane

	float3 pullVector; // Pulling direction (perpendicular to plane)
    float d; // plane displacement from (0,0,0)
    float d0; // initial plane displacement from (0,0,0)
    float cant_d; // cantilever displacement

	float3 planeCoord0; // Initial plane coordinate
	float3 planeCoord;  // Current (actual) plane coordinate
	float3 cantCoord;   // Position of cantilever

	float4* h_extForces;
	float4* d_extForces;

private:
    int blockNum, blockSize;
};

class PullingPlaneUpdater : public SOPUpdater {
public:
    PullingPlaneUpdater(PullingPlanePotential *pullingPlane);
    virtual ~PullingPlaneUpdater() { }
    virtual void update();
private:
    PullingPlanePotential *pullingPlane;
};

void createPullingPlanePotential();

