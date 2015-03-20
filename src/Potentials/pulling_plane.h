/*
 * pulling_plane.cuh
 *
 *  Created on: Jan 17, 2014
 *      Author: alekseenko
 */

#pragma once

#include "../Util/parameters.h"

PARAMETER(pullingPlane, bool, false, "true/false", "...")

PARAMETER(pullingPlaneDeltax, float, 0.0f, "?", "...")
PARAMETER(pullingPlaneKs, float, 0.05f, "?", "...")

PARAMETER_MANDATORY(plane_fixed, std::vector<int>, "", "...")
PARAMETER_MANDATORY(plane_pulled, std::vector<int>, "", "...")

PARAMETER_MANDATORY(pullingPlaneDir, float3, "?", "...")
PARAMETER_MANDATORY(pullingPlanePos, float3, "?", "...")

PARAMETER(pullingPlaneOutput, std::string, "pullplane.<name>_<run>.dat", "path", "...")
PARAMETER(pullingPlaneFreq, int, parameters::nav, "steps", "...")

class PullingPlanePotential : public SOPPotential {
public:
    PullingPlanePotential();
    virtual ~PullingPlanePotential() { }
    virtual void compute();

	float Ks;
	float deltax;
    float mass;

    std::vector<int> fixed;
    std::vector<int> pulled;

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

