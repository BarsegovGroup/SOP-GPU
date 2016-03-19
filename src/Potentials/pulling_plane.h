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
PARAMETER(pullingPlaneKsCant, float, 0.05f, "?", "...")
PARAMETER(pullingPlaneZeta, float, 200.0f, "?", "...")

PARAMETER_MANDATORY(plane_fixed, std::vector<int>, "", "...")
PARAMETER_MANDATORY(plane_pulled, std::vector<int>, "", "...")

PARAMETER_MANDATORY(pullingPlaneDir, float3, "?", "...")
PARAMETER_MANDATORY(pullingPlanePos, float3, "?", "...")

PARAMETER(pullingPlaneOutput, std::string, "pullplane.<name>_<run>.dat", "path", "...")
PARAMETER(pullingPlaneUpdateFreq, int, parameters::nav, "steps", "...")
PARAMETER(pullingPlaneOutputFreq, int, parameters::nav, "steps", "...")

class PullingPlanePotential : public SOPPotential {
public:
	PullingPlanePotential();
	virtual ~PullingPlanePotential() { }
	virtual void compute();

	int outputWidth;
	int printRuns;

	float Ks;
	float KsCant;
	float deltax;
	float zeta;

	std::vector<int> fixed;
	std::vector<int> pulled;

	float extForce; // Force acting on plane from particles
	float cantForce; // Force acting on cantilever from plane

	float3 pullVector; // Pulling direction (perpendicular to plane)
	float* h_planeDispl; // plane displacement from (0,0,0)
	float* d_planeDispl; // plane displacement from (0,0,0)
	float cantDispl; // cantilever displacement

	float3 planeCoord0; // Initial plane coordinate

	float* h_extForces;
	float* d_extForces;

	int* h_masks;  // 1 for fixed, 2 for pulled beads
	int* d_masks;

	float* h_r0; // Initial distances from the planes
	float* d_r0;

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
	int outputFreq;
};

void createPullingPlanePotential();

