/*
 * native.cuh
 *
 *  Created on: Mar 11, 2010
 *      Author: zhmurov
 */
#pragma once

#include "../Util/parameters.h"

PARAMETER(R_limit_bond, float, 8.0f, "A", "...")
PARAMETER(max_native, int, 128, "", "...")
#define NATIVE_R_LIMIT_BOND_STRING		"R_limit_bond"
#define MAX_NATIVE_STRING				"max_native"

#define DEFAULT_NATIVE_R_LIMIT_BOND			8.0f
#define DEFAULT_MAX_NATIVE					128

struct __align__(8) GNativeParameters{
	float r02;
	float minus12ehOverR02;
};

void createNativePotential();

class NativePotential : public SOPPotential{
public:
    NativePotential();
    virtual ~NativePotential() { }
    virtual void compute();
	virtual int getEnergiesCount();
	virtual float* computeEnergy(int id);
	virtual float getEnergy(int traj, int id);

private:
    void buildMap();
    void updateParametersOnGPU();
	int max_native;
    float R_limit_bond;

	int* h_native;   // Map of native interactions
	int* h_nativeCount;
	GNativeParameters* h_nativeParameters;

	int* d_native;
	int* d_nativeCount;
	GNativeParameters* d_nativeParameters;

	float* h_energies;
	float* d_energies;
	float* energies;

	int totalNative;

	int blockSize;
	int blockNum;

};

