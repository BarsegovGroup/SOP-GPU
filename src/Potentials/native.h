/*
 * native.cuh
 *
 *  Created on: Mar 11, 2010
 *      Author: zhmurov
 */
#pragma once

#define NATIVE_R_LIMIT_BOND_STRING		"R_limit_bond"
#define NATIVE_DESOLVATION_STRING		"desolvation"
#define NATIVE_R_WATER_STRING			"rWater"
#define MAX_NATIVE_STRING				"max_native"
#define NATIVE_BLOCK_SIZE_STRING		"block_size_native"

#define DEFAULT_NATIVE_R_LIMIT_BOND			8.0f
#define DEFAULT_NATIVE_DESOLVATION			0
#define DEFAULT_NATIVE_R_WATER				3.0f
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
	virtual void computeEnergy();

    float R_limit_bond; // Used by Updaters/output_manager.cpp
private:
    void buildMap();
    void updateParametersOnGPU();
	int max_native;

    // TODO: next two parameters are completely unused
	int desolvation;
	float rWater;

	int* h_native;   // Map of native interactions
	int* h_nativeCount;
	GNativeParameters* h_nativeParameters;

	int* d_native;
	int* d_nativeCount;
	GNativeParameters* d_nativeParameters;

	int totalNative;

	int blockSize;
	int blockNum;

};

