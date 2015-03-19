/*
 * covalent.cuh
 *
 *  Created on: Mar 11, 2010
 *      Author: zhmurov
 */

#pragma once

#include "../Util/parameters.h"

PARAMETER(kspring_cov, float, 20.0f, "kcal/mol", "Covalent spring constant (FENE)")
PARAMETER(R_limit, float, 2.0f, "A", "FENE parameter")
PARAMETER(max_covalent, int, 8, "", "Max number of covalent bonds per particle")

/*
 * Data structure for covalent bond on GPU (int for map + float for r0)
 */
struct __align__(8) GCovalentBond{
	int i2; // Index of the second particle
	float r0; // Equilibrium distance
};

/*
 * All data neded to compute covalent potential
 */

void createCovalentPotential();

class CovalentPotential: public SOPPotential{
public:
    CovalentPotential();
    virtual ~CovalentPotential() { }
    void updateParametersOnGPU();
    virtual void compute();
	virtual int getEnergiesCount() const;
	virtual float* computeEnergy(int id);
	virtual float getEnergy(int traj, int id);

	float R_limit; // FENE parameter // Used by Updaters/output_manager.cpp
private:
    void buildMap();

	int max_covalent; // Max number of covalent bonds per particle
	float kspring_cov; // Covalent spring constant (FENE)
	float R_limit_sq; // Same, squared

	GCovalentBond* h_bonds; // Map of covalent bonds (host)
	int* h_covalentCount; // Numbers of covalent bonds per particle (host)

	GCovalentBond* d_bonds; // Same on device
	int* d_covalentCount;

	float* h_energies; // Energies per particle (on CPU)
	float* d_energies; // Energies per particle (on GPU)
	float* energies; // Energies per trajectory (on CPU)

	int blockSize; // CUDA block size and number of blocks for covalent bonds kernel evaluation
	int blockNum; // Number of CUDA blocks
};

