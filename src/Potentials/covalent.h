/*
 * covalent.cuh
 *
 *  Created on: Mar 11, 2010
 *      Author: zhmurov
 */

#ifndef COVALENT_CUH_
#define COVALENT_CUH_

#define COVALENT_KS_STRING 			"kspring_cov"
#define COVALENT_R_LIMIT_STRING		"R_limit"
#define COVALENT_BLOCK_SIZE_STRING	"block_size_covalent"
#define MAX_COVALENT_STRING			"max_covalent" //Maximum number of covalent bonds per particle

#define DEFAULT_COVALENT_KS			20.0f // Default value for Ks (in kcal/mol)
#define DEFAULT_COVALENT_R_LIMIT	2.0f // Default value for Rc (in A)
#define DEFAULT_MAX_COVALENT		8 // Default maximum of covalent bonds per particle

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
    void buildMap();
    virtual void compute();
	virtual int getEnergiesCount();
	virtual float* computeEnergy(int id);
	virtual float getEnergy(int traj, int id);

	float R_limit; // FENE parameter // Used by Updaters/output_manager.cpp
private:
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

#endif /* COVALENT_CUH_ */
