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
#define MAX_COVALENT_STRING			"max_covalent"

#define DEFAULT_COVALENT_KS			20.0f
#define DEFAULT_COVALENT_R_LIMIT	2.0f
#define DEFAULT_MAX_COVALENT		4

/*
 * Data structure for covalent bond on GPU (int for map + float for r0)
 */
struct __align__(8) GCovalentBond{
	int i2;
	float r0;
};

/*
 * All data neded to compute covalent potential
 */

void createCovalentPotential();

class CovalentPotential: public SOPPotential{
public:
    CovalentPotential();
    virtual ~CovalentPotential() { }
    virtual void compute();
    virtual void computeEnergy();

	int max_covalent; // Max number of covalent bonds per particle
	float kspring_cov; // Covalent spring constant (FENE)
	float R_limit; // FENE parameter
	float R_limit_sq; // Same, squared

	GCovalentBond* h_bonds; // Map of covalent bonds (host)
	int* h_covalentCount; // Numbers of covalent bonds per particle (host)

	GCovalentBond* d_bonds; // Same on device
	int* d_covalentCount;

	int blockSize; // Block size and number of blocks for covalent bonds kernel evaluation
	int blockNum;
};

#endif /* COVALENT_CUH_ */
