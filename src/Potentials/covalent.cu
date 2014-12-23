/*
 * covalent.cu
 *
 *  Created on: Jun 5, 2009
 *      Author: zhmurov
 */

#include "../gsop.cuh"
#include "covalent.h"
#include "../IO/configreader.h"
#include "../Util/wrapper.h"

#include "covalent_kernel.cu"

/*
 * Covalent potential. Adds itself into the list of potentials and creates timer.
 */
void createCovalentPotential(){
	potentials[potentialsCount] = new CovalentPotential();
	potentialsCount++;
}

/*
 * Construct map of covalent bonds, copy it to the device memory
 */
CovalentPotential::CovalentPotential(){
	this->name = "Covalent";

	printf("Building a map of covalent bounds...\n");

	// Reading parameters from configuration file
	this->kspring_cov = getFloatParameter(COVALENT_KS_STRING, DEFAULT_COVALENT_KS, 1);
	this->R_limit = getFloatParameter(COVALENT_R_LIMIT_STRING, DEFAULT_COVALENT_R_LIMIT, 1);
	this->R_limit_sq = this->R_limit*this->R_limit;
	this->max_covalent = getIntegerParameter(MAX_COVALENT_STRING, DEFAULT_MAX_COVALENT, 1);
	this->blockSize = getIntegerParameter(COVALENT_BLOCK_SIZE_STRING, gsop.blockSize, 1);
	this->blockNum = gsop.aminoCount/this->blockSize + 1;

	// Allocating memory
	this->h_bonds = (GCovalentBond*)calloc(gsop.aminoCount*this->max_covalent, sizeof(GCovalentBond));
	this->h_covalentCount = (int*)calloc(gsop.aminoCount, sizeof(int));
	cudaMalloc((void**)&this->d_bonds, gsop.aminoCount*this->max_covalent*sizeof(GCovalentBond));
	cudaMalloc((void**)&this->d_covalentCount, gsop.aminoCount*sizeof(int));
	this->h_energies = (float*)calloc(gsop.aminoCount, sizeof(float));
	cudaMalloc((void**)&this->d_energies, gsop.aminoCount*sizeof(float));
	this->energies = (float*)calloc(gsop.Ntr, sizeof(float));
	checkCUDAError();

	this->buildMap();

	// Copying data to the device
	cudaMemcpy(this->d_bonds, this->h_bonds, this->max_covalent*gsop.aminoCount*sizeof(GCovalentBond), cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_covalentCount, this->h_covalentCount, gsop.aminoCount*sizeof(int), cudaMemcpyHostToDevice);

	this->updateParametersOnGPU();

	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(covalent_kernel, cudaFuncCachePreferL1);
		cudaFuncSetCacheConfig(covalentEnergy_kernel, cudaFuncCachePreferL1);
	}
}

/*
 * Constructing map of covalent bonds from sop-model
 */
void CovalentPotential::buildMap(){
	int totalCovalent = 0;
	int i, j, k;
	printf("Constructing map of covalent bonds.\n");
	for(k = 0; k < sop.bondCount; k++){
		i = sop.bonds[k].i;
		j = sop.bonds[k].j;
		this->h_bonds[this->h_covalentCount[i]*gsop.aminoCount + i].i2 = j;
		this->h_bonds[this->h_covalentCount[j]*gsop.aminoCount + j].i2 = i;
		this->h_bonds[this->h_covalentCount[i]*gsop.aminoCount + i].r0 = sop.bonds[k].r0;
		this->h_bonds[this->h_covalentCount[j]*gsop.aminoCount + j].r0 = sop.bonds[k].r0;
		this->h_covalentCount[i] ++;
		this->h_covalentCount[j] ++;
		if(this->h_covalentCount[i] > this->max_covalent || this->h_covalentCount[j] > this->max_covalent){
			DIE("Maximum number of covalent bonds exceeded the limit of %d.", this->max_covalent);
		}
		totalCovalent++;
	}

	// Multiplying data for many trajectories (Many-runs-per-GPU)
	int traj;
	for(traj = 1; traj < gsop.Ntr; traj++){
		for(i = 0; i < sop.aminoCount; i++){
			for(k = 0; k < this->max_covalent; k++){
				this->h_bonds[traj*sop.aminoCount + i + k*gsop.aminoCount].i2 =
							this->h_bonds[i + k*gsop.aminoCount].i2 + traj*sop.aminoCount;
				this->h_bonds[traj*sop.aminoCount + i + k*gsop.aminoCount].r0 =
							this->h_bonds[i + k*gsop.aminoCount].r0;
			}
			this->h_covalentCount[traj*sop.aminoCount + i] = this->h_covalentCount[i];
		}
	}

	#ifdef DEBUG
	printf("Covalent bounds (number of bounds, #s of beads, zeros):\n");
	for(int i = 0; i < gsop.aminoCount; i++){
		printf("%d(%d): ", i, gsop.h_covalent[i*width + 0]);
		for(int j = 0; j < max_covalent; j++){
			printf("%d(%3.2f) ", gsop.h_covalent[i*width + j + 1], gsop.h_covalentR0[i*max_covalent + j]);
		}
		printf("\n");
	}
	#endif
	printf("Total number of covalent bonds: %d\n", totalCovalent);
}

/*
 * Organizing data on GPU
 */
void CovalentPotential::updateParametersOnGPU(){
	// All values are organized into a constant for easy access
	hc_covalent.d_bonds = this->d_bonds;
	hc_covalent.d_covalentCount = this->d_covalentCount;
	hc_covalent.d_energies = this->d_energies;
	hc_covalent.kspring_cov = this->kspring_cov;
	hc_covalent.R_limit_sq = this->R_limit_sq;
	cudaMemcpyToSymbol(c_covalent, &hc_covalent, sizeof(CovalentConstant), 0, cudaMemcpyHostToDevice);
	checkCUDAError();
}

/*
 * Call of the computational kernel
 */
void CovalentPotential::compute(){
	covalent_kernel<<<this->blockNum, this->blockSize>>>();
	checkCUDAError();
}

/*
 * Number of energy outputs (Covalent potential returns only FENE energy)
 */
int CovalentPotential::getEnergiesCount(){
	return 1;
}
/*
 * Compute energy for .dat output
 */
float* CovalentPotential::computeEnergy(int id){
	if(id == 0){
		covalentEnergy_kernel<<<this->blockNum, this->blockSize>>>();
		cudaMemcpy(this->h_energies, this->d_energies, gsop.aminoCount*sizeof(float), cudaMemcpyDeviceToHost);
		SOPPotential::sumEnergies(this->h_energies, this->energies);
		checkCUDAError();
		return this->energies;
	} else {
		DIE("Only one energy term is computed by covalent potential");
		return NULL;
	}
}

/*
 * Get energy for a specific trajectory. Should be called after energy is computed.
 */

float CovalentPotential::getEnergy(int traj, int id){
	if(traj < gsop.Ntr && id == 0){
		return this->energies[traj];
	} else {
		DIE("Either trajectory or energy index is out of boundary");
		return 0.0f;
	}
}

