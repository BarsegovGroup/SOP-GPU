/*
 * pairs.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "pairs.h"


#include "pairs_kernel.cu"

void createPairsPotential(){
	potentials[potentialsCount] = new PairsPotential();
	potentialsCount++;
}

PairsPotential::PairsPotential(){
    this->name = "Long Range";
	printf("Initializing pairlist...");
	this->blockSize = gsop.blockSize;
	this->blockNum = gsop.aminoCount/this->blockSize + 1;
	this->a = parameters::a.get();
	this->a2 = this->a*this->a;
	this->el = parameters::el.get();
	this->minus6elovera2 = -6.0f*this->el/this->a2;
	this->max_pairs = parameters::max_pairs.get();
	this->pairsCutoff = parameters::pairs_cutoff.get();
	this->pairsCutoff2 = this->pairsCutoff*this->pairsCutoff;

	// Allocating memory
	this->h_pairs = (int*)calloc(gsop.aminoCount*this->max_pairs, sizeof(int));
	this->h_pairsCount = (int*)calloc(gsop.aminoCount, sizeof(int));
	cudaMalloc((void**)&this->d_pairs, gsop.aminoCount*this->max_pairs*sizeof(int));
	checkCUDAError();
	cudaMalloc((void**)&this->d_pairsCount, gsop.aminoCount*sizeof(int));
    checkCUDAError();
    this->h_energies = (float*)calloc(gsop.aminoCount, sizeof(float));
    cudaMalloc((void**)&this->d_energies, gsop.aminoCount*sizeof(float));
    this->energies = (float*)calloc(gsop.Ntr, sizeof(float));

    this->updateParametersOnGPU();
	printf("done.\n");

	if(deviceProp.major == 2){ // TODO: >= 2 // TODO: do we really need it?
		cudaFuncSetCacheConfig(pairs_kernel, cudaFuncCachePreferL1);
		cudaFuncSetCacheConfig(pairsEnergy_kernel, cudaFuncCachePreferL1);
	}
}

void PairsPotential::updateParametersOnGPU(){
    hc_pairs.pairsCutoff2 = this->pairsCutoff2;
    hc_pairs.a2 = this->a2;
    hc_pairs.el = this->el;
    hc_pairs.minus6elovera2 = this->minus6elovera2;
    hc_pairs.d_pairs = this->d_pairs;
    hc_pairs.d_pairsCount = this->d_pairsCount;
    hc_pairs.max_pairs = this->max_pairs;
    hc_pairs.d_energies = this->d_energies;
	cudaMemcpyToSymbol(c_pairs, &hc_pairs, sizeof(PairsConstant), 0, cudaMemcpyHostToDevice);
	checkCUDAError();
}

void PairsPotential::compute(){
	pairs_kernel<<<this->blockNum, this->blockSize>>>();
	checkCUDAError();
}

int PairsPotential::getEnergiesCount(){
	return 1;
}

float* PairsPotential::computeEnergy(int id){
	if(id == 0){
		pairsEnergy_kernel<<<this->blockNum, this->blockSize>>>();
		cudaMemcpy(this->h_energies, this->d_energies, gsop.aminoCount*sizeof(float), cudaMemcpyDeviceToHost);
		SOPPotential::sumEnergies(this->h_energies, this->energies);
		checkCUDAError();
		return this->energies;
	} else {
		DIE("Pairs potential returns only one energy term");
		return NULL;
	}
}

float PairsPotential::getEnergy(int traj, int id){
	if(traj < gsop.Ntr && id == 0){
		return this->energies[traj];
	} else {
		DIE("Either trajectory or energy index is out of bounds");
		return 0.0f;
	}
}

