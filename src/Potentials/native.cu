/*
 * native.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "../Util/wrapper.h"
#include "native.h"
//#define DEBUG1

#include "native_kernel.cu"

void createNativePotential(){
	potentials[potentialsCount] = new NativePotential();
	potentialsCount++;
}

NativePotential::NativePotential(){
	this->name = "Native";

	printf("Building map of native contacts...\n");

	this->R_limit_bond = getFloatParameter(NATIVE_R_LIMIT_BOND_STRING, DEFAULT_NATIVE_R_LIMIT_BOND, 1);
	this->desolvation = getYesNoParameter(NATIVE_DESOLVATION_STRING, DEFAULT_NATIVE_DESOLVATION, 1);
	if(this->desolvation){
		this->rWater = getFloatParameter(NATIVE_R_WATER_STRING, DEFAULT_NATIVE_R_WATER, 1);;
	}
	this->max_native = getIntegerParameter(MAX_NATIVE_STRING, DEFAULT_MAX_NATIVE, 1);
	this->blockSize = getIntegerParameter(NATIVE_BLOCK_SIZE_STRING, gsop.blockSize, 1);
	this->blockNum = gsop.aminoCount/this->blockSize + 1;


	this->h_native = (int*)calloc(gsop.aminoCount*this->max_native, sizeof(int));
	cudaMalloc((void**)&this->d_native, gsop.aminoCount*this->max_native*sizeof(int));
	this->h_nativeCount = (int*)calloc(gsop.aminoCount, sizeof(int));
	cudaMalloc((void**)&this->d_nativeCount, gsop.aminoCount*sizeof(int));
	this->h_nativeParameters = (GNativeParameters*)calloc(this->max_native*gsop.aminoCount, sizeof(GNativeParameters));
	cudaMalloc((void**)&this->d_nativeParameters, this->max_native*gsop.aminoCount*sizeof(GNativeParameters));
	// Building map
	int totalNative = 0;
	int i, j, k;
	for(k = 0; k < sop.nativeCount; k++){
		i = sop.natives[k].i;
		j = sop.natives[k].j;
		this->h_native[this->h_nativeCount[i]*gsop.aminoCount + i] = j;
		this->h_native[this->h_nativeCount[j]*gsop.aminoCount + j] = i;
		this->h_nativeParameters[this->h_nativeCount[i]*gsop.aminoCount + i].r02 = sop.natives[k].r0*sop.natives[k].r0;
		this->h_nativeParameters[this->h_nativeCount[j]*gsop.aminoCount + j].r02 = sop.natives[k].r0*sop.natives[k].r0;
		this->h_nativeParameters[this->h_nativeCount[i]*gsop.aminoCount + i].minus12ehOverR02 =
														-12.0*sop.natives[k].eh/(sop.natives[k].r0*sop.natives[k].r0);
		this->h_nativeParameters[this->h_nativeCount[j]*gsop.aminoCount + j].minus12ehOverR02 =
														-12.0*sop.natives[k].eh/(sop.natives[k].r0*sop.natives[k].r0);
		this->h_nativeCount[i] ++;
		this->h_nativeCount[j] ++;
		if(this->h_nativeCount[i] > this->max_native || this->h_nativeCount[j] > this->max_native){
			DIE("ERROR: Maximum number of native contacts exceeded the limit of %d.", this->max_native);
		}
		totalNative++;
	}

	for(j = 1; j < gsop.Ntr; j++){
		for(i = 0; i < sop.aminoCount; i++){
			for(k = 0; k < this->max_native; k++){
				this->h_native[j*sop.aminoCount + i + k*gsop.aminoCount] =
							this->h_native[i + k*gsop.aminoCount] + j*sop.aminoCount;
				this->h_nativeParameters[j*sop.aminoCount + i + k*gsop.aminoCount] =
							this->h_nativeParameters[i + k*gsop.aminoCount];
			}
			this->h_nativeCount[j*sop.aminoCount + i] = this->h_nativeCount[i];
		}
	}


	#ifdef DEBUG1
	printf("Native contacts (number of contacts, #s of beads, zeros):\n");
	for(int i = 0; i < gsop.aminoCount; i++){
		printf("%d (%d): ", i, this->h_nativeCount[i]);
		for(int j = 0; j < this->h_nativeCount[i]; j++){
			printf("%d(%3.1f,%3.1f)  ", this->h_native[j*gsop.aminoCount + i],
					this->h_nativeParameters[j*gsop.aminoCount + i].x,
					this->h_nativeParameters[j*gsop.aminoCount + i].y);
		}
		printf("\n");
	}
	#endif
	// Copying to device
	cudaMemcpy(this->d_native, this->h_native, this->max_native*gsop.aminoCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_nativeCount, this->h_nativeCount, gsop.aminoCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_nativeParameters, this->h_nativeParameters, this->max_native*gsop.aminoCount*sizeof(float2), cudaMemcpyHostToDevice);

    hc_native.d_native = this->d_native;
    hc_native.d_nativeCount = this->d_nativeCount;
    hc_native.d_nativeParameters = this->d_nativeParameters;
	cudaMemcpyToSymbol(c_native, &hc_native, sizeof(NativeConstant), 0, cudaMemcpyHostToDevice);

	printf("Total number of native contacts: %d \n", totalNative);

	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(native_kernel, cudaFuncCachePreferL1);
		cudaFuncSetCacheConfig(nativeEnergy_kernel, cudaFuncCachePreferL1);
	}
}

void NativePotential::compute(){
	native_kernel<<<this->blockNum, this->blockSize>>>();
	checkCUDAError();
}

void NativePotential::computeEnergy(){
	nativeEnergy_kernel<<<this->blockNum, this->blockSize>>>();
	checkCUDAError();
}

