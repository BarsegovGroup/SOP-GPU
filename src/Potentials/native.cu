/*
 * native.cu
 *
 *  Created on: Jun 16, 2009
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "native.h"
//#define DEBUG1

Native native;
__device__ __constant__ Native c_native;
SOPPotential nativePotential;

#include "native_kernel.cu"

void createNativePotential(){
	sprintf(nativePotential.name, "Native");
	nativePotential.compute = &computeNative;
	nativePotential.computeEnergy = &computeNativeEnergy;
	potentials[potentialsCount] = &nativePotential;
	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(native_kernel, cudaFuncCachePreferL1);
		cudaFuncSetCacheConfig(nativeEnergy_kernel, cudaFuncCachePreferL1);
	}
	potentialsCount++;
	initNative();
}

void initNative(){

	printf("Building map of native contacts...\n");

	native.R_limit_bond = getFloatParameter(NATIVE_R_LIMIT_BOND_STRING, DEFAULT_NATIVE_R_LIMIT_BOND, 1);
	native.desolvation = getYesNoParameter(NATIVE_DESOLVATION_STRING, DEFAULT_NATIVE_DESOLVATION, 1);
	if(native.desolvation){
		native.rWater = getFloatParameter(NATIVE_R_WATER_STRING, DEFAULT_NATIVE_R_WATER, 1);;
	}
	native.max_native = getIntegerParameter(MAX_NATIVE_STRING, DEFAULT_MAX_NATIVE, 1);
	native.blockSize = getIntegerParameter(NATIVE_BLOCK_SIZE_STRING, gsop.blockSize, 1);
	native.blockNum = gsop.aminoCount/native.blockSize + 1;


	native.h_native = (int*)calloc(gsop.aminoCount*native.max_native, sizeof(int));
	cudaMalloc((void**)&native.d_native, gsop.aminoCount*native.max_native*sizeof(int));
	native.h_nativeCount = (int*)calloc(gsop.aminoCount, sizeof(int));
	cudaMalloc((void**)&native.d_nativeCount, gsop.aminoCount*sizeof(int));
	native.h_nativeParameters = (GNativeParameters*)calloc(native.max_native*gsop.aminoCount, sizeof(GNativeParameters));
	cudaMalloc((void**)&native.d_nativeParameters, native.max_native*gsop.aminoCount*sizeof(GNativeParameters));
	// Building map
	int totalNative = 0;
	int i, j, k;
	for(k = 0; k < sop.nativeCount; k++){
		i = sop.natives[k].i;
		j = sop.natives[k].j;
		native.h_native[native.h_nativeCount[i]*gsop.aminoCount + i] = j;
		native.h_native[native.h_nativeCount[j]*gsop.aminoCount + j] = i;
		native.h_nativeParameters[native.h_nativeCount[i]*gsop.aminoCount + i].r02 = sop.natives[k].r0*sop.natives[k].r0;
		native.h_nativeParameters[native.h_nativeCount[j]*gsop.aminoCount + j].r02 = sop.natives[k].r0*sop.natives[k].r0;
		native.h_nativeParameters[native.h_nativeCount[i]*gsop.aminoCount + i].minus12ehOverR02 =
														-12.0*sop.natives[k].eh/(sop.natives[k].r0*sop.natives[k].r0);
		native.h_nativeParameters[native.h_nativeCount[j]*gsop.aminoCount + j].minus12ehOverR02 =
														-12.0*sop.natives[k].eh/(sop.natives[k].r0*sop.natives[k].r0);
		native.h_nativeCount[i] ++;
		native.h_nativeCount[j] ++;
		if(native.h_nativeCount[i] > native.max_native || native.h_nativeCount[j] > native.max_native){
			printf("ERROR: Maximum number of native contacts exceeded the limit of %d.\n", native.max_native);
			exit(-1);
		}
		totalNative++;
	}

	for(j = 1; j < gsop.Ntr; j++){
		for(i = 0; i < sop.aminoCount; i++){
			for(k = 0; k < native.max_native; k++){
				native.h_native[j*sop.aminoCount + i + k*gsop.aminoCount] =
							native.h_native[i + k*gsop.aminoCount] + j*sop.aminoCount;
				native.h_nativeParameters[j*sop.aminoCount + i + k*gsop.aminoCount] =
							native.h_nativeParameters[i + k*gsop.aminoCount];
			}
			native.h_nativeCount[j*sop.aminoCount + i] = native.h_nativeCount[i];
		}
	}


	#ifdef DEBUG1
	printf("Native contacts (number of contacts, #s of beads, zeros):\n");
	for(int i = 0; i < gsop.aminoCount; i++){
		printf("%d (%d): ", i, native.h_nativeCount[i]);
		for(int j = 0; j < native.h_nativeCount[i]; j++){
			printf("%d(%3.1f,%3.1f)  ", native.h_native[j*gsop.aminoCount + i],
					native.h_nativeParameters[j*gsop.aminoCount + i].x,
					native.h_nativeParameters[j*gsop.aminoCount + i].y);
		}
		printf("\n");
	}
	#endif
	// Copying to device
	cudaMemcpy(native.d_native, native.h_native, native.max_native*gsop.aminoCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(native.d_nativeCount, native.h_nativeCount, gsop.aminoCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(native.d_nativeParameters, native.h_nativeParameters, native.max_native*gsop.aminoCount*sizeof(float2), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_native, &native, sizeof(Native), 0, cudaMemcpyHostToDevice);

	printf("Total number of native contacts: %d \n", totalNative);
}

inline void computeNative(){
	native_kernel<<<native.blockNum, native.blockSize>>>();
	checkCUDAError();
}

inline void computeNativeEnergy(){
	nativeEnergy_kernel<<<native.blockNum, native.blockSize>>>();
	checkCUDAError();
}
