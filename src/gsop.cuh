#pragma once

/*
 * gsop.cuh
 *
 *  Created on: Jun 4, 2009
 *      Author: zhmurov
 */

#include "gsop.h"

__device__ __constant__ GSOP c_gsop;

extern cudaDeviceProp deviceProp;
extern int BLOCK_SIZE;

#ifndef NOTEXTURE
texture<float4, 1, cudaReadModeElementType> t_coord; // Coordinates
#endif

void __checkCUDAError(const char *file, int line);
void copyCoordDeviceToHost();
#define checkCUDAError() __checkCUDAError(__FILE__, __LINE__)

