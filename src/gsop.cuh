#ifndef GSOP_CUH
#define GSOP_CUH
/*
 * gsop.cuh
 *
 *  Created on: Jun 4, 2009
 *      Author: zhmurov
 */

#include "gsop.h"

__device__ __constant__ GSOP c_gsop;

#ifndef NOTEXTURE
texture<float4, 1, cudaReadModeElementType> t_coord; // Coordinates
#endif

extern void __checkCUDAError(const char *file, int line);
extern void copyCoordDeviceToHost();
#define checkCUDAError() __checkCUDAError(__FILE__, __LINE__)

#endif
