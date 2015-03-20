#pragma once

/*
 * gsop.cuh
 *
 *  Created on: Jun 4, 2009
 *      Author: zhmurov
 */

#include "gsop.h"
#include "Util/parameters.h"

PARAMETER(block_size, int, 256, "", "...")
PARAMETER(finalcoord, std::string, "<name>_<run>_final.pdb", "path", "...")

__device__ __constant__ GSOP c_gsop;

#ifndef NOTEXTURE
texture<float4, 1, cudaReadModeElementType> t_coord; // Coordinates
#endif

