#include "../gsop.h"
#include "output_manager.h"

__global__ void reset_temperature(){
	int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_gsop.aminoCount){
		c_gsop.d_energies[d_i].w = 0.0f;
    }
}

void OutputManager::resetTemperatureCounter() {
    const int blockSize = gsop.blockSize;
    const int blockNum = gsop.aminoCount/blockSize + 1;
    reset_temperature<<<blockNum, blockSize>>>();
}

