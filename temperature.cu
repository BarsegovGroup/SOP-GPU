/*
 * temperature.cu
 *
 *  Created on: Jun 26, 2009
 *      Author: zhmurov
 */

#include "gsop.cuh"

void initTemperature(){
	if(stage == heat_stage){
		temp = initialT;
	}
	gsop.var = sqrtf(2.0f * temp * zeta/dh);
}

void increaseTemperature(){
	if(step % tempFreq == 0){
		temp = initialT + (step/tempFreq)*deltaT;
		gsop.var = sqrtf(2.0f* temp * zeta/dh);
		cudaMemcpyToSymbol(c_gsop, &gsop, sizeof(GSOP), 0, cudaMemcpyHostToDevice);
	}
}
