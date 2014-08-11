/*
 * langevin.cu
 *
 *  Created on: Apr 27, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "ht.cu"
#include "langevin.cuh"
#include "langevin_kernel.cu"

void createLangevinIntegrator(){
	sprintf(langevinIntegrator.name, "Langevin");
	langevinIntegrator.integrate = &integrateLangevin;
	langevinIntegrator.destroy = &deleteLangevinIntegrator;
	integrator = &langevinIntegrator;
	initLangevinIntegrator();
}

void initLangevinIntegrator(){
	langevin.h = getFloatParameter(LANGEVIN_TIMESTEP_STRING, 0, 0);
	langevinIntegrator.h = langevin.h;
	langevin.zeta = getFloatParameter(LANGEVIN_ZETA_STRING, DEFAULT_LANGEVIN_ZETA, 1);
	langevin.hOverZeta = langevin.h/langevin.zeta;
	langevin.tempNorm = langevin.zeta/(6.0*langevin.h);
    int seed = getIntegerParameter("seed", time(NULL), 1) + gsop.firstrun + getIntegerParameter("run", -1, 1); // TODO: do we need `run` here?
	initRand(seed, gsop.aminoCount);
	createTemperatureControl();
	cudaMemcpyToSymbol(c_langevin, &langevin, sizeof(Langevin), 0, cudaMemcpyHostToDevice);
	if(gsop.deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(integrateLangevin_kernel, cudaFuncCachePreferL1);
	}
}

void integrateLangevin(){
	integrateLangevin_kernel<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
}

void deleteLangevinIntegrator(){

}

void createTemperatureControl(){
	int heating = getYesNoParameter(TC_HEATING_STRING, 0, 1);
	if(heating == 1 || gsop.heatingOn == 1){
		gsop.heatingOn = 1;
		printf("Initializing heating protocol...\n");
		sprintf(temperatureUpdater.name, "TemperatureControl");
		temperatureUpdater.update = &updateTemperature;
		temperatureUpdater.destroy = &deleteTemperatureUpdater;
		temperatureUpdater.frequency = getIntegerParameter("tempFreq", 0, 0);
		updaters[updatersCount] = &temperatureUpdater;
		updatersCount++;
	} else {
		printf("Simulations will be held at constant temperature.\n");
	}
	initTemperatureControl();
}

void initTemperatureControl(){
	if(gsop.heatingOn){
		langevin.initialT = getFloatParameter(TC_INITIAL_T_STRING, 0, 0);
		langevin.temp = langevin.initialT;
		langevin.deltaT = getFloatParameter(TC_DELTA_T_STRING, 0, 0);
	} else {
		langevin.temp = getFloatParameter(TC_TEMPERATURE_STRING, DEFAULT_TEMPERATURE, 1);
	}
	langevin.var = sqrtf(2.0f* langevin.temp * langevin.zeta/langevin.h);
}

void updateTemperature(){
	if(step % temperatureUpdater.frequency == 0){
		langevin.temp = langevin.initialT + (step/temperatureUpdater.frequency)*langevin.deltaT;
		langevin.var = sqrtf(2.0f* langevin.temp * langevin.zeta/langevin.h);
		cudaMemcpyToSymbol(c_langevin, &langevin, sizeof(Langevin), 0, cudaMemcpyHostToDevice);
		printf("Heating: %ld\t%f\n", step, langevin.temp);
	}
}

void deleteTemperatureUpdater(){

}

