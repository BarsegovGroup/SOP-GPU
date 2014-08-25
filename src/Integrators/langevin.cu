/*
 * langevin.cu
 *
 *  Created on: Apr 27, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "ht.cu"
#include "langevin.h"

#include "langevin_kernel.cu"

void createLangevinIntegrator(){
	integrator = new LangevinIntegrator();
}

LangevinIntegrator::LangevinIntegrator(){
	this->name = "Langevin";
	this->h = getFloatParameter(LANGEVIN_TIMESTEP_STRING, 0, 0);
	this->zeta = getFloatParameter(LANGEVIN_ZETA_STRING, DEFAULT_LANGEVIN_ZETA, 1);
	this->hOverZeta = this->h/this->zeta;
	this->tempNorm = this->zeta/(6.0*this->h);
    int seed = getIntegerParameter("seed", time(NULL), 1) + gsop.firstrun + getIntegerParameter("run", -1, 1); // TODO: do we need `run` here?
	initRand(seed, gsop.aminoCount);

	int heating = getYesNoParameter(TC_HEATING_STRING, 0, 1);
	if(heating == 1 || gsop.heatingOn == 1){
		gsop.heatingOn = 1;
		printf("Initializing heating protocol...\n");
		updaters[updatersCount] = new TemperatureUpdater(this);
		updatersCount++;
		this->temp = getFloatParameter(TC_INITIAL_T_STRING, 0, 0);
	} else {
		printf("Simulations will be held at constant temperature.\n");
		this->temp = getFloatParameter(TC_TEMPERATURE_STRING, DEFAULT_TEMPERATURE, 1);
	}
	this->var = sqrtf(2.0f* this->temp * this->zeta/this->h);

    hc_langevin.var = this->var;
    hc_langevin.hOverZeta = this->hOverZeta;
    hc_langevin.tempNorm = this->tempNorm;
	cudaMemcpyToSymbol(c_langevin, &hc_langevin, sizeof(LangevinConstant), 0, cudaMemcpyHostToDevice);

	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(integrateLangevin_kernel, cudaFuncCachePreferL1);
	}
}

void LangevinIntegrator::integrate(){
	integrateLangevin_kernel<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
}

LangevinIntegrator::~LangevinIntegrator(){

}

TemperatureUpdater::TemperatureUpdater(LangevinIntegrator* langevin){
    this->name = "TemperatureUpdater";
	this->frequency = getIntegerParameter("tempFreq", 0, 0);
    this->langevin = langevin;
	this->initialT = getFloatParameter(TC_INITIAL_T_STRING, 0, 0);
	this->deltaT = getFloatParameter(TC_DELTA_T_STRING, 0, 0);
}

void TemperatureUpdater::update(){
	if(step % this->frequency == 0){
		langevin->temp = this->initialT + (step/this->frequency)*this->deltaT;
		langevin->var = sqrtf(2.0f* langevin->temp * langevin->zeta/langevin->h);

        hc_langevin.var = langevin->var;
		cudaMemcpyToSymbol(c_langevin, &hc_langevin, sizeof(LangevinConstant), 0, cudaMemcpyHostToDevice);
		printf("Heating: %ld\t%f\n", step, langevin->temp);
	}
}

TemperatureUpdater::~TemperatureUpdater(){

}

