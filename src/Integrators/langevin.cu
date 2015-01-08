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
	this->h = parameters::timestep.get();
	this->zeta = parameters::zeta.get();
	this->hOverZeta = this->h/this->zeta;
	this->tempNorm = this->zeta/(6.0*this->h);
    int seed = parameters::seed.get() + gsop.firstrun + getIntegerParameter("run", -1, 1); // TODO: do we need `run` here?
	initRand(seed, gsop.aminoCount);

	int heating = parameters::heating.get();
	if(heating == 1 || gsop.heatingOn == 1){ // TODO: Remove these double-checks from everywhere!
		gsop.heatingOn = 1;
		printf("Initializing heating protocol...\n");
		updaters[updatersCount] = new TemperatureUpdater(this);
		updatersCount++;
		this->setTemperature(parameters::initialT.get());
	} else {
		printf("Simulations will be held at constant temperature.\n");
		this->setTemperature(parameters::temperature.get());
	}

	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(integrateLangevin_kernel, cudaFuncCachePreferL1);
	}
    this->updateParametersOnGPU();
}

void LangevinIntegrator::integrate(){
	integrateLangevin_kernel<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
}

void LangevinIntegrator::setTemperature(float temp){
    this->temp = temp;
    this->var = sqrtf(2.0f* this->temp * this->zeta/this->h);
    this->updateParametersOnGPU();
}

void LangevinIntegrator::updateParametersOnGPU(){
    hc_langevin.var = this->var;
    hc_langevin.hOverZeta = this->hOverZeta;
    hc_langevin.tempNorm = this->tempNorm;
	cudaMemcpyToSymbol(c_langevin, &hc_langevin, sizeof(LangevinConstant), 0, cudaMemcpyHostToDevice);
    checkCUDAError();
}

LangevinIntegrator::~LangevinIntegrator(){

}

TemperatureUpdater::TemperatureUpdater(LangevinIntegrator* langevin){
    this->name = "TemperatureUpdater";
	this->frequency = parameters::tempFreq.get();
    this->langevin = langevin;
	this->initialT = parameters::initialT.get();
	this->deltaT = parameters::deltaT.get();
}

void TemperatureUpdater::update(){
	if(gsop.step % this->frequency == 0){
        float temp = this->initialT + (gsop.step/this->frequency)*this->deltaT;
		langevin->setTemperature( temp );
		printf("Heating: %ld\t%f\n", gsop.step, temp);
	}
}

TemperatureUpdater::~TemperatureUpdater(){

}

