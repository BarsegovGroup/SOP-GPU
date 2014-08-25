/*
 * langevin.cuh
 *
 *  Created on: Apr 27, 2010
 *      Author: zhmurov
 */

#ifndef LANGEVIN_CUH_
#define LANGEVIN_CUH_
#include "../gsop.cuh"

#define LANGEVIN_TIMESTEP_STRING	"timestep"
#define LANGEVIN_ZETA_STRING		"zeta"

#define DEFAULT_LANGEVIN_ZETA		50.0

#define DEFAULT_TEMPERATURE 0.6

#define TC_TEMPERATURE_STRING	"temperature"
#define TC_HEATING_STRING		"heating"
#define TC_INITIAL_T_STRING 	"initialT"
#define TC_DELTA_T_STRING		"deltaT"

struct Langevin {
	float h;
	float zeta;
	float var;
	float hOverZeta;
	float tempNorm;
	float temp;
	float initialT;
	float deltaT;
};

extern Langevin langevin;

void createLangevinIntegrator();
void initLangevinIntegrator();
void integrateLangevin();
void deleteLangevinIntegrator();

void createTemperatureControl();
void initTemperatureControl();
void updateTemperature();
void deleteTemperatureUpdater();

#endif /* LANGEVIN_CUH_ */
