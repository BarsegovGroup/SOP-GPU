/*
 * langevin.cuh
 *
 *  Created on: Apr 27, 2010
 *      Author: zhmurov
 */

#ifndef LANGEVIN_CUH_
#define LANGEVIN_CUH_
#include "../gsop.h"

#define LANGEVIN_TIMESTEP_STRING	"timestep"
#define LANGEVIN_ZETA_STRING		"zeta"

#define DEFAULT_LANGEVIN_ZETA		50.0

#define DEFAULT_TEMPERATURE 0.6

#define TC_TEMPERATURE_STRING	"temperature"
#define TC_HEATING_STRING		"heating"
#define TC_INITIAL_T_STRING 	"initialT"
#define TC_DELTA_T_STRING		"deltaT"

class LangevinIntegrator : public SOPIntegrator {
public:
    LangevinIntegrator();
    virtual ~LangevinIntegrator();
    virtual void integrate();

	float zeta;
	float var;
	float hOverZeta;
	float tempNorm;
	float temp;
};

class TemperatureUpdater : public SOPUpdater {
public:
    TemperatureUpdater(LangevinIntegrator *langevin);
    virtual ~TemperatureUpdater();
    virtual void update();
private:
    LangevinIntegrator *langevin;
	float initialT;
	float deltaT;
};

void createLangevinIntegrator();

#endif /* LANGEVIN_CUH_ */
