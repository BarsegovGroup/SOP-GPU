/*
 * output_manager.cuh
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#pragma once

#include "../gsop.h"

#include "../Potentials/native.h"
#include "../Potentials/covalent.h"
#include "../Integrators/bdhitea.h"

#include <vector>
#include <string>

#define OUTPUT_FREQUENCY_STRING			"outputtiming"
#define OUTPUT_FILENAME_STRING			"outputname"

#define DEFAULT_OUTPUT_FREQUENCY		10000
#define DEFAULT_OUTPUT_FILENAME			"energy.<name>_<run>_<stage>.dat"

struct OutputData{
	long int step;
	float tempav;
	float epot_tot;
	float epot_native;
	float epot_longrange;
	float epot_LJ;
	float epot_fene;
	int nat_num;
	float rg;
	float tea_eps;
};

class OutputManager : public SOPUpdater{
public:
    OutputManager();
    virtual ~OutputManager();
    virtual void update();
private:
    void computeEnergies(int traj);
    void computeNativeNumber(int traj);
    void computeRg(int traj);
    void computeTEAeps(int traj);
    void printDataToScreen() const;
    void printDataToFile(FILE* dat_file) const;
    void resetTemperatureCounter();

    CovalentPotential *covalent;
    NativePotential *native;
    TeaIntegrator *tea;

    OutputData outputData;
    std::vector<std::string> dat_filenames;
};

void createOutputManager();

