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
#include "../Util/parameters.h"

#include <vector>
#include <string>


PARAMETER(outputtiming, int, 10000, "steps", "...")
PARAMETER(outputname, std::string, "energy.<name>_<run>_<stage>.dat", "path", "...")
PARAMETER(outputcolwidth, int, 16, "symbols", "...")
PARAMETER(printruns, int, 10, "", "...")
PARAMETER(computeRg, bool, true, "true/false", "...")

struct OutputData{
	float tea_eps;
};

class OutputManager : public SOPUpdater{
public:
    OutputManager();
    virtual ~OutputManager();
    virtual void update();
private:
    void computeNativeCounts();
    void resetTemperatures();
    void computeTemperatures();
    void computeRgs();
    float getRg(int traj);
    void computeTEAeps(int traj);
    void printTimeEstimates();
    void printDataTitleToScreen() const;
    void printDataToScreen(int traj) const;
    void printDataLdotsToScreen() const;
    void printDataToFile(FILE* dat_file, int traj) const;

    long long int initialTime;
    long long int lastTime;

    OutputData outputData;
    std::vector<std::string> dat_filenames;

    // Gyration radii
    float* d_rgs; //For all particles (on GPU)
    float* h_rgs; //For all particles (on CPU)
    float* rgs; //For all trajectories

    float* temperatures; // Temperatures for all trajectories
    int* nativeCounts; // Native contacts for all trajectories

    // Cutoff values for native contacts counting
    float R_limit;
    float R_limit_bond;

    int printRuns; // How many runs data is to be printed into standard output
    int outputWidth; // Width of the columns in standard output
    bool computeRgFlag; // Are gyration radii needed
};

void createOutputManager();

