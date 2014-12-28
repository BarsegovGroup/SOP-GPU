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
#define OUTPUT_DATA_WIDTH_STRING		"outputcolwidth"
#define OUTPUT_PRINT_RUNS_STRING		"printruns"
#define OUTPUT_COMPUTE_RG_STRING		"coputeRg"

#define DEFAULT_OUTPUT_FREQUENCY		10000
#define DEFAULT_OUTPUT_FILENAME			"energy.<name>_<run>_<stage>.dat"
#define DEFAULT_OUTPUT_DATA_WIDTH		16
#define DEFAULT_OUTPUT_PRINT_RUNS		10
#define DEFAULT_OUTPUT_COMPUTE_RG		1

const int MODE_CAPSID = 1;

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

