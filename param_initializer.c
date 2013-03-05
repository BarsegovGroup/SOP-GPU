/*
 * parameters_initializer.c
 *
 *  Created on: Jan 16, 2009
 *      Author: zhmurov
 */

#include "def_param.h"
#include "config_reader.h"

#define buf_size 2048

int device;
int mode;

int BLOCK_SIZE = 256;

long int numsteps;
int DCDfreq;
int restartfreq;

float R_sigma;
float eh;
int nav;

int pairs_freq;
int possiblepairs_freq;

int restart;

int seed;
int run;
int Ntr;
int firstrun;

char protein_name[100];
char stageString[100];

int readInitialCoord = 0;
char top_filename[100];
char coord_filename[100];
char ref_filename[100];
char final_filename[100];

void initParameters(char* configFile);
void makeRestartFilename(char* filename);

void initParameters(char* configFile, int createTopology){

	parseFile(configFile);

	getParameter(protein_name, "name", "unnamed", 1);
	printf("Initializing simulations for '%s'\n", protein_name);

	if(createTopology){

	} else {

		device = getIntegerParameter("device", 0, 1);

		char modeString[100];
		getParameter(modeString, "mode", "", 1);
		if(strcmp(modeString, "capsid") == 0){
			mode = MODE_CAPSID;
		} else {
			mode = 0;
		}

		BLOCK_SIZE = getIntegerParameter("block_size", BLOCK_SIZE, 1);

		run = getIntegerParameter("run", -1, 1);
		if(run == -1){
			Ntr = getIntegerParameter("runnum", 0, 0);
			firstrun = getIntegerParameter("firstrun", 0 ,0);
		} else {
			Ntr = 1;
			firstrun = run;
		}

		getMaskedParameter(stageString, "stage", "", 0);
		heatingOn = 0;
		pullingOn = 0;
		indentationOn = 0;
		minimizationOn = 0;
		if(strcmp(stageString, "heat") == 0){
			printf("Preparing heating simulations.\n");
			heatingOn = 1;
		} else if(strcmp(stageString, "minim") == 0){
			printf("Preparing minimization simulations.\n");
			minimizationOn = 1;
		} else if(strcmp(stageString, "pull") == 0){
			printf("Preparing pulling simulations.\n");
			pullingOn = 1;
		} else if(strcmp(stageString, "indent") == 0){
			printf("Preparing indentation simulations.\n");
			indentationOn = 1;
		}

		seed = getIntegerParameter("seed", time(NULL), 1) + run + firstrun;
		srand(seed);

		numsteps = getLongIntegerParameter("numsteps", 0, 0);
		DCDfreq = getIntegerParameter("dcdfreq", 10000, 1);
		restartfreq = getIntegerParameter("restartfreq", 100000, 1);

		R_sigma = getFloatParameter("R_sigma", 16.0f, 1);

		nav = getIntegerParameter("nav", 1000, 1);

		pairs_freq = getIntegerParameter("pairs_freq", 1000, 1);
		possiblepairs_freq = getIntegerParameter("possiblepairs_freq", 100000, 1);

		getMaskedParameter(top_filename, "topology", "", 0);
		getMaskedParameter(coord_filename, "coordinates", "", 0);

		getMaskedParameter(ref_filename, "reffilename", "<name>.ref.pdb", 1);
		//printf("Coordinates will be saved into '%s'\n", dcd_filename);
		getMaskedParameter(final_filename, "finalcoord", "<name>_<author><run>_<stage>_final.pdb", 1);


	}

}

/*void initRestartParameters(char* configFile){
	printf("Reading restart parameters from '%s'.\n", configFile);
	parseFile(configFile);
	step = getLongIntegerParameter("timeStep", 0, 0);
	pullVectorX = getFloatParameter("pullVector.x", 0.0f, 0);
	pullVectorY = getFloatParameter("pullVector.y", 0.0f, 0);
	pullVectorZ = getFloatParameter("pullVector.z", 0.0f, 0);
	xt = getFloatParameter("extension", 0.0f, 0);
	chipCoordX = getFloatParameter("chipCoord.x", 0.0f, 0);
	chipCoordY = getFloatParameter("chipCoord.y", 0.0f, 0);
	chipCoordZ = getFloatParameter("chipCoord.z", 0.0f, 0);
	makeRestartFilename(ref_filename);
	makeRestartFilename(dat_filename);
	makeRestartFilename(dcd_filename);
	printf("Coordinates will be saved into '%s'\n", dcd_filename);
	makeRestartFilename(final_filename);
	seed += step;
	srand(seed);
	printf("Resuming simulations from step %ld.\n", step);
	FILE* src = fopen(configFile, "r");
	char bckFilename[100];
	sprintf(bckFilename, "%s.%ld.bck", configFile, step);
	FILE* bck = fopen(bckFilename, "w");
	char buffer[buf_size];
	while(fgets(buffer, buf_size, src) != NULL){
		fputs(buffer, bck);
	}
	fclose(bck);
	fclose(src);
}*/

void makeRestartFilename(char* filename){
	char tmp[100];
	strcpy(tmp, filename);
	char* pch = strtok(tmp, ".");
	sprintf(filename, "%s.restart%ld.%s", pch, step, strtok(NULL, ""));
	printf("%s\n", filename);
}
