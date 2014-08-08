/*
 * parameters_initializer.c
 *
 *  Created on: Jan 16, 2009
 *      Author: zhmurov
 */

#include "def_param.h"
#include "IO/configreader.h"
#include "gsop.h"

#define buf_size 2048

int mode;

int BLOCK_SIZE = 256;

long int numsteps;

int nav;

char top_filename[100];
char coord_filename[100];
char ref_filename[100];
char final_filename[100];

void initParameters(const char* configFile){

	parseParametersFile(configFile);

    char protein_name[100];
	getParameter(protein_name, "name", "unnamed", 1);
	printf("Initializing simulations for '%s'\n", protein_name);


	char modeString[100];
	getParameter(modeString, "mode", "", 1);
	if(strcmp(modeString, "capsid") == 0){
		mode = MODE_CAPSID;
	} else {
		mode = 0;
	}

	BLOCK_SIZE = getIntegerParameter("block_size", BLOCK_SIZE, 1);

	int run = getIntegerParameter("run", -1, 1);
	if(run == -1){
		gsop.Ntr = getIntegerParameter("runnum", 0, 0);
		gsop.firstrun = getIntegerParameter("firstrun", 0 ,0);
	} else {
		gsop.Ntr = 1;
		gsop.firstrun = run;
	}

	numsteps = getLongIntegerParameter("numsteps", 0, 0);

	nav = getIntegerParameter("nav", 1000, 1);
	
    gsop.deviceId = getIntegerParameter("device", 0, 1);

	getMaskedParameter(top_filename, "topology", "", 0);
	getMaskedParameter(coord_filename, "coordinates", "", 0);

	getMaskedParameter(ref_filename, "reffilename", "<name>.ref.pdb", 1);
	//printf("Coordinates will be saved into '%s'\n", dcd_filename);
	getMaskedParameter(final_filename, "finalcoord", "<name>_<author><run>_<stage>_final.pdb", 1);
}

