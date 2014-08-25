/*
 * parameters_initializer.c
 *
 *  Created on: Jan 16, 2009
 *      Author: zhmurov
 */

#include "def_param.h"
#include "IO/configreader.h"
#include "gsop.h"

int BLOCK_SIZE = 256;

long int numsteps;

int nav;

void initParameters(const char* configFile){

	parseParametersFile(configFile);

    char protein_name[100];
	getParameter(protein_name, "name", "unnamed", 1);
	printf("Initializing simulations for '%s'\n", protein_name);

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
}

