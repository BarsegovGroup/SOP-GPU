#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "gsop.h"

#include "IO/configreader.h"
#include "Util/wrapper.h"

SOP sop;

void initParameters(const char* configFile);

int main(int argc, char *argv[]){

	time_t initialTime = time(NULL);

	printf("===========================\n");
	printf("SOP-GPU version %s\n", VERSION);
	printf("===========================\n");

#ifdef DEBUG
	printf("Running in DEBUG mode.\n");
#endif
	gsop.step = 0;

	int i;

	if(argc < 2){
		DIE("ERROR: Configuration file should be specified.");
	}

	// Read main parameters of the simulation (number of steps, etc.)
	initParameters(argv[1]);

    char top_filename[100];
	getMaskedParameter(top_filename, "topology", "", 0);
	sop.load(top_filename);

	initGPU(); // Set device, allocate memory for the common variables (coordinates/forces)
    char coord_filename[100];
	for(i = 0; i < gsop.Ntr; i++){
        getMaskedParameterWithReplacementT(coord_filename, "coordinates", i+gsop.firstrun, "<run>");
		readCoord(coord_filename, sop); // Read coordinates from a file for all (or single) trajectories
		copyCoordinatesTrajectory(i); // Copy coordinates to a specific location in a coordinates array
	}

	copyCoordinates(); // Copy all coordinates to a gpu
	initFF(); // Initialize all potentials and updaters

    char ref_filename[100];
    getMaskedParameter(ref_filename, "reffilename", "<name>.ref.pdb", 1);
    printf("Saving reference PDB: '%s'.\n", ref_filename);
	savePDB(ref_filename, sop); // Save reference PDB at the begining of the trajectory

	printf("System initialization took %.f second(s).\n", difftime(time(NULL), initialTime));

	printf("Starting %d runs (%ld steps).\n", gsop.Ntr, gsop.numsteps);
	runGPU(); // Strart simulations
	//savePDB(final_filename);
	// All done. Release memory/close files.
	//cleanup();
}

void initParameters(const char* configFile){

	parseParametersFile(configFile);

    char protein_name[100];
	getParameter(protein_name, "name", "unnamed", 1);
	printf("Initializing simulations for '%s'\n", protein_name);


	int run = getIntegerParameter("run", -1, 1);
	if(run == -1){
		gsop.Ntr = getIntegerParameter("runnum", 0, 0);
		gsop.firstrun = getIntegerParameter("firstrun", 0 ,0);
	} else {
		gsop.Ntr = 1;
		gsop.firstrun = run;
	}

	gsop.numsteps = getLongIntegerParameter("numsteps", 0, 0);

	gsop.nav = getIntegerParameter("nav", 1000, 1);

    gsop.deviceId = getIntegerParameter("device", 0, 1);
}

