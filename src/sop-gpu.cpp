#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "gsop.h"

#include "Util/wrapper.h"
#include "Util/parameters.h"

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

    std::string top_filename = parameters::topology.get();
	sop.load(top_filename.c_str());

	initGPU(); // Set device, allocate memory for the common variables (coordinates/forces)
    std::string coord_filename;
	for(i = 0; i < gsop.Ntr; i++){
        coord_filename = parameters::coordinates.replace("<run>", i+gsop.firstrun);
		readCoord(coord_filename.c_str(), sop); // Read coordinates from a file for all (or single) trajectories
		copyCoordinatesTrajectory(i); // Copy coordinates to a specific location in a coordinates array
	}

	copyCoordinates(); // Copy all coordinates to a gpu
	initFF(); // Initialize all potentials and updaters

    std::string ref_filename = parameters::reffilename.get();
    printf("Saving reference PDB: '%s'.\n", ref_filename.c_str());
	savePDB(ref_filename.c_str(), sop); // Save reference PDB at the begining of the trajectory

	printf("System initialization took %.f second(s).\n", difftime(time(NULL), initialTime));

	printf("Starting %d runs (%ld steps).\n", gsop.Ntr, gsop.numsteps);
	runGPU(); // Strart simulations
	//savePDB(final_filename);
	// All done. Release memory/close files.
	//cleanup();
}

void initParameters(const char* configFile){

    parameters::_initialize(configFile);

	printf("Initializing simulations for '%s'\n", parameters::name.get().c_str());


	int run = parameters::run.get();
	if(run == -1){
		gsop.Ntr = parameters::runnum.get();
		gsop.firstrun = parameters::firstrun.get();
	} else {
		gsop.Ntr = 1;
		gsop.firstrun = run;
	}

	gsop.numsteps = parameters::numsteps.get();

	gsop.nav = parameters::nav.get();

    gsop.deviceId = parameters::device.get();
}

