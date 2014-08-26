#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "def_param.h"
#include "gsop.h"

#include "param_initializer.h"
#include "IO/configreader.h"
#include "Util/wrapper.h"

long long int initialTime;
long long int lastTime;

long int step;

SOP sop;

int main(int argc, char *argv[]){

	printf("===========================\n");
	printf("SOP-GPU version %s\n", VERSION);
	printf("===========================\n");

#ifdef DEBUG
	printf("Running in DEBUG mode.\n");
#endif
	initialTime = time(NULL);
	step = 0;

	int i;

	if(argc < 1){
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

    char ref_filename[100];
    getMaskedParameter(ref_filename, "reffilename", "<name>.ref.pdb", 1);
    printf("Saving reference PDB: '%s'.\n", ref_filename);
	savePDB(ref_filename, sop); // Save reference PDB at the begining of the trajectory

	copyCoordinates(); // Copy all coordinates to a gpu
	initFF(); // Initialize all potentials and updaters

	lastTime = time(NULL); // Reset timer

	printf("Starting %d runs (%ld steps).\n", gsop.Ntr, numsteps);
	runGPU(); // Strart simulations
	//savePDB(final_filename);
	// All done. Release memory/close files.
	//cleanup();
}

