#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "def_param.h"
#include "gsop.h"
#include "Util/wrapper.h"
#define buf_size 2048

long long int initialTime;
long long int lastTime;

long int step;

PDB pdbdata;
SOP sop;

extern void initGPU();
extern void initFF();
extern void copyCoordinates();
extern void copyCoordinatesTrajectory(int traj);
extern void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace);
extern void runGPU();
extern void cleanup();

extern void initParameters(const char* configFile);

void getOutputFilnames();

int main(int argc, char *argv[]){

	printf("==========================\n");
	printf("gSOP version %s\n", version);
	printf("==========================\n");

#ifdef DEBUG
	printf("Running in DEBUG mode.\n");
#endif
	initialTime = time(NULL);
	step = 0;

	int i;

	if(argc < 1){
		DIE("ERROR: Configuration file should be specified.");
	}

	int restart = 0;

	for(i = 1; i < argc; i++){
		if(strcmp(argv[i], "--restart") == 0){
			//printf("Reading restart parameters...\n");
			restart = 1;
			DIE("Restart is not currently supported. Sorry.");
		}
	}

	// Read main parameters of the simulation (number of steps, etc.)
	initParameters(argv[1]);

   	if(top_filename == NULL){
		DIE("Please, use sop-top utility to create topology.");
	} else {
		sop.load(top_filename);
	}

	initGPU(); // Set device, allocate memory for the common variables (coordinates/forces)
	if(restart == 0){ //Restart feature is not working at a time (condition always satisfied)
		for(i = 0; i < gsop.Ntr; i++){
			char trajnum[10];
			char trajCoordFilename[100];
			sprintf(trajnum, "%d", i+gsop.firstrun);
			replaceString(trajCoordFilename, coord_filename, trajnum, "<run>");
			readCoord(trajCoordFilename, sop); // Read coordinates from a file for all (or single) trajectories
			copyCoordinatesTrajectory(i); // Copy coordinates to a specific location in a coordinates array
		}
	} else {
		FILE* src = fopen(restartpdb_filename, "r");
		char bckFilename[100];
		sprintf(bckFilename, "%s.%ld.bck", restartpdb_filename, step);
		FILE* bck = fopen(bckFilename, "w");
		char buffer[buf_size];
		while(fgets(buffer, buf_size, src) != NULL){
			fputs(buffer, bck);
		}
		fclose(bck);
		fclose(src);
		readCoord(restartpdb_filename, sop);
	}
	copyCoordinates(); // Copy all coordinates to a gpu
	initFF(); // Initialize all potentials and updaters

	printf("Saving reference PDB: '%s'.\n", ref_filename);
	savePDB(ref_filename, sop); // Save reference PDB at the begining of the trajectory
	lastTime = time(NULL); // Reset timer

	/*if(restart == 1){
		printf("Reading restart coordinates from %s.\n", restartpdb_filename);
		readCoord(restartpdb_filename);
		copyCoordinates();
	}*/

	printf("Starting %d runs (%ld steps).\n", gsop.Ntr, numsteps);
	runGPU(); // Strart simulations
	//savePDB(final_filename);
	// All done. Release memory/close files.
	//cleanup();
}

