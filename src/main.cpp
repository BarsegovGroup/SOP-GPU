#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/prctl.h>

#include "def_param.h"
#include "Util/wrapper.h"
#define buf_size 2048

long long int initialTime;
long long int lastTime;
long long int switchTime;

long int step;

float forcex, forcey, forcez, xt, x0;
float epot_fene, epot_LJ, tempav, epot_LJ_att, epot_LJ_nei, epot_LJ_rep;
float x_R, R, rg;

int pullingOn = 0;
int heatingOn = 0;
int indentationOn = 0;
int minimizationOn = 0;

PDB pdbdata;
SOP sop;


extern void initGPU();
extern void initFF();
extern void copyCoordinates();
extern void copyCoordinatesTrajectory(int traj);
extern void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace);
extern void runGPU();
extern void cleanup();

extern void initParameters(const char* configFile, int createTopology);

void getOutputFilnames();

extern FILE* dcd_open_read(char *dcd_filename);
extern void dcd_read_header(FILE* dcd_file, char *dcd_filename, int* N, int* NFILE, int* NPRIV, int* NSAVC, float* DELTA);
extern int dcd_read_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z);
extern void dcd_close(FILE* dcd_file);

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
	char   *str;

	if(argc < 1){
		DIE("ERROR: Configuration file should be specified.");
	}

	restart = 0;
	int createTopology = 0;

	for(i = 1; i < argc; i++){
		if(strcmp(argv[i], "--make-top") == 0){
			createTopology = 1;
		}
		if(strcmp(argv[i], "--restart") == 0){
			//printf("Reading restart parameters...\n");
			restart = 1;
			DIE("Restart is not currently supported. Sorry.");
		}
	}

	str = argv[1];
	// Read main parameters of the simulation (number of steps, etc.)
	initParameters(str, createTopology);

    if(createTopology){
		DIE("Please, use sop-top utility to create topology.");
	}

	if(top_filename == NULL){
		DIE("Please, use sop-top utility to create topology.");
	} else {
		sop.load(top_filename);
	}
	char processNameFull[100];
	char processName[100];
	if(run != -1){
		sprintf(processNameFull, "gsop_%d_%d", run, device);
	} else {
		sprintf(processNameFull, "gsop_%d_%d_%d", firstrun, firstrun+Ntr, device);
	}
	strncpy(processName, processNameFull, 15);
	processName[15] = '\0';

	prctl(PR_SET_NAME, processName);

	initGPU(); // Set device, allocate memory for the common variables (coordinates/forces)
	if(restart == 0){ //Restart feature is not working at a time (condition always satisfied)
		for(i = 0; i < Ntr; i++){
			char trajnum[10];
			char trajCoordFilename[100];
			sprintf(trajnum, "%d", i+firstrun);
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

	if(restart == 0){
		xt = 0.0;
	}
	printf("Starting %d runs (%ld steps).\n", Ntr, numsteps);
	runGPU(); // Strart simulations
	//savePDB(final_filename);
	// All done. Release memory/close files.
	//cleanup();

}

