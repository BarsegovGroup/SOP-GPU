/*
 * dcd_manager.cu
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */
#include "dcd_manager.cuh"

char dcd_filename[100];
char restartpdb_filename[100];
char restartconf_filename[100];
int restartfreq;

char** dcd_filenames;
DCD dcd;
float* X;
float* Y;
float* Z;

void createDCDOutputManager(){
	printf("Initializing dcd output manager...\n");
	sprintf(dcdOutputManager.name, "DCD output");
	dcdOutputManager.update = &saveCoord;
	dcdOutputManager.frequency = getIntegerParameter(DCD_FREQUENCY_STRING, 10000, 1);;
	updaters[updatersCount] = &dcdOutputManager;
	updatersCount++;
	restartfreq = getIntegerParameter("restartfreq", 100000, 1);
	initDCD();
	printf("Done initializing dcd output manager...\n");
}


/*
 * Initialization of coordinates output
 */
void initDCD(){

	getMaskedParameter(dcd_filename, "DCDfile", "<name>_<author><run>_<stage>.dcd", 1);
	//printf("Coordinates will be saved into '%s'\n", dcd_filename);
	char restart_filename[100];
	getMaskedParameter(restart_filename, "restartname", "<name>_<author><run>_restart", 1);
	sprintf(restartpdb_filename, "%s.pdb", restart_filename);
	sprintf(restartconf_filename, "%s.conf", restart_filename);
	getMaskedParameter(final_filename, "finalcoord", "<name>_<author><run>_<stage>_final.pdb", 1);

	int traj;
	dcd_filenames = (char**)calloc(gsop.Ntr, sizeof(char*));
	for(traj = 0; traj < gsop.Ntr; traj++){
		//printf("%d.1\n", traj);
		dcd_filenames[traj] = (char*)calloc(100, sizeof(char));
		char trajnum[10];
		sprintf(trajnum, "%d", traj+gsop.firstrun);
		replaceString(dcd_filenames[traj], dcd_filename, trajnum, "<run>");
		//printf("%d.2\n", traj);
		int particleCount = sop.aminoCount;
		if(sop.additionalAminosCount > 0){
			particleCount = particleCount + sop.additionalAminosCount;
		}
		//printf("Coordinates will be saved as '%s'.\n", dcd_filenames[traj]);
		dcd.open_write(dcd_filenames[traj]);
        dcd.N = particleCount;
        dcd.NFILE = dcd.NPRIV = dcd.NSAVC = dcdOutputManager.frequency;
        dcd.DELTA = 0.0; // The integrator is not initialized yet, so we consider timestep to be equal to zero
		dcd.write_header();
		dcd.close();
	}
	// Allocate memory for temporary data
	int size = sop.aminoCount*sizeof(float);
	if(sop.additionalAminosCount > 0){
		size = size + sop.additionalAminosCount*sizeof(float);
	}
	X = (float*) malloc(size);
	Y = (float*) malloc(size);
	Z = (float*) malloc(size);
	printf("Coordinates will be saved in '%s'.\n", dcd_filename);
}

/*
 * Saving coordinates to DCD
 */
void saveCoord(){
	int i, traj;
	if(step % dcdOutputManager.frequency == 0){
		printf("Saving coordinates into dcd...");
		copyCoordDeviceToHost();
		int particleCount = sop.aminoCount;
		for(traj = 0; traj < gsop.Ntr; traj++){
			for(i = 0; i < sop.aminoCount; i++){
				X[i] = gsop.h_coord[sop.aminoCount*traj + i].x;
				Y[i] = gsop.h_coord[sop.aminoCount*traj + i].y;
				Z[i] = gsop.h_coord[sop.aminoCount*traj + i].z;
			}
			if(sop.additionalAminosCount > 0){
				particleCount = particleCount + sop.additionalAminosCount;
				for(i = 0; i < sop.additionalAminosCount; i++){
					X[i+sop.aminoCount] = sop.additionalAminos[i].x;
					Y[i+sop.aminoCount] = sop.additionalAminos[i].y;
					Z[i+sop.aminoCount] = sop.additionalAminos[i].z;
				}
			}
			dcd.open_append(dcd_filenames[traj]);
			dcd.write_frame(X, Y, Z);
			dcd.close();
		}
		printf("done.\n");
		long long int timer = time(NULL) - initialTime;
		int days = timer / (3600*24);
		int hours = timer / 3600 - days * 24;
		int minutes = timer / 60 - hours * 60 - days * 24 * 60;
		int seconds = timer - hours * 3600 - days * 24 * 3600 - minutes * 60;
		printf("Computation time: %dd %dh %dm %ds (%f s/step current, %f ms/step overall)\n", days, hours, minutes, seconds, ((float)time(NULL) - lastTime)/((float)restartfreq), 1000.0f*((float)timer)/((float)(step)));
		if(step != 0){
			timer = ((time(NULL)-initialTime) * (numsteps - step)) / step;
			days = timer / (3600*24);
			hours = timer / 3600 - days * 24;
			minutes = timer / 60 - hours * 60 - days * 24 * 60;
			printf("Estimated time left: %dd %dh %dm\n", days, hours, minutes);
		}
		lastTime = time(NULL); // Timer stop-line
	}

	if(step % restartfreq == 0){
		for(int traj = 0; traj < gsop.Ntr; traj++){
			char trajnum[10];
			char tempRestartFilename[100];
			sprintf(trajnum, "%d", traj+gsop.firstrun);
			replaceString(tempRestartFilename, restartpdb_filename, trajnum, "<run>");
			for(i = 0; i < sop.aminoCount; i++){
				sop.aminos[i].x = gsop.h_coord[sop.aminoCount*traj + i].x;
				sop.aminos[i].y = gsop.h_coord[sop.aminoCount*traj + i].y;
				sop.aminos[i].z = gsop.h_coord[sop.aminoCount*traj + i].z;
			}
			savePDB(tempRestartFilename, sop);
			/*replaceString(tempRestartFilename, restartconf_filename, trajnum, "<run>");
			FILE* param_file = fopen(tempRestartFilename, "w");
			fprintf(param_file, "stage %s\n", stageString);
			fprintf(param_file, "timeStep %ld\n", step);
			fprintf(param_file, "extension %f\n", xt);
			fprintf(param_file, "pullVector.x %f\npullVector.y %f\npullVector.z %f\n",
						pull_vector[traj].x, pull_vector[traj].y, pull_vector[traj].z);
			fprintf(param_file, "chipCoord.x %f\nchipCoord.y %f\nchipCoord.z %f\n",
						cantilever_coord[traj].x, cantilever_coord[traj].y, cantilever_coord[traj].z);
			fclose(param_file);*/
		}
		printf("Saving restart coordinates into '%s'.\n", restartpdb_filename);
	//	printf("Saving restart parameters into '%s'.\n", restartconf_filename);
	}
	checkCUDAError();
}

