/*
 * dcd_manager.cu
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#include "dcd_manager.h"
#include "../IO/dcdio.h"
#include "../IO/pdbio.h"
#include "../IO/configreader.h"

char dcd_filename[100];
char restartpdb_filename[100];
int restartfreq;

char** dcd_filenames;
DCD dcd;
float* X;
float* Y;
float* Z;

void createDCDOutputManager(){
	updaters[updatersCount] = new DcdOutputManager();
	updatersCount++;
}

/*
 * Initialization of coordinates output
 */
DcdOutputManager::DcdOutputManager(){
	printf("Initializing dcd output manager...\n");
	this->name = "DCD output";
	this->frequency = getIntegerParameter(DCD_FREQUENCY_STRING, 10000, 1);;
	restartfreq = getIntegerParameter("restartfreq", 100000, 1);

	getMaskedParameter(dcd_filename, "DCDfile", "<name>_<author><run>_<stage>.dcd", 1);
	//printf("Coordinates will be saved into '%s'\n", dcd_filename);
	char restart_filename[100];
	getMaskedParameter(restart_filename, "restartname", "<name>_<author><run>_restart", 1);
	sprintf(restartpdb_filename, "%s.pdb", restart_filename);

	int traj;
	dcd_filenames = (char**)calloc(gsop.Ntr, sizeof(char*));
	for(traj = 0; traj < gsop.Ntr; traj++){
		dcd_filenames[traj] = (char*)calloc(100, sizeof(char));
		char trajnum[10];
		sprintf(trajnum, "%d", traj+gsop.firstrun);
		replaceString(dcd_filenames[traj], dcd_filename, trajnum, "<run>");
		int particleCount = sop.aminoCount;
		if(sop.additionalAminosCount > 0){
			particleCount = particleCount + sop.additionalAminosCount;
		}
		//printf("Coordinates will be saved as '%s'.\n", dcd_filenames[traj]);
		dcd.open_write(dcd_filenames[traj]);
        dcd.N = particleCount;
        dcd.NFILE = dcd.NPRIV = dcd.NSAVC = this->frequency;
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
	printf("Done initializing dcd output manager...\n");
}

/*
 * Saving coordinates to DCD
 */
void DcdOutputManager::update(){
	int i, traj;
	if(gsop.step % this->frequency == 0){
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
	}

	if(gsop.step % restartfreq == 0){
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
		}
		printf("Saving restart coordinates into '%s'.\n", restartpdb_filename);
	}
	checkCUDAError();
}

