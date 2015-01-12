/*
 * dcd_manager.cu
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#include "dcd_manager.h"
#include "../IO/dcdio.h"
#include "../IO/pdbio.h"

std::vector<std::string> dcd_filenames;
std::vector<std::string> restart_filenames;

int restartfreq;

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
	this->frequency = parameters::dcdfreq.get();
	restartfreq = parameters::restartfreq.get();

	int traj;
	dcd_filenames.resize(gsop.Ntr);
	restart_filenames.resize(gsop.Ntr);
	for(traj = 0; traj < gsop.Ntr; traj++){
        dcd_filenames[traj] = parameters::DCDfile.replace("<run>", traj);
        restart_filenames[traj] = parameters::restartname.replace("<run>", traj);
		int particleCount = sop.aminoCount;

		if(sop.additionalAminosCount > 0){
			particleCount = particleCount + sop.additionalAminosCount;
		}
		//printf("Coordinates will be saved as '%s'.\n", dcd_filenames[traj]);
		dcd.open_write(dcd_filenames[traj].c_str());
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
	printf("Coordinates will be saved in '%s'.\n", parameters::DCDfile.get().c_str());
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
			dcd.open_append(dcd_filenames[traj].c_str());
			dcd.write_frame(X, Y, Z);
			dcd.close();
		}
		printf("done.\n");
	}

	if(gsop.step % restartfreq == 0){
		for(int traj = 0; traj < gsop.Ntr; traj++){
			for(i = 0; i < sop.aminoCount; i++){
				sop.aminos[i].x = gsop.h_coord[sop.aminoCount*traj + i].x;
				sop.aminos[i].y = gsop.h_coord[sop.aminoCount*traj + i].y;
				sop.aminos[i].z = gsop.h_coord[sop.aminoCount*traj + i].z;
			}
			savePDB(restart_filenames[traj].c_str(), sop);
		}
		printf("Saving restart coordinates into '%s'.\n", parameters::restartname.get().c_str());
	}
	checkCUDAError();
}

