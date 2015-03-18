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
		int particleCount = sop.aminos.size();

		if(sop.additionalAminos.size() > 0){
			particleCount = particleCount + sop.additionalAminos.size();
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
	int size = sop.aminos.size()*sizeof(float);
	if(sop.additionalAminos.size() > 0){
		size = size + sop.additionalAminos.size()*sizeof(float);
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
		int particleCount = sop.aminos.size();
		for(traj = 0; traj < gsop.Ntr; traj++){
			for(i = 0; i < sop.aminos.size(); i++){
				X[i] = gsop.h_coord[sop.aminos.size()*traj + i].x;
				Y[i] = gsop.h_coord[sop.aminos.size()*traj + i].y;
				Z[i] = gsop.h_coord[sop.aminos.size()*traj + i].z;
			}
			if(sop.additionalAminos.size() > 0){
				particleCount = particleCount + sop.additionalAminos.size();
				for(i = 0; i < sop.additionalAminos.size(); i++){
					X[i+sop.aminos.size()] = sop.additionalAminos[i].x;
					Y[i+sop.aminos.size()] = sop.additionalAminos[i].y;
					Z[i+sop.aminos.size()] = sop.additionalAminos[i].z;
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
			for(i = 0; i < sop.aminos.size(); i++){
				sop.aminos[i].x = gsop.h_coord[sop.aminos.size()*traj + i].x;
				sop.aminos[i].y = gsop.h_coord[sop.aminos.size()*traj + i].y;
				sop.aminos[i].z = gsop.h_coord[sop.aminos.size()*traj + i].z;
			}
			savePDB(restart_filenames[traj].c_str(), sop);
		}
		printf("Saving restart coordinates into '%s'.\n", parameters::restartname.get().c_str());
	}
	checkCUDAError();
}

