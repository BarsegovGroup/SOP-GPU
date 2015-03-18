/*
 * gsop.cu
 *
 *  Created on: Jun 4, 2009
 *      Author: zhmurov
 */

#include "gsop.cuh"

#include "IO/topio.h"
#include "IO/dcdio.h"
#include "IO/pdbio.h"
#include "Util/wrapper.h"
#include "Util/parameters.h"

#include "Potentials/covalent.cu"
#include "Potentials/native.cu"
#include "Potentials/pairs.cu"
#include "Potentials/indentation.cu"
#include "Potentials/pulling.cu"
#include "Potentials/pulling_plane.cu"
#include "Integrators/langevin.cu"
#include "Integrators/bdhitea.cu"
#include "Updaters/pairlist.cu"
#include "Updaters/possiblepairlist.cu"
#include "Updaters/output_manager.cu"
#include "Updaters/dcd_manager.h"

long int lastStepCoordCopied = -1;

GSOP gsop;
SOPPotential** potentials;
SOPUpdater** updaters;
SOPIntegrator* integrator;

cudaDeviceProp deviceProp;

int potentialsCount;
int updatersCount;
int integratorTea;

/*
 * Prepare data on Host, copy it into Device, bind textures
 */
void initGPU(){

	//initDCD();
	cudaSetDevice(gsop.deviceId);
	cudaGetDeviceProperties(&deviceProp, gsop.deviceId);
	printf("Using device %d: \"%s\"\n", gsop.deviceId, deviceProp.name);
	printf("CUDA Capability revision number: %d.%d\n", deviceProp.major, deviceProp.minor);
	printf("CUDA PCIe ID: %x:%x:%x\n", deviceProp.pciDomainID, deviceProp.pciBusID, deviceProp.pciDeviceID);
	gsop.aminoCount = gsop.Ntr*sop.aminos.size();
	gsop.width = gsop.aminoCount; // It's only used for mica lists in indentation; do we really need it?
	while(gsop.width % 8 != 0){
		gsop.width++;
	}
	printf("Will align structures to width of %d.\n", gsop.width);
	gsop.blockSize = parameters::block_size.get();

	initCoordinates(); // Allocate memory for coordinates
	initForces(); // Allocate memory for forces
	initTemperature(); //Allocate memory for temperature
}

void initFF(){

	potentialsCount = 0;
	updatersCount = 0;

    std::string stageString = parameters::stage.get();
	gsop.minimizationOn = (stageString == "minim" );
	gsop.heatingOn      = (stageString == "heat"  );
	gsop.pullingOn      = (stageString == "pull"  );
	gsop.indentationOn  = (stageString == "indent");

	// Allocating memory for the model
	potentials = (SOPPotential**)calloc(max_potentials, sizeof(SOPPotential*));
	updaters = (SOPUpdater**)calloc(max_updaters, sizeof(SOPUpdater*));

	// Creating model
	createCovalentPotential(); // FENE
	createNativePotential(); // Full LJ
	createPairsPotential(); // Repulsive LJ
	createIndentationPotential(); // Indentations sphere and surface
	createPullingPotential(); // External force
	createPullingPlanePotential();

	if(gsop.Ntr == 1){
		createPossiblepairlistUpdater(); // Updates the list of all pairs (for Verlet list)
	}
	createOutputManager(); // Save dat output
	createPairlistUpdater(); // Verlet list
	createDCDOutputManager(); // Save coordinates (dcd + pdb restart)

	// Create integrator
	if (parameters::tea_on.get()){
		createTeaIntegrator();
		integratorTea = 1;
	}else{
		createLangevinIntegrator();
		integratorTea = 0;
	}

	cudaMemcpyToSymbol(c_gsop, &gsop, sizeof(GSOP), 0, cudaMemcpyHostToDevice);
	bindTextures();
	checkCUDAError();
}

/*
 *
 */
void runGPU(){
	printf("Starting simulations.\n");
	/*if(stage == pull_stage){
		engageCantileverTip();
	}*/

	int i, j, p, u;

	// Leave internal loop only when updater execution is needed (most frequent updater)
	int stride = gsop.nav;
	for(u = 0; u < updatersCount; u ++){
		if(updaters[u]->frequency < stride){
			stride = updaters[u]->frequency;
		}
	}
    // TODO: WTF is this: vvvv
	// generatePairlist(); // WUT? Isn't it called in the loop below?
    // I don't know why was it necessary, so I will recreate it for now
    updaterByName("Pairlist")->update();

	// External loop
	for(i = gsop.step/stride; i <= gsop.numsteps/stride; i++){
		// Run all updaters
		for(u = 0; u < updatersCount; u++){
			updaters[u]->update();
		}
		// Internal loop
		for(j = 0; j < stride; j++){
			// Compute all potentials
			for(p = 0; p < potentialsCount; p++){
				potentials[p]->compute();
			}
			integrator->integrate(); // Integrate
			checkCUDAError();
		}

		//cudaMemcpy(gsop.h_coord, gsop.d_coord, size, cudaMemcpyDeviceToHost);
		gsop.step += stride;

		checkCUDAError();
	}
	copyCoordDeviceToHost();


	int traj;
	for(traj = 0; traj < gsop.Ntr; traj++){
        std::string trajCoordFilename = parameters::finalcoord.replace("<run>", traj+gsop.firstrun);
		for(i = 0; i < sop.aminos.size(); i++){
			sop.aminos[i].x = gsop.h_coord[sop.aminos.size()*traj + i].x;
			sop.aminos[i].y = gsop.h_coord[sop.aminos.size()*traj + i].y;
			sop.aminos[i].z = gsop.h_coord[sop.aminos.size()*traj + i].z;
		}
		savePDB(trajCoordFilename.c_str(), sop);
	}

}

void initCoordinates(){
	printf("Initializing coordinates (%d particles)...\n", gsop.aminoCount);
	int size = gsop.aminoCount*sizeof(float4);
	cudaMallocHost((void**)&gsop.h_coord, size);
	cudaMalloc((void**)&gsop.d_coord, size);
}

void copyCoordinates(){
	printf("Copying coordinates...\n");
	int size = gsop.aminoCount*sizeof(float4);
	cudaMemcpy(gsop.d_coord, gsop.h_coord, size, cudaMemcpyHostToDevice);
}

void copyCoordDeviceToHost(){
	if(gsop.step != lastStepCoordCopied){
		cudaMemcpy(gsop.h_coord, gsop.d_coord, gsop.aminoCount*sizeof(float4), cudaMemcpyDeviceToHost);
		lastStepCoordCopied = gsop.step;
	}
}

void copyCoordinatesTrajectory(int traj){
	int i;
	for(i = 0; i < sop.aminos.size(); i++){
		gsop.h_coord[traj*sop.aminos.size() + i].x = sop.aminos[i].x;
		gsop.h_coord[traj*sop.aminos.size() + i].y = sop.aminos[i].y;
		gsop.h_coord[traj*sop.aminos.size() + i].z = sop.aminos[i].z;
		gsop.h_coord[traj*sop.aminos.size() + i].w = 0;
	}
#ifdef DEBUG
	printf("Coordinates for run #%d:\n", traj+gsop.firstrun);
	for(i = 0; i < sop.aminos.size(); i++){
		printf("%d:\t%f\t%f\t%f\n", i, sop.aminos[i].x, sop.aminos[i].y, sop.aminos[i].z);
	}
#endif
}

/*
 * Initialize array of forces
 */
void initForces(){
	// Allocating memory
	printf("Allocating memory for forces...\n");
	int size = gsop.aminoCount*sizeof(float4);
	gsop.h_forces = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMalloc((void**)&gsop.d_forces, size);
	// Copying to the Device
	cudaMemcpy(gsop.d_forces, gsop.h_forces, size, cudaMemcpyHostToDevice);
}

/*
 * Initialize array for temperature output
 */
void initTemperature(){
	// Allocating memory
	printf("Allocating memory for energies...\n");
	int size = gsop.aminoCount*sizeof(float);
	//gsop.h_energies = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMallocHost((void**)&gsop.h_T, size);
	int i;
	for(i = 0; i < gsop.aminoCount; i++){
		gsop.h_T[i] = 0.0f;
	}
	cudaMalloc((void**)&gsop.d_T, size);
	// Copying to the Device
	cudaMemcpy(gsop.d_T, gsop.h_T, size, cudaMemcpyHostToDevice);
}

void copyTemperatureDeviceToHost(){
	int size = gsop.aminoCount*sizeof(float);
	cudaMemcpy(gsop.h_T, gsop.d_T, size, cudaMemcpyDeviceToHost);
}

SOPPotential* potentialByName(const char *name){
    int i;
    for(i = 0; i < potentialsCount; i++){
        if (potentials[i]->name == name)
            return potentials[i];
    }
    DIE("Unable to find requested potential '%s'", name);
}

SOPUpdater* updaterByName(const char *name){
    int i;
    for(i = 0; i < updatersCount; i++){
        if (updaters[i]->name == name)
            return updaters[i];
    }
    DIE("Unable to find requested updater '%s'", name);
}

void SOPPotential::sumEnergies(const float *h_energies, float *energies) {
	int traj, i;
	for(traj = 0; traj < gsop.Ntr; traj++){
		double energy = 0.0;
		for(i = 0; i < sop.aminos.size(); i++){
			int index = traj*sop.aminos.size() + i;
			energy += h_energies[index];
			if(isinf(h_energies[index])){
				printf("WARNING: Bead #%d in trajectory #%d has infinite %s energy\n", i, traj, this->name.c_str());
			}
			if(isnan(h_energies[index])){
				printf("WARNING: The %s energy of bead #%d in trajectory #%d is NaN\n", this->name.c_str(), i, traj);
			}
		}
		energies[traj] = energy*0.5f;
	}
}

void bindTextures(){
#ifndef NOTEXTURE
	cudaBindTexture(0, t_coord, gsop.d_coord, gsop.aminoCount*sizeof(float4));
#endif
}

void __checkCUDAError(const char *file, int line){
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess){
		DIE("CUDA error: %s [%s:%d]\n", cudaGetErrorString(error), file, line);
	}
}

