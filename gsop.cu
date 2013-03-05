/*
 * gsop.cu
 *
 *  Created on: Jun 4, 2009
 *      Author: zhmurov
 */
#include "def_param.h"
#include "gsop.cuh"
#include "Potentials/covalent.cu"
#include "Potentials/native.cu"
#include "Potentials/pairs.cu"
#include "Potentials/indentation.cu"
#include "Potentials/pulling.cu"
#include "Updaters/pairlist.cu"
#include "Updaters/possiblepairlist.cu"
#include "Updaters/output_manager.cu"
#include "Updaters/dcd_manager.cu"
#include "Integrators/langevin.cu"
//#include "externalForce.cu"

void initGPU();
void initFF();
void runGPU();
void bindTextures();
void checkCUDAError();

extern void savePDB(char* pdb_filename);

void initCoordinates();
void copyCoordinates();
void copyCoordinatesTrajectory(int traj);
void initForces();
void initEnergies();

extern void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace);
extern void initTemperature();
extern void initCovalent();
extern void initNative();
extern void initPairs();
extern uint4* initRandom(int seed, int N);
extern void shiftRandomSeeds(uint4* d_newseeds, int N);

unsigned int cpuTimer;
double gpuTime = 0.0;
double cpuTime = 0.0;
extern double pairlistTime = 0.0;
extern double covalentTime = 0.0;
extern double nativeTime = 0.0;
extern double pairsTime = 0.0;

long int lastStepCoordCopied = -1;

GSOP gsop;
SOPPotential** potentials;
SOPUpdater** updaters;
SOPIntegrator* integrator;

int potentialsCount;
int updatersCount;

/*
 * Prepare data on Host, copy it into Device, bind textures
 */
void initGPU(){

	//initDCD();
	cudaSetDevice(device);
	cudaGetDeviceProperties(&gsop.deviceProp, device);
	printf("Using device %d: \"%s\"\n", device, gsop.deviceProp.name);
	printf("CUDA Capability revision number: %d.%d\n", gsop.deviceProp.major, gsop.deviceProp.minor);
	gsop.aminoCount = Ntr*sop.aminoCount;
	gsop.width = gsop.aminoCount;
	while(gsop.width % 8 != 0){
		gsop.width++;
	}
	printf("Will align structures to width of %d.\n", gsop.width);

	initCoordinates(); // Alocate memory for coordinates
	initForces(); // Allocate memory for forces
}

void initFF(){

	potentialsCount = 0;
	updatersCount = 0;
	gsop.minimizationOn = minimizationOn;
	gsop.pullingOn = pullingOn;

	// Allocating memory for the model
	int i;
	potentials = (SOPPotential**)calloc(max_potentials, sizeof(SOPPotential*));
	for(i = 0; i < max_potentials; i++){
		potentials[i] = (SOPPotential*)malloc(sizeof(SOPPotential));
	}
	updaters = (SOPUpdater**)calloc(max_updaters, sizeof(SOPUpdater*));
	for(i = 0; i < max_updaters; i++){
		updaters[i] = (SOPUpdater*)malloc(sizeof(SOPUpdater));
	}
	integrator = (SOPIntegrator*)malloc(sizeof(SOPIntegrator));

	// Creating model
	createCovalentPotential(); // FENE
	createNativePotential(); // Full LJ
	createPairsPotential(); // Repulsive LJ
	createIndentationPotential(); // Indentations sphere and surface
	createPullingPotential(); // External force

	if(Ntr == 1){
		createPossiblepairlistUpdater(); // Updates the list of all pairs (for Verlet list)
	}
	createOutputManager(); // Save dat output
	createPairlistUpdater(); // Verlet list
	createDCDOutputManager(); // Save coordinates (dcd + pdb restart)

	createLangevinIntegrator(); // Create integrator

	initEnergies(); // Allocate memory for energy output (move to initGPU() ?)

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
	int stride = nav;
	for(u = 0; u < updatersCount; u ++){
		if(updaters[u]->frequency < stride){
			stride = updaters[u]->frequency;
		}
	}
	int size = gsop.aminoCount*sizeof(float4);
	generatePairlist();

	// External loop
	for(i = step/stride; i <= numsteps/stride; i++){
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
		step += stride;

		checkCUDAError();
	}
	copyCoordDeviceToHost();


	int traj;
	for(traj = 0; traj < Ntr; traj++){
		char trajnum[10];
		char trajCoordFilename[100];
		sprintf(trajnum, "%d\0", traj+firstrun);
		replaceString(trajCoordFilename, final_filename, trajnum, "<run>");
		for(i = 0; i < sop.aminoCount; i++){
			sop.aminos[i].x = gsop.h_coord[sop.aminoCount*traj + i].x;
			sop.aminos[i].y = gsop.h_coord[sop.aminoCount*traj + i].y;
			sop.aminos[i].z = gsop.h_coord[sop.aminoCount*traj + i].z;
		}
		savePDB(trajCoordFilename);
	}

}

void initCoordinates(){
	printf("Initializing coordinates (%d particles)...\n", gsop.aminoCount);
	int size = gsop.aminoCount*sizeof(float4);
	//gsop.h_coord = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMallocHost((void**)&gsop.h_coord, size);
	cudaMalloc((void**)&gsop.d_coord, size);
	//cudaMalloc((void**)&gsop.d_coordToSave, size);
}

void copyCoordinates(){
	printf("Copying coordinates...\n");
	int size = gsop.aminoCount*sizeof(float4);
	cudaMemcpy(gsop.d_coord, gsop.h_coord, size, cudaMemcpyHostToDevice);
	//cudaMemcpy(gsop.d_coordToSave, gsop.h_coord, size, cudaMemcpyHostToDevice);
}

void copyCoordDeviceToHost(){
	if(step != lastStepCoordCopied){
		cudaMemcpy(gsop.h_coord, gsop.d_coord, gsop.aminoCount*sizeof(float4), cudaMemcpyDeviceToHost);
		lastStepCoordCopied = step;
	}
}

void copyCoordinatesTrajectory(int traj){
	int i;
	for(i = 0; i < sop.aminoCount; i++){
		gsop.h_coord[traj*sop.aminoCount + i].x = sop.aminos[i].x;
		gsop.h_coord[traj*sop.aminoCount + i].y = sop.aminos[i].y;
		gsop.h_coord[traj*sop.aminoCount + i].z = sop.aminos[i].z;
		gsop.h_coord[traj*sop.aminoCount + i].w = 0;
	}
#ifdef DEBUG
	printf("Coordinates for run #%d:\n", traj+firstrun);
	for(i = 0; i < sop.aminoCount; i++){
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
 * Initialize array of energies for output
 */
void initEnergies(){
	// Allocating memory
	printf("Allocating memory for energies...\n");
	int size = gsop.aminoCount*sizeof(float4);
	//gsop.h_energies = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMallocHost((void**)&gsop.h_energies, size);
	int i;
	for(i = 0; i < gsop.aminoCount; i++){
		gsop.h_energies[i].x = 0.0f;
		gsop.h_energies[i].y = 0.0f;
		gsop.h_energies[i].z = 0.0f;
		gsop.h_energies[i].w = 0.0f;
	}
	cudaMalloc((void**)&gsop.d_energies, size);
	//cudaMalloc((void**)&gsop.d_energiesToSave, size);
	// Copying to the Device
	cudaMemcpy(gsop.d_energies, gsop.h_energies, size, cudaMemcpyHostToDevice);
	//cudaMemcpy(gsop.d_energiesToSave, gsop.h_energies, size, cudaMemcpyHostToDevice);
}

void bindTextures(){
#ifndef NOTEXTURE
	cudaBindTexture(0, t_coord, gsop.d_coord, gsop.aminoCount*sizeof(float4));
#endif
}

void checkCUDAError(){
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess){
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		exit(0);
	}
}
