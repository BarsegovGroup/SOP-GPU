/*
 * pulling_plane.cu
 *
 *  Created on: Jan 17, 2014
 *      Author: alekseenko
 */
#include "../gsop.cuh"
#include "../Integrators/langevin.h"
#include "../IO/configreader.h"
#include "pulling_plane.h"

char pullingPlaneFilename[500];
FILE* pullingPlaneFile;

#include "pulling_plane_kernel.cu"

void createPullingPlanePotential(){
	if(getYesNoParameter(PULLINGPLANE_ON_STRING, 0, 1) == 1){
		gsop.pullingPlaneOn = 1;

        PullingPlanePotential *pot;
		potentials[potentialsCount] = pot = new PullingPlanePotential();
		potentialsCount++;

		updaters[updatersCount] = new PullingPlaneUpdater(pot);
		updatersCount++;
	}
}

PullingPlaneUpdater::PullingPlaneUpdater(PullingPlanePotential *pullingPlane){
	this->name = "Pulling Plane";
    this->frequency = getIntegerParameter(PULLINGPLANE_FREQ, nav, 1);
    this->pullingPlane = pullingPlane;
}

PullingPlanePotential::PullingPlanePotential(){
    this->name = "Pulling Plane";
	printf("Initializing pulling plane protocol...\n");
    if (gsop.Ntr != 1) {
        DIE("Pulling plane can only run in single-trajectory-per-GPU mode (runnum 1)\n");
    }
	int i, j;
	this->h_extForces = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMalloc((void**)&this->d_extForces, gsop.aminoCount*sizeof(float4));

	this->deltax = getFloatParameter(PULLINGPLANE_DELTAX_STRING, 0, 1);
	this->Ks = getFloatParameter(PULLINGPLANE_KS_STRING, DEFAULT_PULLINGPLANE_KS, 1);

	if(this->deltax == 0){
		printf("ERROR: '%s' parameter should be specified to initiate this->\n", PULLINGPLANE_DELTAX_STRING);
		exit(-1);
	}

	this->fixedCount = getIntegerParameter(PULLINGPLANE_FIXED_COUNT_STRING, 0, 0);
	this->fixed = (int*)malloc(this->fixedCount*sizeof(int));
	this->pulledCount = getIntegerParameter(PULLINGPLANE_PULLED_COUNT_STRING, 0, 0);
	this->pulled = (int*)malloc(this->pulledCount*sizeof(int));
	printf("%d resid(s) fixed, %d pulled.\n", this->fixedCount, this->pulledCount);
	char paramName[50];
	for(i = 0; i < this->fixedCount; i++){
		sprintf(paramName, "%s%d", PULLINGPLANE_FIXED_STRING, i+1);
		this->fixed[i] = getIntegerParameter(paramName, 0, 0);
		if(this->fixed[i] < 0 || this->fixed[i] >= gsop.aminoCount){
			DIE("ERROR: Fixed bead %s %d not exists. Protein has only %d amino-acids. Bead numbers should start with zero.\n", paramName, this->fixed[i], gsop.aminoCount);
		}
		printf("Resid %d is fixed.\n", this->fixed[i]);
	}
	for(i = 0; i < this->pulledCount; i++){
		sprintf(paramName, "%s%d", PULLINGPLANE_PULLED_STRING, i+1);
		this->pulled[i] = getIntegerParameter(paramName, 0, 0);
		if(this->pulled[i] < 0 || this->pulled[i] >= gsop.aminoCount){
			DIE("ERROR: Pulled bead %s %d not exists. Protein has only %d amino-acids. Bead numbers should start with zero.\n", paramName, this->pulled[i], gsop.aminoCount);
		}
		printf("Pulling resid %d.\n", this->pulled[i]);
	}

    getVectorParameter(PULLINGPLANE_PULLVECTOR_STRING, &this->pullVector.x, &this->pullVector.y, &this->pullVector.z, 0, 0, 0, 0);
    double t = sqrt(this->pullVector.x*this->pullVector.x + this->pullVector.y*this->pullVector.y + this->pullVector.z*this->pullVector.z);
    this->pullVector.x /= t;
    this->pullVector.y /= t;
    this->pullVector.z /= t;
    getVectorParameter(PULLINGPLANE_ZEROVECTOR_STRING, &this->planeCoord0.x, &this->planeCoord0.y, &this->planeCoord0.z, 0, 0, 0, 0);
    this->d = - (this->planeCoord.x*this->pullVector.x + this->planeCoord.y*this->pullVector.y + this->planeCoord.z*this->pullVector.z);
    this->cant_d = this->d;

	for(i = 0; i < sop.aminoCount; i++){
		this->h_extForces[i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
	}
	for(j = 0; j < this->pulledCount; j++){
		i = this->pulled[j];
		printf("Pulling bead #%d (%s %d chain %c).\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
		this->h_extForces[i] = make_float4(0.0, 0.0, 0.0, 2.0);
	}
	for(j = 0; j < this->fixedCount; j++){
		i = this->fixed[j];
		printf("Fixing bead #%d (%s %d chain %c).\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
		this->h_extForces[i] = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
	}

	for(i = 0; i < sop.aminoCount; i++){
		sop.aminos[i].beta = this->h_extForces[i].w;
	}
	cudaMemcpy(this->d_extForces, this->h_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
    checkCUDAError();

    hc_pullingPlane.d_extForces = this->d_extForces;
    hc_pullingPlane.pullVector = this->pullVector;
    hc_pullingPlane.d = this->d;
    hc_pullingPlane.Ks = this->Ks;
	cudaMemcpyToSymbol(c_pullingPlane, &hc_pullingPlane, sizeof(PullingPlaneConstant), 0, cudaMemcpyHostToDevice);
	checkCUDAError();

    char tmpstr[512];
	getMaskedParameter(tmpstr, PULLINGPLANE_FILENAME, "", 0);
	char trajnum[10];
	sprintf(trajnum, "%d", gsop.firstrun);
	replaceString(pullingPlaneFilename, tmpstr, trajnum, "<run>");
	pullingPlaneFile = fopen(pullingPlaneFilename, "w");
	fclose(pullingPlaneFile);
	printf("PullingPlane data will be saved in '%s'.\n", pullingPlaneFilename);
	
    if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(pullingPlane_kernel, cudaFuncCachePreferL1);
	}

    this->blockSize = getIntegerParameter(COVALENT_BLOCK_SIZE_STRING, gsop.blockSize, 1);
	this->blockNum = gsop.aminoCount/this->blockSize + 1;

	printf("Done initializing pulling plane protocol...\n");
}

void PullingPlanePotential::compute(){
	pullingPlane_kernel<<<this->blockNum, this->blockSize>>>();
	checkCUDAError();
}

void PullingPlanePotential::computeEnergy(){
    // Perl has operator "..."
    // It should be used inside functions to be implemented
    // Unlike just empty functions, the warning is generated when "..." is used
    // This comment is completely useless and barely related, because this function is empty by design
}

void PullingPlaneUpdater::update(){

	//copyCoordDeviceToHost();
	int j;
	// Increasing the force
	float xt = pullingPlane->deltax * step / this->frequency;
    /*
	pullingPlane->cantCoord.x = pullingPlane->planeCoord0.x + xt * pullingPlane->pullVector.x;
	pullingPlane->cantCoord.y = pullingPlane->planeCoord0.y + xt * pullingPlane->pullVector.y;
	pullingPlane->cantCoord.z = pullingPlane->planeCoord0.z + xt * pullingPlane->pullVector.z;
    */
    pullingPlane->cant_d = pullingPlane->d0 + xt;
	checkCUDAError();
	cudaMemcpy(pullingPlane->h_extForces, pullingPlane->d_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyDeviceToHost);
	checkCUDAError();
	if(step % this->frequency == 0){
		if(step % 100000 == 0){
			printf("Cantilever coordinates for run #%d: %f, %f, %f\n",
					gsop.firstrun, pullingPlane->cantCoord.x, pullingPlane->cantCoord.y, pullingPlane->cantCoord.z);
			printf("Plane coordinates for run #%d: %f, %f, %f\n",
					gsop.firstrun, pullingPlane->planeCoord.x, pullingPlane->planeCoord.y, pullingPlane->planeCoord.z);
		}
        pullingPlaneFile = fopen(pullingPlaneFilename, "a");
        float3 extForce = make_float3(0.f, 0.f, 0.f);
        for (j = 0; j < gsop.aminoCount; j++){
            extForce.x += pullingPlane->h_extForces[j].x;
            extForce.y += pullingPlane->h_extForces[j].y;
            extForce.z += pullingPlane->h_extForces[j].z;
        }

        float extForceProj = extForce.x*pullingPlane->pullVector.x + extForce.y*pullingPlane->pullVector.y + extForce.z*pullingPlane->pullVector.z;

        float totForce = pullingPlane->Ks * (pullingPlane->cant_d - pullingPlane->d) - extForceProj;

        pullingPlane->d += totForce * this->frequency * ((LangevinIntegrator*)integrator)->hOverZeta; // TODO: fix dependency

        hc_pullingPlane.d = pullingPlane->d;
        cudaMemcpyToSymbol(c_pullingPlane, &hc_pullingPlane, sizeof(PullingPlaneConstant), 0, cudaMemcpyHostToDevice);

        fprintf(pullingPlaneFile, "%12ld\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",
                step, pullingPlane->d, pullingPlane->cant_d, extForceProj,
                extForce.x, extForce.y, extForce.z);

        fflush(pullingPlaneFile);
        fclose(pullingPlaneFile);
    }

	checkCUDAError();
}

