/*
 * pulling.cu
 *
 *  Created on: May 25, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "../IO/configreader.h"
#include "../Util/mystl.h"
#include "../Util/wrapper.h"
#include "../Util/vector_helpers.h"
#include "pulling.h"

#include "pulling_kernel.cu"

void createPullingPotential(){
	if(gsop.pullingOn == 1 || getYesNoParameter(PULLING_ON_STRING, 0, 1) == 1){
		gsop.pullingOn = 1;

        PullingPotential *pot;
		potentials[potentialsCount] = pot = new PullingPotential();
        potentialsCount++;

		if(getFloatParameter(PULLING_DELTAX_STRING, 0, 1) != 0){
			updaters[updatersCount] = new PullingUpdater(pot);
			updatersCount++;
		}
	}
}

PullingUpdater::PullingUpdater(PullingPotential *pulling){
	this->name = "Pulling";
    this->frequency = getIntegerParameter(PULLING_FREQ, gsop.nav, 1);
    this->pulling = pulling;
}

PullingPotential::PullingPotential(){
    this->name = "Pulling";
	printf("Initializing pulling protocol...\n");
	int i, j, traj;
	this->pullVector = (float3*)calloc(gsop.Ntr, sizeof(float3));
	this->cantCoord0 = (float3*)calloc(gsop.Ntr, sizeof(float3));
	this->cantCoord = (float3*)calloc(gsop.Ntr, sizeof(float3));
	this->extForce = (float3*)calloc(gsop.Ntr, sizeof(float3));
	this->h_extForces = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMalloc((void**)&this->d_extForces, gsop.aminoCount*sizeof(float4));

	this->deltax = getFloatParameter(PULLING_DELTAX_STRING, 0, 1);
	this->Ks = getFloatParameter(PULLING_KS_STRING, DEFAULT_PULLING_KS, 1);
	this->fconst = getFloatParameter(PULLING_FCONST_STRING, 0, 1);

	if(this->deltax == 0 && this->fconst == 0){
		DIE("ERROR: Either 'deltax' or 'fconst' parameter should be specified to initiate pulling\n");
	}

	this->fixedCount = getIntegerParameter(PULLING_FIXED_COUNT_STRING, 0, 0);
	this->fixed = (int*)malloc(this->fixedCount*sizeof(int));
	this->pulledCount = getIntegerParameter(PULLING_PULLED_COUNT_STRING, 0, 0);
	this->pulled = (int*)malloc(this->pulledCount*sizeof(int));
	printf("%d resid(s) fixed, %d pulled.\n", this->fixedCount, this->pulledCount);
	char paramName[10];
	for(i = 0; i < this->fixedCount; i++){
		sprintf(paramName, "%s%d", PULLING_FIXED_STRING, i+1);
		this->fixed[i] = getIntegerParameter(paramName, 0, 0);
		printf("Resid %d is fixed.\n", this->fixed[i]);
	}
	for(i = 0; i < this->pulledCount; i++){
		sprintf(paramName, "%s%d", PULLING_PULLED_STRING, i+1);
		this->pulled[i] = getIntegerParameter(paramName, 0, 0);
		printf("Pulling resid %d.\n", this->pulled[i]);
	}

	char pullDirection[30];
	getMaskedParameter(pullDirection, PULLING_DIRECTION_STRING, DEFAULT_PULLING_DIRECTION, 1);
	float3 pullVector;
	if(strcmp(pullDirection, PULLING_DIRECTION_VECTOR_STRING) == 0){
		this->fixedEnd = getIntegerParameter(PULLING_FIXED_END_STRING, 1, 1);
		this->pulledEnd = getIntegerParameter(PULLING_PULLED_END_STRING, 2, 1);
		getVectorParameter(PULLING_VECTOR_STRING, &pullVector.x, &pullVector.y, &pullVector.z, 0, 0, 0, 0);
		for(traj = 0; traj < gsop.Ntr; traj++){
			this->pullVector[traj] = pullVector;
			this->cantCoord0[traj] = make_float3(gsop.h_coord[traj*sop.aminoCount + this->pulledEnd]);
			this->cantCoord[traj]  = make_float3(gsop.h_coord[traj*sop.aminoCount + this->pulledEnd]);
		}
		printf("Pulling in direction of vector (%f, %f, %f).\n", pullVector.x, pullVector.y, pullVector.z);
	} else if(strcmp(pullDirection, PULLING_DIRECTION_ENDTOEND_STRING) == 0){
		this->fixedEnd = getIntegerParameter(PULLING_FIXED_END_STRING, 0, 0);
		this->pulledEnd = getIntegerParameter(PULLING_PULLED_END_STRING, 0, 0);
		for(traj = 0; traj < gsop.Ntr; traj++){
			this->pullVector[traj] = make_float3(
                    gsop.h_coord[traj*sop.aminoCount + this->pulledEnd]
					- gsop.h_coord[traj*sop.aminoCount + this->fixedEnd]);
			normalize(this->pullVector[traj]);
			this->cantCoord0[traj] = make_float3(gsop.h_coord[traj*sop.aminoCount + this->pulledEnd]);
			this->cantCoord[traj]  = make_float3(gsop.h_coord[traj*sop.aminoCount + this->pulledEnd]);
		}
		printf("Pulling in end-to-end direction: %c%d(%s) - %c%d(%s).\n",
				sop.aminos[this->fixedEnd].chain, sop.aminos[this->fixedEnd].resid, sop.aminos[this->fixedEnd].resName,
				sop.aminos[this->pulledEnd].chain, sop.aminos[this->pulledEnd].resid, sop.aminos[this->pulledEnd].resName);
		printf("Pulling vector is (%5.2f, %5.2f, %5.2f).\n",
				this->pullVector[0].x,
				this->pullVector[0].y,
				this->pullVector[0].z);
	} else {
		DIE("ERROR: 'pullDirection' parameter should be set to 'endToEnd' or 'vector'.\n");
	}

	for(traj = 0; traj < gsop.Ntr; traj++){
		for(i = 0; i < sop.aminoCount; i++){
			this->h_extForces[traj*sop.aminoCount + i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
		}
		for(j = 0; j < this->pulledCount; j++){
			i = this->pulled[j];
			if(traj == 0){
				printf("Pulling bead #%d (%s %d chain %c).\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
			}
			this->h_extForces[traj*sop.aminoCount + i] = make_float4(
					this->pullVector[traj] * this->fconst, 2.0f);
		}
		for(j = 0; j < this->fixedCount; j++){
			i = this->fixed[j];
			if(traj == 0){
				printf("Fixing bead #%d (%s %d chain %c).\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
			}
			this->h_extForces[traj*sop.aminoCount + i] = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
		}
	}

	for(i = 0; i < sop.aminoCount; i++){
		sop.aminos[i].beta = this->h_extForces[i].w;
	}
	cudaMemcpy(this->d_extForces, this->h_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
    checkCUDAError();

    this->updateParametersOnGPU();

	if(this->deltax != 0.0f){
		pullFilenames.resize(gsop.Ntr);
        std::string pullFilename;
        pullFilename = getMaskedParameterAs<std::string>(PULLING_FILENAME, DEFAULT_PULLING_FILENAME);
		for(traj = 0; traj < gsop.Ntr; traj++){
			pullFilenames[traj] = string_replace(pullFilename, "<run>", traj+gsop.firstrun);
			FILE* pullFile = safe_fopen(pullFilenames[traj].c_str(), "w");
			fclose(pullFile);
		}
		printf("Pulling data will be saved in '%s'.\n", pullFilename.c_str());
	}
	
    this->blockSize = gsop.blockSize;
	this->blockNum = gsop.aminoCount/this->blockSize + 1;

	printf("Done initializing pulling protocol...\n");
}

void PullingPotential::compute(){
	// The force is added in integrator
}

void PullingPotential::updateForces(float xt){
	copyCoordDeviceToHost();
	int traj, j;
	for(traj = 0; traj < gsop.Ntr; traj++){
		this->extForce[traj] = this->computeForce(gsop.h_coord[sop.aminoCount*traj + this->pulledEnd], traj);
		// Increasing the force'
		this->cantCoord[traj] = this->cantCoord0[traj] + xt * this->pullVector[traj];
		for(j = 0; j < this->pulledCount; j++){
			this->h_extForces[traj*sop.aminoCount + this->pulled[j]] =
					make_float4(this->extForce[traj], 2.0f);
		}
	}
	cudaMemcpy(this->d_extForces, this->h_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
}

void PullingPotential::savePullingData(){
    int traj;
	if(gsop.step % 100000 == 0){
		printf("'Cantilever chip' coordinates for run #%d: %f, %f, %f\n",
				gsop.firstrun, this->cantCoord[0].x, this->cantCoord[0].y, this->cantCoord[0].z);
		printf("'Cantilever tip' coordinates for run #%d: %f, %f, %f\n",
				gsop.firstrun, gsop.h_coord[this->pulledEnd].x, gsop.h_coord[this->pulledEnd].y, gsop.h_coord[this->pulledEnd].z);
	}
	for(traj = 0; traj < gsop.Ntr; traj++){
		FILE* pullFile = safe_fopen(pullFilenames[traj].c_str(), "a");

        float3 endToEnd_vector = make_float3(gsop.h_coord[traj*sop.aminoCount + this->pulledEnd] - gsop.h_coord[traj*sop.aminoCount + this->fixedEnd]);
		float endToEnd_x = dot( endToEnd_vector , this->pullVector[traj] );
		float endToEnd = abs( endToEnd_vector );
        float f = dot(this->extForce[traj] , this->pullVector[traj]);

		fprintf(pullFile, "%12ld\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",
				gsop.step, endToEnd, endToEnd_x, f,
				this->extForce[traj].x, this->extForce[traj].y, this->extForce[traj].z);

		fclose(pullFile);
	}
	checkCUDAError();
}

void PullingPotential::updateParametersOnGPU(){
    hc_pulling.d_extForces = this->d_extForces;
	cudaMemcpyToSymbol(c_pulling, &hc_pulling, sizeof(PullingConstant), 0, cudaMemcpyHostToDevice);
	checkCUDAError();
}

void PullingUpdater::update(){
    float xt = pulling->deltax*(gsop.step / this->frequency);
    pulling->updateForces(xt);
	if(gsop.step % this->frequency == 0){
        pulling->savePullingData();
	}
}

float3 PullingPotential::computeForce(const float4 &coordN, int traj) const{
	float3 f;

	f.x = this->Ks * (this->cantCoord[traj].x - coordN.x);
	f.y = this->Ks * (this->cantCoord[traj].y - coordN.y);
	f.z = this->Ks * (this->cantCoord[traj].z - coordN.z);

	return f;
}

