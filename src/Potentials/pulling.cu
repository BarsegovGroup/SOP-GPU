/*
 * pulling.cu
 *
 *  Created on: May 25, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "../Util/mystl.h"
#include "../Util/wrapper.h"
#include "../Util/vector_helpers.h"
#include "pulling.h"

#include "pulling_kernel.cu"

void createPullingPotential(){
	if(parameters::pulling.get()){
		gsop.pullingOn = 1;

        PullingPotential *pot;
		potentials[potentialsCount] = pot = new PullingPotential();
        potentialsCount++;

		if(parameters::deltax.get() != 0){
			updaters[updatersCount] = new PullingUpdater(pot);
			updatersCount++;
		}
	}
}

PullingUpdater::PullingUpdater(PullingPotential *pulling){
	this->name = "Pulling";
    this->frequency = parameters::pullFreq.get();
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

	this->deltax = parameters::deltax.get();
	this->Ks = parameters::k_trans.get();
	this->fconst = parameters::fconst.get();

	if(this->deltax == 0 && this->fconst == 0){
		DIE("ERROR: Either 'deltax' or 'fconst' parameter should be specified to initiate pulling\n");
	}

    if (parameters::_is_defined("fixed_beads") || parameters::_is_defined("pulled_beads")) {
        DIE("'fixed_beads' and 'pulled_beads' parameters are deprecated. Use 'fixed' and 'pulled' instead");
    }
	this->fixed = parameters::fixed.get();
	this->pulled = parameters::pulled.get();
	printf("%ld resid(s) fixed, %ld pulled.\n", this->fixed.size(), this->pulled.size());

    std::string pullDirection = parameters::pullDirection.get();
	float3 pullVector;
	if(pullDirection == "vector"){
		this->fixedEnd = parameters::fixedEnd.get();
		this->pulledEnd = parameters::pulledEnd.get();
        pullVector = parameters::pullVector.get();
		for(traj = 0; traj < gsop.Ntr; traj++){
			this->pullVector[traj] = pullVector;
			this->cantCoord0[traj] = make_float3(gsop.h_coord[traj*sop.aminoCount + this->pulledEnd]);
			this->cantCoord[traj]  = make_float3(gsop.h_coord[traj*sop.aminoCount + this->pulledEnd]);
		}
		printf("Pulling in direction of vector (%f, %f, %f).\n", pullVector.x, pullVector.y, pullVector.z);
	} else if(pullDirection == "endToEnd"){
		this->fixedEnd = parameters::fixedEnd.get();
		this->pulledEnd = parameters::pulledEnd.get();
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
		for(j = 0; j < this->pulled.size(); j++){
			i = this->pulled[j];
			if (i > sop.aminoCount || i < 0) {
				DIE("Invalid pulling bead #%d; should be >=0 and <%d", i, sop.aminoCount);
			}
			if(traj == 0){
				printf("Pulling bead #%d (%s %d chain %c).\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
			}
			this->h_extForces[traj*sop.aminoCount + i] = make_float4(
					this->pullVector[traj] * this->fconst, 2.0f);
		}
		for(j = 0; j < this->fixed.size(); j++){
			i = this->fixed[j];
			if (i > sop.aminoCount || i < 0) {
				DIE("Invalid fixed bead #%d; should be >=0 and <%d", i, sop.aminoCount);
			}
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
        std::string pullFilename = parameters::pullOutput.get();
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
		for(j = 0; j < this->pulled.size(); j++){
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

