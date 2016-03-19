/*
 * pulling_plane.cu
 *
 *  Created on: Jan 17, 2014
 *      Author: alekseenko
 */
#include "../gsop.cuh"
#include "../Integrators/langevin.h"
#include "pulling_plane.h"

FILE* pullingPlaneFile;
std::vector<std::string> pullingPlaneFilename;

#include "pulling_plane_kernel.cu"

void createPullingPlanePotential(){
	if(parameters::pullingPlane.get()){
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
	this->frequency = parameters::pullingPlaneUpdateFreq.get();
	this->outputFreq = parameters::pullingPlaneOutputFreq.get();
	this->pullingPlane = pullingPlane;
}

PullingPlanePotential::PullingPlanePotential(){
	this->name = "Pulling Plane";
	printf("Initializing pulling plane protocol...\n");

	// Parameters for std out
	this->outputWidth = parameters::outputcolwidth.get();
	this->printRuns = parameters::printruns.get();

	/*if (gsop.Ntr != 1) {
        	DIE("Pulling plane can only run in single-trajectory-per-GPU mode (runnum 1)\n");
	}*/
	int i, j;
	this->h_extForces = (float*)calloc(gsop.aminoCount, sizeof(float));
	cudaMalloc((void**)&this->d_extForces, gsop.aminoCount*sizeof(float));
	this->h_masks = (int*)calloc(gsop.aminoCount, sizeof(int));
	cudaMalloc((void**)&this->d_masks, gsop.aminoCount*sizeof(int));
	this->h_r0 = (float*)calloc(gsop.aminoCount, sizeof(float));
	cudaMalloc((void**)&this->d_r0, gsop.aminoCount*sizeof(float));
	this->h_planeDispl = (float*)calloc(gsop.Ntr, sizeof(float));
	cudaMalloc((void**)&this->d_planeDispl, gsop.Ntr*sizeof(float));

	this->deltax = parameters::pullingPlaneDeltax.get();
	this->Ks = parameters::pullingPlaneKs.get();
	this->KsCant = parameters::pullingPlaneKsCant.get();
	this->zeta = parameters::pullingPlaneZeta.get();

	if(this->deltax == 0){
		DIE("'%s' parameter should be specified to initiate this->\n", parameters::pullingPlaneDeltax.name().c_str());
	}

	if (parameters::_is_defined("plane_fixed_beads") || parameters::_is_defined("plane_pulled_beads")) {
		DIE("'plane_fixed_beads' and 'plane_pulled_beads' parameters are deprecated. Use 'plane_fixed' and 'plane_pulled' instead");
	}
	this->fixed = parameters::plane_fixed.get();
	this->pulled = parameters::plane_pulled.get();
	printf("%ld resid(s) fixed, %ld pulled.\n", this->fixed.size(), this->pulled.size());

	this->pullVector = parameters::pullingPlaneDir.get();
	double t = abs(this->pullVector);
	this->pullVector /= t;
	this->planeCoord0 = parameters::pullingPlanePos.get();
	this->cantDispl = 0;

	for(i = 0; i < gsop.aminoCount; i++){
		float3 dr;
		dr.x = gsop.h_coord[i].x - this->planeCoord0.x;
		dr.y = gsop.h_coord[i].y - this->planeCoord0.y;
		dr.z = gsop.h_coord[i].z - this->planeCoord0.z;
		this->h_r0[i] = dr.x*this->pullVector.x + dr.y*this->pullVector.y + dr.z*this->pullVector.z;
		this->h_extForces[i] = 0.0f;
	}
	int traj;
	for(traj = 0; traj < gsop.Ntr; traj++){
		for(j = 0; j < this->pulled.size(); j++){
			i = this->pulled[j];
			printf("Pulling bead #%lu (%s %d chain %c).\n", i + traj*sop.aminos.size(), sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
			this->h_masks[i + traj*sop.aminos.size()] = 2;
		}
		for(j = 0; j < this->fixed.size(); j++){
			i = this->fixed[j];
			printf("Fixing bead #%lu (%s %d chain %c).\n", i + traj*sop.aminos.size(), sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
			this->h_masks[i + traj*sop.aminos.size()] = 1;
		}
	}

	for(i = 0; i < sop.aminos.size(); i++){
		sop.aminos[i].beta = (float)this->h_masks[i];
	}
	cudaMemcpy(this->d_extForces, this->h_extForces, gsop.aminoCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_masks, this->h_masks, gsop.aminoCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_r0, this->h_r0, gsop.aminoCount*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(this->d_planeDispl, this->h_planeDispl, gsop.Ntr*sizeof(float), cudaMemcpyHostToDevice);
	checkCUDAError();

	hc_pullingPlane.d_extForces = this->d_extForces;
	hc_pullingPlane.d_masks = this->d_masks;
	hc_pullingPlane.d_r0 = this->d_r0;
	hc_pullingPlane.pullVector = this->pullVector;
	hc_pullingPlane.planeCoord0 = this->planeCoord0;
	hc_pullingPlane.d_planeDispl = this->d_planeDispl;
	hc_pullingPlane.Ks = this->Ks;
	cudaMemcpyToSymbol(c_pullingPlane, &hc_pullingPlane, sizeof(PullingPlaneConstant), 0, cudaMemcpyHostToDevice);
	checkCUDAError();

	for(traj = 0; traj < gsop.Ntr; traj++){
		pullingPlaneFilename.push_back(parameters::pullingPlaneOutput.replace("<run>", gsop.firstrun+traj));
		pullingPlaneFile = safe_fopen(pullingPlaneFilename.at(traj).c_str(), "w");
		fclose(pullingPlaneFile);
		printf("PullingPlane data will be saved in '%s'.\n", pullingPlaneFilename.at(traj).c_str());
	}
	
	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(pullingPlane_kernel, cudaFuncCachePreferL1);
	}

	this->blockSize = gsop.blockSize;
	this->blockNum = gsop.aminoCount/this->blockSize + 1;

	printf("Done initializing pulling plane protocol...\n");
}

void PullingPlanePotential::compute(){
	pullingPlane_kernel<<<this->blockNum, this->blockSize>>>();
	checkCUDAError();
}

/*int PullingPlanePotential::getEnergiesCount(){
	return 0;
}

float* PullingPlanePotential::computeEnergy(int id){
	DIE("No energy terms are returned by pulling plane potential.");
	return NULL;
    // Perl has operator "..."
    // It should be used inside functions to be implemented
    // Unlike just empty functions, the warning is generated when "..." is used
    // This comment is completely useless and barely related, because this function is empty by design
}*/

void PullingPlaneUpdater::update(){

	//copyCoordDeviceToHost();
	int i;
	// Increasing the force
	float xt = pullingPlane->deltax * gsop.step / this->frequency;
	pullingPlane->cantDispl = xt;

	cudaMemcpy(pullingPlane->h_extForces, pullingPlane->d_extForces, gsop.aminoCount*sizeof(float), cudaMemcpyDeviceToHost);
	if(gsop.step % this->outputFreq == 0){
		printf("%*s%*s%*s%*s%*s%*s%*s\n",
					pullingPlane->outputWidth, "Run", 
					pullingPlane->outputWidth, "Timestep", 
					pullingPlane->outputWidth, "Plane displ.", 
					pullingPlane->outputWidth, "Cant. displ.", 
					pullingPlane->outputWidth, "F cant.", 
					pullingPlane->outputWidth, "F atomic", 
					pullingPlane->outputWidth, "End-to-end");
	}
	int traj;
	for(traj = 0; traj < gsop.Ntr; traj++){
		float extForce = 0.0f;
		for(i = traj*sop.aminos.size(); i < (traj+1)*sop.aminos.size(); i++){
			if(pullingPlane->h_masks[i] == 2){
				extForce += pullingPlane->h_extForces[i];
			}
		}
		//extForce = 0.0;
		float cantForce = (xt - pullingPlane->h_planeDispl[traj])*pullingPlane->KsCant;
		pullingPlane->h_planeDispl[traj] += (this->frequency*integrator->h/pullingPlane->zeta)*(cantForce - extForce);

		checkCUDAError();
		checkCUDAError();
		/*if(gsop.step % this->outputFreq == 0){
			float3 cantCoord;
			cantCoord.x = pullingPlane->planeCoord0.x + xt * pullingPlane->pullVector.x;
			cantCoord.y = pullingPlane->planeCoord0.y + xt * pullingPlane->pullVector.y;
			cantCoord.z = pullingPlane->planeCoord0.z + xt * pullingPlane->pullVector.z;
			printf("Cantilever coordinates for run #%d: %f, %f, %f (displacement: %f)\n",
					gsop.firstrun+traj, cantCoord.x, cantCoord.y, cantCoord.z, xt);
			float3 planeCoords;
			planeCoords.x = pullingPlane->planeCoord0.x + pullingPlane->h_planeDispl[traj] * pullingPlane->pullVector.x;
			planeCoords.y = pullingPlane->planeCoord0.y + pullingPlane->h_planeDispl[traj] * pullingPlane->pullVector.y;
			planeCoords.z = pullingPlane->planeCoord0.z + pullingPlane->h_planeDispl[traj] * pullingPlane->pullVector.z;
			printf("Plane coordinates for run #%d: %f, %f, %f (displacement: %f)\n",
					gsop.firstrun+traj, planeCoords.x, planeCoords.y, planeCoords.z, pullingPlane->h_planeDispl[traj]);
		}*/
	
		if(gsop.step % this->outputFreq == 0){
			copyCoordDeviceToHost();
			float3 cmfix = make_float3(0.0, 0.0, 0.0);
			float3 cmpull = make_float3(0.0, 0.0, 0.0);
			int fixCount = 0;
			int pullCount = 0;
			for(i = traj*sop.aminos.size(); i < (traj+1)*sop.aminos.size(); i++){
			   	if(pullingPlane->h_masks[i] == 1){
					cmfix.x += gsop.h_coord[i].x;
					cmfix.y += gsop.h_coord[i].y;
					cmfix.z += gsop.h_coord[i].z;
					fixCount ++;
				}
				if(pullingPlane->h_masks[i] == 2){
					cmpull.x += gsop.h_coord[i].x;
					cmpull.y += gsop.h_coord[i].y;
					cmpull.z += gsop.h_coord[i].z;
					pullCount ++;
				}
			}
			cmfix.x /= fixCount;
			cmfix.y /= fixCount;
			cmfix.z /= fixCount;
			cmpull.x /= pullCount;
			cmpull.y /= pullCount;
			cmpull.z /= pullCount;
			float dx = cmpull.x - cmfix.x;
			float dy = cmpull.y - cmfix.y;
			float dz = cmpull.z - cmfix.z;
			float endToEnd = sqrtf(dx*dx + dy*dy + dz*dz);
		
			pullingPlaneFile = safe_fopen(pullingPlaneFilename.at(traj).c_str(), "a");
			fprintf(pullingPlaneFile, "%12ld\t%f\t%f\t%f\t%f\t%f\n",
				gsop.step, pullingPlane->h_planeDispl[traj], pullingPlane->cantDispl, cantForce, extForce, endToEnd);
			if(traj < pullingPlane->printRuns){
				printf("%*d%*ld%*f%*f%*f%*f%*f\n",
					pullingPlane->outputWidth, gsop.firstrun+traj, 
					pullingPlane->outputWidth, gsop.step, 
					pullingPlane->outputWidth, pullingPlane->h_planeDispl[traj], 
					pullingPlane->outputWidth, pullingPlane->cantDispl, 
					pullingPlane->outputWidth, cantForce, 
					pullingPlane->outputWidth, extForce, 
					pullingPlane->outputWidth, endToEnd);
			}
			if(traj == pullingPlane->printRuns){
				printf("%*s%*s%*s%*s%*s%*s%*s\n",
					pullingPlane->outputWidth, "...", 
					pullingPlane->outputWidth, "...", 
					pullingPlane->outputWidth, "...", 
					pullingPlane->outputWidth, "...", 
					pullingPlane->outputWidth, "...", 
					pullingPlane->outputWidth, "...", 
					pullingPlane->outputWidth, "...");
			}

			fclose(pullingPlaneFile);
		}
	}

	cudaMemcpy(pullingPlane->d_planeDispl, pullingPlane->h_planeDispl, gsop.Ntr*sizeof(float), cudaMemcpyHostToDevice);
       //cudaMemcpyToSymbol(c_pullingPlane, &hc_pullingPlane, sizeof(PullingPlaneConstant), 0, cudaMemcpyHostToDevice);

	checkCUDAError();
}

