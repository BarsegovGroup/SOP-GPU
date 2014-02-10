/*
 * pulling_plane.cu
 *
 *  Created on: Jan 17, 2014
 *      Author: alekseenko
 */
#include "../gsop.cuh"
#include "../Integrators/langevin.cuh"
#include "pulling_plane.cuh"
#include "pulling_plane_kernel.cu"

char pullingPlaneFilename[500];
FILE* pullingPlaneFile;

extern void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace);

void createPullingPlanePotential(){
	if(getYesNoParameter(PULLINGPLANE_ON_STRING, 0, 1) == 1){
		gsop.pullingPlaneOn = 1;

		sprintf(pullingPlanePotential.name, "Pulling Plane");
		pullingPlanePotential.compute = &computePullingPlane;
		pullingPlanePotential.computeEnergy = &computePullingPlaneEnergy;
		potentials[potentialsCount] = &pullingPlanePotential;
		if(gsop.deviceProp.major == 2){
			cudaFuncSetCacheConfig(pullingPlane_kernel, cudaFuncCachePreferL1);
		}
		potentialsCount++;

		sprintf(pullingPlaneUpdater.name, "Pulling Plane");
		pullingPlaneUpdater.update = &updatePullingPlane;
		pullingPlaneUpdater.frequency = getIntegerParameter(PULLINGPLANE_FREQ, nav, 1);
		updaters[updatersCount] = &pullingPlaneUpdater;
		updatersCount++;

		initPullingPlane();
	}
}

void initPullingPlane(){
	printf("Initializing pulling plane protocol...\n");
    if (Ntr != 1) {
        printf("Pulling plane can only run in single-trajectory-per-GPU mode (runnum 1)\n");
        exit(-1);
    }
	int i, j;
	pullingPlane.h_extForces = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMalloc((void**)&pullingPlane.d_extForces, gsop.aminoCount*sizeof(float4));

	pullingPlane.deltax = getFloatParameter(PULLINGPLANE_DELTAX_STRING, 0, 1);
	pullingPlane.Ks = getFloatParameter(PULLINGPLANE_KS_STRING, DEFAULT_PULLINGPLANE_KS, 1);
	pullingPlane.mass = getFloatParameter(PULLINGPLANE_MASS_STRING, 0, 0); // in "mass of single AA residue" units

	if(pullingPlane.deltax == 0){
		printf("ERROR: '%s' parameter should be specified to initiate pullingPlane.\n", PULLINGPLANE_DELTAX_STRING);
		exit(-1);
	}

	pullingPlane.fixedCount = getIntegerParameter(PULLINGPLANE_FIXED_COUNT_STRING, 0, 0);
	pullingPlane.fixed = (int*)malloc(pullingPlane.fixedCount*sizeof(int));
	pullingPlane.pulledCount = getIntegerParameter(PULLINGPLANE_PULLED_COUNT_STRING, 0, 0);
	pullingPlane.pulled = (int*)malloc(pullingPlane.pulledCount*sizeof(int));
	printf("%d resid(s) fixed, %d pulled.\n", pullingPlane.fixedCount, pullingPlane.pulledCount);
	char paramName[50];
	for(i = 0; i < pullingPlane.fixedCount; i++){
		sprintf(paramName, "%s%d\0", PULLINGPLANE_FIXED_STRING, i+1);
		pullingPlane.fixed[i] = getIntegerParameter(paramName, 0, 0);
		if(pullingPlane.fixed[i] < 0 || pullingPlane.fixed[i] >= gsop.aminoCount){
			printf("ERROR: Fixed bead %s %d not exists. Protein has only %d amino-acids. Bead numbers should start with zero.\n", paramName, pullingPlane.fixed[i], gsop.aminoCount);
			exit(-1);
		}
		printf("Resid %d is fixed.\n", pullingPlane.fixed[i]);
	}
	for(i = 0; i < pullingPlane.pulledCount; i++){
		sprintf(paramName, "%s%d\0", PULLINGPLANE_PULLED_STRING, i+1);
		pullingPlane.pulled[i] = getIntegerParameter(paramName, 0, 0);
		if(pullingPlane.pulled[i] < 0 || pullingPlane.pulled[i] >= gsop.aminoCount){
			printf("ERROR: Pulled bead %s %d not exists. Protein has only %d amino-acids. Bead numbers should start with zero.\n", paramName, pullingPlane.pulled[i], gsop.aminoCount);
			exit(-1);
		}
		printf("Pulling resid %d.\n", pullingPlane.pulled[i]);
	}

    getVectorParameter(PULLINGPLANE_PULLVECTOR_STRING, &pullingPlane.pullVector.x, &pullingPlane.pullVector.y, &pullingPlane.pullVector.z, 0, 0, 0, 0);
    double t = sqrt(pullingPlane.pullVector.x*pullingPlane.pullVector.x + pullingPlane.pullVector.y*pullingPlane.pullVector.y + pullingPlane.pullVector.z*pullingPlane.pullVector.z);
    pullingPlane.pullVector.x /= t;
    pullingPlane.pullVector.y /= t;
    pullingPlane.pullVector.z /= t;
    getVectorParameter(PULLINGPLANE_ZEROVECTOR_STRING, &pullingPlane.planeCoord.x, &pullingPlane.planeCoord.y, &pullingPlane.planeCoord.z, 0, 0, 0, 0);
    pullingPlane.d = - (pullingPlane.planeCoord.x*pullingPlane.pullVector.x + pullingPlane.planeCoord.y*pullingPlane.pullVector.y + pullingPlane.planeCoord.z*pullingPlane.pullVector.z);
    pullingPlane.cant_d = pullingPlane.d;

	for(i = 0; i < sop.aminoCount; i++){
		pullingPlane.h_extForces[i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
	}
	for(j = 0; j < pullingPlane.pulledCount; j++){
		i = pullingPlane.pulled[j];
		printf("Pulling bead #%d (%s %d chain %c).\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
		pullingPlane.h_extForces[i] = make_float4(0.0, 0.0, 0.0, 2.0);
	}
	for(j = 0; j < pullingPlane.fixedCount; j++){
		i = pullingPlane.fixed[j];
		printf("Fixing bead #%d (%s %d chain %c).\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
		pullingPlane.h_extForces[i] = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
	}

	for(i = 0; i < sop.aminoCount; i++){
		sop.aminos[i].beta = pullingPlane.h_extForces[i].w;
	}
	cudaMemcpy(pullingPlane.d_extForces, pullingPlane.h_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_pullingPlane, &pullingPlane, sizeof(PullingPlane), 0, cudaMemcpyHostToDevice);
	checkCUDAError();

    char tmpstr[512];
	getMaskedParameter(tmpstr, PULLINGPLANE_FILENAME, "", 0);
	char trajnum[10];
	sprintf(trajnum, "%d\0", firstrun);
	replaceString(pullingPlaneFilename, tmpstr, trajnum, "<run>");
	pullingPlaneFile = fopen(pullingPlaneFilename, "w");
	fclose(pullingPlaneFile);
	printf("PullingPlane data will be saved in '%s'.\n", pullingPlaneFilename);

	printf("Done initializing pulling plane protocol...\n");
}

inline void computePullingPlane(){
	pullingPlane_kernel<<<covalent.blockNum, covalent.blockSize>>>();
	checkCUDAError();
}

inline void computePullingPlaneEnergy(){
    // Perl has operator "..."
    // It should be used inside functions to be implemented
    // Unlike just empty functions, the warning is generated when "..." is used
    // This comment is completely useless and barely related, because this function is empty by design
}

inline void updatePullingPlane(){

	//copyCoordDeviceToHost();
	int j;
	// Increasing the force
	xt = pullingPlane.deltax * step / pullingPlaneUpdater.frequency;
    /*
	pullingPlane.cantCoord.x = pullingPlane.planeCoord0.x + xt * pullingPlane.pullVector.x;
	pullingPlane.cantCoord.y = pullingPlane.planeCoord0.y + xt * pullingPlane.pullVector.y;
	pullingPlane.cantCoord.z = pullingPlane.planeCoord0.z + xt * pullingPlane.pullVector.z;
    */
    pullingPlane.cant_d = pullingPlane.d0 + xt;
	checkCUDAError();
	cudaMemcpy(pullingPlane.h_extForces, pullingPlane.d_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyDeviceToHost);
	checkCUDAError();
	if(step % pullingPlaneUpdater.frequency == 0){
		if(step % 100000 == 0){
			printf("Cantilever coordinates for run #%d: %f, %f, %f\n",
					firstrun, pullingPlane.cantCoord.x, pullingPlane.cantCoord.y, pullingPlane.cantCoord.z);
			printf("Plane coordinates for run #%d: %f, %f, %f\n",
					firstrun, pullingPlane.planeCoord.x, pullingPlane.planeCoord.y, pullingPlane.planeCoord.z);
		}
        pullingPlaneFile = fopen(pullingPlaneFilename, "a");
        float3 extForce = make_float3(0.f, 0.f, 0.f);
        for (j = 0; j < gsop.aminoCount; j++){
            extForce.x += pullingPlane.h_extForces[j].x;
            extForce.y += pullingPlane.h_extForces[j].y;
            extForce.z += pullingPlane.h_extForces[j].z;
        }

        float extForceProj = extForce.x*pullingPlane.pullVector.x + extForce.y*pullingPlane.pullVector.y + extForce.z*pullingPlane.pullVector.z;

        float totForce = pullingPlane.Ks * (pullingPlane.cant_d - pullingPlane.d) - extForceProj;

        pullingPlane.d += totForce * pullingPlaneUpdater.frequency * langevin.hOverZeta / pullingPlane.mass;
	
        cudaMemcpyToSymbol(c_pullingPlane, &pullingPlane, sizeof(PullingPlane), 0, cudaMemcpyHostToDevice);

        fprintf(pullingPlaneFile, "%12ld\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",
                step, pullingPlane.d, pullingPlane.cant_d, extForceProj,
                extForce.x, extForce.y, extForce.z);

        fflush(pullingPlaneFile);
        fclose(pullingPlaneFile);
    }

	checkCUDAError();
}

