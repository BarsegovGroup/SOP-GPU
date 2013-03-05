/*
 * pulling.cu
 *
 *  Created on: May 25, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "pulling.cuh"
#include "pulling_kernel.cu"

char** pullFilenames;
FILE* pullFile;

float3 computeForce(float4 coordN, int traj);

extern void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace);

void createPullingPotential(){
	if(pullingOn == 1 || getYesNoParameter(PULLING_ON_STRING, 0, 1) == 1){
		pullingOn = 1;
		gsop.pullingOn = 1;
		//sprintf(pullingPotential.name, "Pulling");
		//pullingPotential.compute = &computePulling;
		//pullingPotential.computeEnergy = &computePullingEnergy;
		//potentials[potentialsCount] = &pullingPotential;
		//if(gsop.deviceProp.major == 2){
		//	cudaFuncSetCacheConfig(pulling_kernel, cudaFuncCachePreferL1);
		//}
		//potentialsCount++;
		if(getFloatParameter(PULLING_DELTAX_STRING, 0, 1) == 0 && getFloatParameter(PULLING_FCONST_STRING, 0, 1) == 0){
			printf("ERROR: Either 'deltax' or 'fconst' parameter should be specified to initiate pulling.\n");
			exit(0);
		}
		if(getFloatParameter(PULLING_DELTAX_STRING, 0, 1) != 0){
			sprintf(pullingUpdater.name, "Pulling");
			pullingUpdater.update = &updatePulling;
			pullingUpdater.frequency = getIntegerParameter(PULLING_FREQ, nav, 1);
			updaters[updatersCount] = &pullingUpdater;
			updatersCount++;
		}
		initPulling();
	}
}

void initPulling(){
	printf("Initializing pulling protocol...\n");
	int i, j, traj;
	pulling.pullVector = (float3*)calloc(Ntr, sizeof(float3));
	pulling.cantCoord0 = (float3*)calloc(Ntr, sizeof(float3));
	pulling.cantCoord = (float3*)calloc(Ntr, sizeof(float3));
	pulling.extForce = (float3*)calloc(Ntr, sizeof(float3));
	pulling.h_extForces = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMalloc((void**)&pulling.d_extForces, gsop.aminoCount*sizeof(float4));

	pulling.deltax = getFloatParameter(PULLING_DELTAX_STRING, 0, 1);
	pulling.Ks = getFloatParameter(PULLING_KS_STRING, DEFAULT_PULLING_KS, 1);
	pulling.fconst = getFloatParameter(PULLING_FCONST_STRING, 0, 1);

	if(pulling.deltax == 0 && pulling.fconst == 0){
		printf("ERROR: Either 'deltax' or 'fconst' parameter should be specified to initiate pulling.\n");
		exit(0);
	}

	pulling.fixedCount = getIntegerParameter(PULLING_FIXED_COUNT_STRING, 0, 0);
	pulling.fixed = (int*)malloc(pulling.fixedCount*sizeof(int));
	pulling.pulledCount = getIntegerParameter(PULLING_PULLED_COUNT_STRING, 0, 0);
	pulling.pulled = (int*)malloc(pulling.pulledCount*sizeof(int));
	printf("%d resid(s) fixed, %d pulled.\n", pulling.fixedCount, pulling.pulledCount);
	char paramName[10];
	for(i = 0; i < pulling.fixedCount; i++){
		sprintf(paramName, "%s%d\0", PULLING_FIXED_STRING, i+1);
		pulling.fixed[i] = getIntegerParameter(paramName, 0, 0);
		/*if(fixed_beads[i] < 0 || fixed_beads[i] >= sop.aminoCount){
			printf("ERROR: Fixed bead %d not exists. Protein has only %d amino-acids. Bead numbers should start with zero.\n", fixed_beads[i], sop.aminoCount);
			exit(0);
		}*/
		printf("Resid %d is fixed.\n", pulling.fixed[i]);
	}
	for(i = 0; i < pulling.pulledCount; i++){
		sprintf(paramName, "%s%d\0", PULLING_PULLED_STRING, i+1);
		pulling.pulled[i] = getIntegerParameter(paramName, 0, 0);
		/*if(pulled_beads[i] < 0 || pulled_beads[i] >= sop.aminoCount){
			printf("ERROR: Pulled bead %d not exists. Protein has only %d amino-acids. Bead numbers should start with zero.\n", pulled_beads[i], sop.aminoCount);
			exit(0);
		}*/
		printf("Pulling resid %d.\n", pulling.pulled[i]);
	}

	char pullDirection[30];
	getMaskedParameter(pullDirection, PULLING_DIRECTION_STRING, DEFAULT_PULLING_DIRECTION, 1);
	float3 pullVector;
	if(strcmp(pullDirection, PULLING_DIRECTION_VECTOR_STRING) == 0){
		pulling.fixedEnd = getIntegerParameter(PULLING_FIXED_END_STRING, 1, 1);
		pulling.pulledEnd = getIntegerParameter(PULLING_PULLED_END_STRING, 2, 1);
		getVectorParameter(PULLING_VECTOR_STRING, &pullVector.x, &pullVector.y, &pullVector.z, 0, 0, 0, 0);
		for(traj = 0; traj < Ntr; traj++){
			pulling.pullVector[traj].x = pullVector.x;
			pulling.pullVector[traj].y = pullVector.y;
			pulling.pullVector[traj].z = pullVector.z;
			pulling.cantCoord0[traj].x = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].x;
			pulling.cantCoord0[traj].y = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].y;
			pulling.cantCoord0[traj].z = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].z;
			pulling.cantCoord[traj].x = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].x;
			pulling.cantCoord[traj].y = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].y;
			pulling.cantCoord[traj].z = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].z;
		}
		printf("Pulling in direction of vector (%f, %f, %f).\n", pullVector.x, pullVector.y, pullVector.z);
	} else if(strcmp(pullDirection, PULLING_DIRECTION_ENDTOEND_STRING) == 0){
		pulling.fixedEnd = getIntegerParameter(PULLING_FIXED_END_STRING, 0, 0);
		pulling.pulledEnd = getIntegerParameter(PULLING_PULLED_END_STRING, 0, 0);
		for(traj = 0; traj < Ntr; traj++){
			pulling.pullVector[traj].x = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].x
					- gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].x;
			pulling.pullVector[traj].y = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].y
					- gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].y;
			pulling.pullVector[traj].z = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].z
					- gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].z;
			float norm = sqrtf(pulling.pullVector[traj].x*pulling.pullVector[traj].x +
					pulling.pullVector[traj].y*pulling.pullVector[traj].y +
					pulling.pullVector[traj].z*pulling.pullVector[traj].z);
			pulling.pullVector[traj].x /= norm;
			pulling.pullVector[traj].y /= norm;
			pulling.pullVector[traj].z /= norm;
			pulling.cantCoord0[traj].x = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].x;
			pulling.cantCoord0[traj].y = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].y;
			pulling.cantCoord0[traj].z = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].z;
			pulling.cantCoord[traj].x = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].x;
			pulling.cantCoord[traj].y = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].y;
			pulling.cantCoord[traj].z = gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].z;
		}
		printf("Pulling in end-to-end direction: %c%d(%s) - %c%d(%s).\n",
				sop.aminos[pulling.fixedEnd].chain, sop.aminos[pulling.fixedEnd].resid, sop.aminos[pulling.fixedEnd].resName,
				sop.aminos[pulling.pulledEnd].chain, sop.aminos[pulling.pulledEnd].resid, sop.aminos[pulling.pulledEnd].resName);
		printf("Pulling vector is (%5.2f, %5.2f, %5.2f).\n",
				pulling.pullVector[0].x,
				pulling.pullVector[0].y,
				pulling.pullVector[0].z);
	} else {
		printf("ERROR: 'pullDirection' parameter should be set to 'endToEnd' or 'vector'.\n", pullVector.x, pullVector.y, pullVector.z);
		exit(0);
	}

	for(traj = 0; traj < Ntr; traj++){
		for(i = 0; i < sop.aminoCount; i++){
			pulling.h_extForces[traj*sop.aminoCount + i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
		}
		for(j = 0; j < pulling.pulledCount; j++){
			i = pulling.pulled[j];
			if(traj == 0){
				printf("Pulling bead #%d (%s %d chain %c).\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
			}
			pulling.h_extForces[traj*sop.aminoCount + i] = make_float4(
					pulling.pullVector[traj].x*pulling.fconst,
					pulling.pullVector[traj].y*pulling.fconst,
					pulling.pullVector[traj].z*pulling.fconst, 2.0);
		}
		for(j = 0; j < pulling.fixedCount; j++){
			i = pulling.fixed[j];
			if(traj == 0){
				printf("Fixing bead #%d (%s %d chain %c).\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
			}
			pulling.h_extForces[traj*sop.aminoCount + i] = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
		}
		/*for(i = 0; i < sop.aminoCount; i++){
			printf("%d - %f\n", i, pulling.h_extForces[traj*sop.aminoCount + i].w);
		}*/
	}

	for(i = 0; i < sop.aminoCount; i++){
		sop.aminos[i].beta = pulling.h_extForces[i].w;
	}
	cudaMemcpy(pulling.d_extForces, pulling.h_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_pulling, &pulling, sizeof(Pulling), 0, cudaMemcpyHostToDevice);
	checkCUDAError();

	if(pulling.deltax != 0.0f){
		pullFilenames = (char**)calloc(Ntr, sizeof(char*));
		char tempFilename[100];

		getMaskedParameter(tempFilename, PULLING_FILENAME, DEFAULT_PULLING_FILENAME, 1);;
		for(traj = 0; traj < Ntr; traj++){
			pullFilenames[traj] = (char*)calloc(100, sizeof(char));
			char trajnum[10];
			sprintf(trajnum, "%d\0", traj+firstrun);
			replaceString(pullFilenames[traj], tempFilename, trajnum, "<run>");
			pullFile = fopen(pullFilenames[traj], "w");
			fclose(pullFile);
		}
		printf("Pulling data will be saved in '%s'.\n", tempFilename);
	}

	printf("Done initializing pulling protocol...\n");
}

inline void computePulling(){
	pulling_kernel<<<covalent.blockNum, covalent.blockSize>>>();
	checkCUDAError();
}

inline void computePullingEnergy(){

}

inline void updatePulling(){

	copyCoordDeviceToHost();
	int traj, j;
	for(traj = 0; traj < Ntr; traj++){
		pulling.extForce[traj] = computeForce(gsop.h_coord[sop.aminoCount*traj + pulling.pulledEnd], traj);
		// Increasing the force
		xt = pulling.deltax*(step / pullingUpdater.frequency);
		pulling.cantCoord[traj].x = pulling.cantCoord0[traj].x + xt * pulling.pullVector[traj].x;
		pulling.cantCoord[traj].y = pulling.cantCoord0[traj].y + xt * pulling.pullVector[traj].y;
		pulling.cantCoord[traj].z = pulling.cantCoord0[traj].z + xt * pulling.pullVector[traj].z;
		for(j = 0; j < pulling.pulledCount; j++){
			pulling.h_extForces[traj*sop.aminoCount + pulling.pulled[j]] =
					make_float4(pulling.extForce[traj].x, pulling.extForce[traj].y, pulling.extForce[traj].z, 2.0);
		}
	}
	cudaMemcpy(pulling.d_extForces, pulling.h_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
	if(step % pullingUpdater.frequency == 0){
		if(step % 100000 == 0){
			printf("'Cantilever chip' coordinates for run #%d: %f, %f, %f\n",
					firstrun, pulling.cantCoord[0].x, pulling.cantCoord[0].y, pulling.cantCoord[0].z);
			printf("'Cantilever tip' coordinates for run #%d: %f, %f, %f\n",
					firstrun, gsop.h_coord[pulling.pulledEnd].x, gsop.h_coord[pulling.pulledEnd].y, gsop.h_coord[pulling.pulledEnd].z);
		}
		for(traj = 0; traj < Ntr; traj++){
			pullFile = fopen(pullFilenames[traj], "a");

			float endToEnd_x = (gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].x - gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].x)*pulling.pullVector[traj].x +
					(gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].y - gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].y)*pulling.pullVector[traj].y +
					(gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].z - gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].z)*pulling.pullVector[traj].z;
			float endToEnd = sqrtf((gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].x-gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].x)*
					(gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].x-gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].x) +
					(gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].y-gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].y)*
					(gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].y-gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].y) +
					(gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].z-gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].z)*
					(gsop.h_coord[traj*sop.aminoCount + pulling.pulledEnd].z-gsop.h_coord[traj*sop.aminoCount + pulling.fixedEnd].z));
			float f = pulling.extForce[traj].x*pulling.pullVector[traj].x + pulling.extForce[traj].y*pulling.pullVector[traj].y + pulling.extForce[traj].z*pulling.pullVector[traj].z;

			fprintf(pullFile, "%12ld\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",
					step, endToEnd, endToEnd_x, f,
					pulling.extForce[traj].x, pulling.extForce[traj].y, pulling.extForce[traj].z);

			fflush(pullFile);
			fclose(pullFile);
		}



		checkCUDAError();
	}
}

float3 computeForce(float4 coordN, int traj){

	float3 f = make_float3(0.0f, 0.0f, 0.0f);

	f.x = pulling.Ks * (pulling.cantCoord[traj].x - coordN.x);
	f.y = pulling.Ks * (pulling.cantCoord[traj].y - coordN.y);
	f.z = pulling.Ks * (pulling.cantCoord[traj].z - coordN.z);

	//printf("Pulling force: (%f, %f, %f)\n", f.x, f.y, f.z);

	return f;

}
