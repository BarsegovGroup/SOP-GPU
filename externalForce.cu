/*
 * externalForce.cu
 *
 *  Created on: Jun 20, 2009
 *      Author: zhmurov
 */
#include "gsop.cuh"
#include "externalForce.cuh"


void engageCantileverTip();
void doPulling();
float3 computeForce(float4 coordN, int traj);

void initExternalForce(){
	int i, j;
	ext_force = (float3*)calloc(Ntr, sizeof(float3));
	pull_vector = (float3*)calloc(Ntr, sizeof(float3));
	cantilever_coord = (float3*)calloc(Ntr, sizeof(float3));
	cantilever_coord0 = (float3*)calloc(Ntr, sizeof(float3));
	gsop.h_extForces = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMalloc((void**)&gsop.d_extForces, gsop.aminoCount*sizeof(float4));
	if(stage == pull_stage){
#ifdef capsid
		float xcm = 0, ycm = 0, zcm = 0;
		for(i = 0; i < sop.aminoCount; i++){
			xcm += gsop.h_coord[i].x;
			ycm += gsop.h_coord[i].y;
			zcm += gsop.h_coord[i].z;
		}
		xcm /= sop.aminoCount;
		ycm /= sop.aminoCount;
		zcm /= sop.aminoCount;
		printf("Center of mass: (%f, %f, %f)\n", xcm, ycm, zcm);
		float rx, ry, rz;
		pulled_count = 0;
		for(i = 0; i < sop.aminoCount; i++){
			rx = gsop.h_coord[i].x - xcm;
			ry = gsop.h_coord[i].y - xcm;
			rz = gsop.h_coord[i].z - xcm;
			if(sqrtf(rx*rx + ry*ry + rz*rz) < 230.0f){
				//printf("%f\n", sqrtf(rx*rx + ry*ry + rz*rz));
				c_gsop.d_extForces[i].w = 2;
				pulled_count ++;
			} else {
				c_gsop.d_extForces[i].w = 0;
			}
		}
#else
		int traj;
		for(traj = 0; traj < Ntr; traj++){
			for(i = 0; i < sop.aminoCount; i++){
				for(j = 0; j < pulled_count; j++){
					if(pulled_beads[j] == i){
						if(traj == 0){
							printf("Pulling bead #%d (%s %d chain %c):\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
						}
						gsop.h_extForces[traj*sop.aminoCount + i] = make_float4(ext_force[traj].x, ext_force[traj].y, ext_force[traj].z, 2.0);
					} else {
						gsop.h_extForces[traj*sop.aminoCount + i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
					}
				}
				for(j = 0; j < fixed_count; j++){
					if(fixed_beads[j] == i){
						if(traj == 0){
							printf("Fixing bead #%d (%s %d chain %c):\n", i, sop.aminos[i].resName, sop.aminos[i].resid, sop.aminos[i].chain);
						}
						gsop.h_extForces[traj*sop.aminoCount + i] = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
					}
				}
				sop.aminos[i].beta = (int)gsop.h_extForces[i].w;
			}
		}
#endif
	} else {
		for(i = 0; i < sop.aminoCount; i++){
			sop.aminos[i].beta = 0.0f;
		}
		sop.aminos[fixed_end].beta = 1.0f;
		sop.aminos[pulled_end].beta = 2.0f;
	}
	gsop.stage = stage;
	cudaMemcpy(gsop.d_extForces, gsop.h_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
	//engageCantileverTip();
}

void engageCantileverTip(){
	int i;
	if(restart == 0){
		if(strcmp(pullDirection, "vector") == 0){
			for(i = 0; i < Ntr; i++){
				pull_vector[i].x = pullVectorX;
				pull_vector[i].y = pullVectorY;
				pull_vector[i].z = pullVectorZ;
			}
		} else if(strcmp(pullDirection, "endToEnd") == 0){
			for(i = 0; i < Ntr; i++){
				pull_vector[i].x = gsop.h_coord[i*sop.aminoCount + pulled_end].x - gsop.h_coord[i*sop.aminoCount + fixed_end].x;
				pull_vector[i].y = gsop.h_coord[i*sop.aminoCount + pulled_end].y - gsop.h_coord[i*sop.aminoCount + fixed_end].y;
				pull_vector[i].z = gsop.h_coord[i*sop.aminoCount + pulled_end].z - gsop.h_coord[i*sop.aminoCount + fixed_end].z;
			}
		}
		if(stage != pull_stage){
			for(i = 0; i < Ntr; i++){
				pull_vector[i].x = gsop.h_coord[i*sop.aminoCount + pulled_end].x - gsop.h_coord[i*sop.aminoCount + fixed_end].x;
				pull_vector[i].y = gsop.h_coord[i*sop.aminoCount + pulled_end].y - gsop.h_coord[i*sop.aminoCount + fixed_end].y;
				pull_vector[i].z = gsop.h_coord[i*sop.aminoCount + pulled_end].z - gsop.h_coord[i*sop.aminoCount + fixed_end].z;
			}
		}
		for(i = 0; i < Ntr; i++){
			float length = sqrtf(pull_vector[i].x*pull_vector[i].x + pull_vector[i].y*pull_vector[i].y + pull_vector[i].z*pull_vector[i].z);
			if(length != 0){
				pull_vector[i].x /= length;
				pull_vector[i].y /= length;
				pull_vector[i].z /= length;
			}
			cantilever_coord0[i].x = gsop.h_coord[i*sop.aminoCount + pulled_end].x;
			cantilever_coord0[i].y = gsop.h_coord[i*sop.aminoCount + pulled_end].y;
			cantilever_coord0[i].z = gsop.h_coord[i*sop.aminoCount + pulled_end].z;
			cantilever_coord[i].x = gsop.h_coord[i*sop.aminoCount + pulled_end].x;
			cantilever_coord[i].y = gsop.h_coord[i*sop.aminoCount + pulled_end].y;
			cantilever_coord[i].z = gsop.h_coord[i*sop.aminoCount + pulled_end].z;
		}
	} else {
		for(i = 0; i < Ntr; i++){
			pull_vector[i].x = pullVectorX;
			pull_vector[i].y = pullVectorY;
			pull_vector[i].z = pullVectorZ;
			cantilever_coord0[i].x = chipCoordX - pullVectorX*xt;
			cantilever_coord0[i].y = chipCoordY - pullVectorY*xt;
			cantilever_coord0[i].z = chipCoordZ - pullVectorZ*xt;
			cantilever_coord[i].x = chipCoordX;
			cantilever_coord[i].y = chipCoordY;
			cantilever_coord[i].z = chipCoordZ;
		}
	}

	printf("Initial 'cantilever chip' coordinates: (%f, %f, %f)\n", cantilever_coord0[0].x, cantilever_coord0[0].y, cantilever_coord0[0].z);
	printf("Pulling vector: (%f, %f, %f)\n", pull_vector[0].x, pull_vector[0].y, pull_vector[0].z);
}

/*
 * Pulling Nth aminoacid
 */
void doPulling(){

#ifdef capsid
	//ext_force.x = (k_trans * step) / ((float)pulled_count);
	//cudaMemcpyToSymbol(c_ext_force, &ext_force, sizeof(float3), 0, cudaMemcpyHostToDevice);
#else
	int i, j;
	copyCoordDeviceToHost();
	for(i = 0; i < Ntr; i++){
		ext_force[i] = computeForce(gsop.h_coord[sop.aminoCount*i + pulled_end], i);

		// Increasing the force
		xt = deltax*(step / nav);
		cantilever_coord[i].x = cantilever_coord0[i].x + xt * pull_vector[i].x;
		cantilever_coord[i].y = cantilever_coord0[i].y + xt * pull_vector[i].y;
		cantilever_coord[i].z = cantilever_coord0[i].z + xt * pull_vector[i].z;
		for(j = 0; j < pulled_count; j++){
			gsop.h_extForces[i*sop.aminoCount + pulled_beads[j]] = make_float4(ext_force[i].x, ext_force[i].y, ext_force[i].z, 2.0);
		}
	}
	cudaMemcpy(gsop.d_extForces, gsop.h_extForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
	if(step % 100000 == 0){
		printf("'Cantilever chip' coordinates: %f, %f, %f\n", cantilever_coord[0].x, cantilever_coord[0].y, cantilever_coord[0].z);
		printf("'Cantilever tip' coordinates: %f, %f, %f\n", gsop.h_coord[pulled_end].x, gsop.h_coord[pulled_end].y, gsop.h_coord[pulled_end].z);
	}

#endif
}

/*
 * Computing force vector based on coordinates of first and last amino acids
 * (Copied from initial (CPU) code)
 */
float3 computeForce(float4 coordN, int traj){

	float3 f = make_float3(0.0f, 0.0f, 0.0f);

	f.x = k_trans * (cantilever_coord[traj].x - coordN.x);
	f.y = k_trans * (cantilever_coord[traj].y - coordN.y);
	f.z = k_trans * (cantilever_coord[traj].z - coordN.z);

	return f;

}
