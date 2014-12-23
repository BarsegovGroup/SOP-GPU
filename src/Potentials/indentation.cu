/*
 * indentation.cu
 *
 *  Created on: Apr 9, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "../Util/wrapper.h"
#include "indentation.h"
#include "../Updaters/output_manager.h"

int showTipMica;

#include "indentation_kernel.cu"

/*
 * Add potential computation and update functions to lists,
 * create timer.
 */

void createIndentationPotential(){
	if(gsop.indentationOn == 1 || getYesNoParameter(INDENTATION_ON_STRING, 0, 1) == 1){
		if(gsop.Ntr > 1){
			DIE("Many-runs-per-GPU mode is not supported by indentation potential.\n");
		}
		gsop.indentationOn = 1;

        IndentationPotential *pot;
		potentials[potentialsCount] = pot = new IndentationPotential();
		potentialsCount++;

		updaters[updatersCount] = new IndentationTipUpdater(pot);
		updatersCount++;
        if (sop.additionalAminos > 0){
            // Reorder indentationTipUpdater and indentationAminoUpdater
            // FIXME: this is ugly
            std::swap(updaters[updatersCount-1], updaters[updatersCount-2]);
        }
	}
}

/*
 * Initialize potential:
 * Tip and surface are both represented by a repulsive LJ-potential.
 * Tip is a sphere, connected to cantilever base by a spring.
 * Motion of a tip is described by Langevin equation but not integrated every step,
 * since tip is much heavier then particles.
 */
IndentationPotential::IndentationPotential(){
	this->name = "Indentation";
	this->discreteSurf = getYesNoParameter(INDENTATION_DISCRETE_SURF_STRING, 0, 1);

	printf("Initializing indentation protocol...\n");
	// Read poteintial from configuration file
	getVectorParameter(INDENTATION_CHIP_POSITION_STRING, &hc_indentation.chipCoord.x, &hc_indentation.chipCoord.y, &hc_indentation.chipCoord.z, 0.0, 0.0, 0.0, 0);
	hc_indentation.chipCoord0 = hc_indentation.chipCoord;
	getVectorParameter(INDENTATION_TIP_POSITION_STRING, &hc_indentation.tipCoord.x, &hc_indentation.tipCoord.y, &hc_indentation.tipCoord.z,
			hc_indentation.chipCoord.x, hc_indentation.chipCoord.y, hc_indentation.chipCoord.z, 1);
	getVectorParameter(INDENTATION_DIRECTION_STRING, &hc_indentation.direction.x, &hc_indentation.direction.y, &hc_indentation.direction.z, 0.0, 0.0, 0.0, 0);
	getVectorParameter(INDENTATION_SURFACE_R0_STRING, &hc_indentation.micaR0.x, &hc_indentation.micaR0.y, &hc_indentation.micaR0.z, 0.0, 0.0, 0.0, 0);
	hc_indentation.micaR = hc_indentation.micaR0;
	getVectorParameter(INDENTATION_SURFACE_N_STRING, &hc_indentation.micaN.x, &hc_indentation.micaN.y, &hc_indentation.micaN.z, 0.0, 0.0, 0.0, 0);
	hc_indentation.tipRadius = getFloatParameter(INDENTATION_TIP_RADIUS_STRING, 0.0, 0);
	hc_indentation.cantileverKs = getFloatParameter(INDENTATION_KS_STRING, 0.0, 0);
	hc_indentation.V = getFloatParameter(INDENTATION_DELTAX_STRING, 0.0, 0);
	hc_indentation.dx = 0;
	hc_indentation.moveSurface = getYesNoParameter(INDENTATION_MOVE_SURFACE, 0, 1);
	hc_indentation.fixTransversal = getYesNoParameter(INDENTATION_FIX_TRANSVERSAL,
			DEFAULT_INDENTATION_FIX_TRANSVERSAL, 1);
	float tempVar;
	tempVar = getFloatParameter(INDENTATION_SIGMA_STRING, 1.0, 1);
	hc_indentation.tipa6 = powf(getFloatParameter(INDENTATION_TIP_SIGMA_STRING, tempVar, 1), 6.0);
	hc_indentation.surfa6 = powf(getFloatParameter(INDENTATION_SURF_SIGMA_STRING, tempVar, 1), 6.0);
	tempVar = getFloatParameter(INDENTATION_EL_STRING, 1.0, 1);
	hc_indentation.tipel = getFloatParameter(INDENTATION_TIP_EL_STRING, tempVar, 1);
	hc_indentation.surfel = getFloatParameter(INDENTATION_SURF_EL_STRING, tempVar, 1);

	float tipA, tipB, surfA, surfB;
	tipA = getFloatParameter(INDENTATION_TIP_A, DEFAULT_INDENTATION_TIP_A, 1);
	tipB = getFloatParameter(INDENTATION_TIP_B, DEFAULT_INDENTATION_TIP_B, 1);
	surfA = getFloatParameter(INDENTATION_SURFACE_A, DEFAULT_INDENTATION_SURFACE_A, 1);
	surfB = getFloatParameter(INDENTATION_SURFACE_B, DEFAULT_INDENTATION_SURFACE_B, 1);
	hc_indentation.tipAprime = 12.0f*hc_indentation.tipel*hc_indentation.tipa6*hc_indentation.tipa6*tipA;
	hc_indentation.tipBprime = 6.0f*hc_indentation.tipel*hc_indentation.tipa6*tipB;
	hc_indentation.surfAprime = 12.0f*hc_indentation.surfel*hc_indentation.surfa6*hc_indentation.surfa6*surfA;
	hc_indentation.surfBprime = 6.0f*hc_indentation.surfel*hc_indentation.surfa6*surfB;

	/*float x = 201;
	FILE* tempFile = safe_fopen("temp.dat", "w");
	while(x < 300){
		float xprime = pow(x - hc_indentation.tipRadius, 6.0);
		float y = hc_indentation.el*((tipA*hc_indentation.a6*hc_indentation.a6)/(xprime*xprime) + (tipB*hc_indentation.a6)/(xprime));
		fprintf(tempFile, "%f\t%f\n", x, y);
		x += 0.1;
	}
	fclose(tempFile);
	exit(-1);*/

	hc_indentation.tipZeta = getFloatParameter(INDENTATION_TIP_ZETA_STRING, 5000.0, 1);
	char INDENTATION_filename[100];
	getMaskedParameter(INDENTATION_filename, INDENTATION_OUTPUT_FILENAME_STRING, DEFAULT_INDENTATION_FILENAME, 1);

	printf("Initial cantilever base coordinates: (%5.2f, %5.2f, %5.2f).\n", hc_indentation.chipCoord.x, hc_indentation.chipCoord.y, hc_indentation.chipCoord.z);

	// Normilize vectors
	float norm = sqrt(hc_indentation.direction.x*hc_indentation.direction.x + hc_indentation.direction.y*hc_indentation.direction.y + hc_indentation.direction.z*hc_indentation.direction.z);
	hc_indentation.direction.x /= norm;
	hc_indentation.direction.y /= norm;
	hc_indentation.direction.z /= norm;
	norm = sqrt(hc_indentation.micaN.x*hc_indentation.micaN.x + hc_indentation.micaN.y*hc_indentation.micaN.y + hc_indentation.micaN.z*hc_indentation.micaN.z);
	hc_indentation.micaN.x /= norm;
	hc_indentation.micaN.y /= norm;
	hc_indentation.micaN.z /= norm;

	// Allocatae memory, copy data to the device
	hc_indentation.h_tipForces = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMalloc((void**)&hc_indentation.d_tipForces, gsop.aminoCount*sizeof(float4));
	cudaMemcpy(hc_indentation.d_tipForces, hc_indentation.h_tipForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
	hc_indentation.out = safe_fopen(INDENTATION_filename, "w");
	hc_indentation.retractionStep = getLongIntegerParameter(INDENTATION_RETRACTION_STEP_STRING, -1, 1);
	checkCUDAError();

	// Add dummy particles to show surface and cantilever in dcd
	showTipMica = getYesNoParameter(INDENTATION_SHOW_TIP_SURFACE_STRING, 0, 1);
	if(!this->discreteSurf){
		printf("Discrete mica representation requires saving mica beads into dcd. Forcing \"showTipSurf\".\n");
		showTipMica = 1;
	}

	if(showTipMica){
		updaters[updatersCount] = new IndentationAminoUpdater();
		updatersCount++;
	} else {
		sop.additionalAminosCount = 0;
	}
	if(this->discreteSurf){

		hc_indentation.pairsCutoff2 = getFloatParameter(INDENTATION_PAIRS_CUTOFF_STRING,
				DEFAULT_INDENTATION_PAIRS_CUTOFF, 1);
		hc_indentation.pairsCutoff2 = hc_indentation.pairsCutoff2*hc_indentation.pairsCutoff2;
		hc_indentation.surfaceBeadsCount = hc_indentation.surfaceSize*hc_indentation.surfaceSize;

		hc_indentation.h_surfacePointsCoord = (float4*)calloc(hc_indentation.surfaceSize*hc_indentation.surfaceSize,
					sizeof(float4));
		int i, j;
		for(i = 0; i < hc_indentation.surfaceSize; i++){
			for(j = 0; j < hc_indentation.surfaceSize; j++){
				hc_indentation.h_surfacePointsCoord[i*hc_indentation.surfaceSize + j].x =
						hc_indentation.surfacePointsR0[i*hc_indentation.surfaceSize + j].x;
				hc_indentation.h_surfacePointsCoord[i*hc_indentation.surfaceSize + j].y =
						hc_indentation.surfacePointsR0[i*hc_indentation.surfaceSize + j].y;
				hc_indentation.h_surfacePointsCoord[i*hc_indentation.surfaceSize + j].z =
						hc_indentation.surfacePointsR0[i*hc_indentation.surfaceSize + j].z;
			}
		}


		printf("Computing maximum pairs in mica pairslist..\n");
		cudaMalloc((void**)&hc_indentation.d_surfacePointsCoord, hc_indentation.surfaceSize*hc_indentation.surfaceSize*sizeof(float4));

		int maxPairs = 0;
		cudaMallocHost((void**)&hc_indentation.h_micaListCounts, gsop.aminoCount*sizeof(int));
		for(i = 0; i < gsop.aminoCount; i++){
			hc_indentation.h_micaListCounts[i] = 0;
		}
		float4 r1, r2;
		for(i = 0; i < gsop.aminoCount; i++){
			for(j = 0; j < hc_indentation.surfaceBeadsCount; j++){
				r1 = gsop.h_coord[i];
				r2 = hc_indentation.h_surfacePointsCoord[j];
				r1.x -= r2.x;
				r1.y -= r2.y;
				r1.z -= r2.z;
				r1.w = r1.x*r1.x + r1.y*r1.y + r1.z*r1.z;
				if(r1.w < hc_indentation.pairsCutoff2){
					hc_indentation.h_micaListCounts[i] ++;
				}
			}
		}
		for(i = 0; i < gsop.aminoCount; i++){
			if(hc_indentation.h_micaListCounts[i] > maxPairs){
				maxPairs = hc_indentation.h_micaListCounts[i];
			}
		}
		printf("Done. Maximum pairs found is %d.\n", maxPairs);
		maxPairs += 256;
		if(maxPairs > hc_indentation.surfaceBeadsCount){
			maxPairs = hc_indentation.surfaceBeadsCount;
		}


		cudaMallocHost((void**)&hc_indentation.h_micaListCounts, gsop.aminoCount*sizeof(int));
		cudaMallocHost((void**)&hc_indentation.h_micaList,	maxPairs*gsop.width*sizeof(int));
		cudaMalloc((void**)&hc_indentation.d_micaListCounts, gsop.aminoCount*sizeof(int));
		cudaMalloc((void**)&hc_indentation.d_micaList, maxPairs*gsop.width*sizeof(int));

		cudaMemcpy(hc_indentation.d_surfacePointsCoord, hc_indentation.h_surfacePointsCoord,
				hc_indentation.surfaceSize*hc_indentation.surfaceSize*sizeof(float4), cudaMemcpyHostToDevice);
	}
	hc_indentation.outputFreq = getIntegerParameter(INDENTATION_OUTPUT_FREQ_STRING, DEFAULT_INDENTATION_OUTPUT_FREQ, 1);
	cudaMemcpyToSymbol(c_indentation, &hc_indentation, sizeof(hc_indentation), 0, cudaMemcpyHostToDevice);

	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(indentation_kernel, cudaFuncCachePreferL1);
	}
}

/*
 * Create dummy particles, showing grid as a surface and a couple of beads for cantilever
 */
IndentationAminoUpdater::IndentationAminoUpdater(){
	this->name = "AdditionalAminos";
	this->frequency = getIntegerParameter(OUTPUT_FREQUENCY_STRING, DEFAULT_OUTPUT_FREQUENCY, 1);

	hc_indentation.surfaceSize = getIntegerParameter(INDENTATION_MICA_SIZE_STRING, 51, 1);
	hc_indentation.surfaceStep = getFloatParameter(INDENTATION_MICA_STEP_STRING, 1.4, 1);
	sop.additionalAminosCount = 2 + hc_indentation.surfaceSize*hc_indentation.surfaceSize;
	sop.additionalAminos = (PDBAtom*)calloc(sop.additionalAminosCount, sizeof(PDBAtom));
	sop.additionalAminos[0].id = sop.aminoCount;
	strcpy(sop.additionalAminos[0].resName, "TIP");
	sop.additionalAminos[0].chain = 'T';
	sop.additionalAminos[0].x = hc_indentation.chipCoord.x;
	sop.additionalAminos[0].y = hc_indentation.chipCoord.y;
	sop.additionalAminos[0].z = hc_indentation.chipCoord.z;
	sop.additionalAminos[0].resid = 1;
	sop.additionalAminos[1].id = sop.aminoCount + 1;
	strcpy(sop.additionalAminos[1].resName, "TIP");
	sop.additionalAminos[1].chain = 'T';
	sop.additionalAminos[1].x = hc_indentation.tipCoord.x;
	sop.additionalAminos[1].y = hc_indentation.tipCoord.y;
	sop.additionalAminos[1].z = hc_indentation.tipCoord.z;
	sop.additionalAminos[1].resid = 2;
	float3 r1, r2, n, r0;
	float3 r;
	float len;
	r0 = hc_indentation.direction;
	if(r0.x != 0.0f){
		r.y = 1.0f;
		r.z = 0.0f;
		r.x = -r.y*r0.y/r0.x;
		len = sqrtf(r.x*r.x + r.y*r.y + r.z*r.z);
		r.x /= len;
		r.y /= len;
		r.z /= len;
	} else {
		r.x = 1.0f;
		r.y = 0.0f;
		r.z = 0.0f;
	}
	len = getFloatParameter(INDENTATION_CANTILEVER_LENGTH, DEFAULT_INDENTATION_CANTILEVER_LENGTH, 1);
	r.x *= len;
	r.y *= len;
	r.z *= len;
	printf("Cantilever representation vector: (%5.3f, %5.3f, %5.3f)\n", r.x, r.y, r.z);
	hc_indentation.cantileverVector = r;
	sop.additionalAminos[0].x += r.x;
	sop.additionalAminos[0].y += r.y;
	sop.additionalAminos[0].z += r.z;

	hc_indentation.surfacePointsR0 = (float3*)calloc(hc_indentation.surfaceSize*hc_indentation.surfaceSize,
			sizeof(float3));
	n = hc_indentation.micaN;
	r0 = hc_indentation.micaR0;
	if(n.x == 0){
		r1.x = 1.0f;
		r1.y = 0;
		r1.z = 0;
	} else if(n.y == 0){
		r1.x = 0;
		r1.y = 1.0f;
		r1.z = 0;
	} else if(n.z == 0){
		r1.x = 0;
		r1.y = 0;
		r1.z = 1.0f;
	} else {
		r1.x = 0;
		r1.y = 1.0f;
		r1.z = - (n.y/n.z)*r1.y;
	}
	float norm = sqrt(r1.x*r1.x + r1.y*r1.y + r1.z*r1.z);
	r1.x /= norm;
	r1.y /= norm;
	r1.z /= norm;
	r2.x = n.y*r1.z - n.z*r1.y;
	r2.y = n.z*r1.x - n.x*r1.z;
	r2.z = n.x*r1.y - n.y*r1.x;
	int i, j;
	float displ = hc_indentation.surfaceStep;
	for(i = 0; i < hc_indentation.surfaceSize; i++){
		for(j = 0; j < hc_indentation.surfaceSize; j++){
			int i1 = i-hc_indentation.surfaceSize/2;
			int i2 = j-hc_indentation.surfaceSize/2;
			r.x = r0.x + i1*displ*r1.x + i2*displ*r2.x;
			r.y = r0.y + i1*displ*r1.y + i2*displ*r2.y;
			r.z = r0.z + i1*displ*r1.z + i2*displ*r2.z;
			sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].id =
					sop.aminoCount + 2 + i*hc_indentation.surfaceSize+j;
			strcpy(sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].resName, "MIC");
			sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].chain = 'M';
			sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].x = r.x;
			sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].y = r.y;
			sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].z = r.z;
			hc_indentation.surfacePointsR0[i*hc_indentation.surfaceSize + j] = r;
		}
	}
	hc_indentation.chipCoordAv.x = 0.0f;
	hc_indentation.chipCoordAv.y = 0.0f;
	hc_indentation.chipCoordAv.z = 0.0f;
	hc_indentation.tipCoordAv.x = 0.0f;
	hc_indentation.tipCoordAv.y = 0.0f;
	hc_indentation.tipCoordAv.z = 0.0f;
	hc_indentation.fav.x = 0.0f;
	hc_indentation.fav.y = 0.0f;
	hc_indentation.fav.z = 0.0f;
	hc_indentation.fav.w = 0.0f;
	char tempFilename[100];
	getMaskedParameter(tempFilename, INDENTATION_VMD_CONNECT_SCRIPT, "connect_mica.vmd", 1);
	FILE* micaConnectFile = safe_fopen(tempFilename, "w");
	printf("Dumping VMD script to connect surface elements into '%s'.\n", tempFilename);
	fprintf(micaConnectFile, "set sel [atomselect top \"index %d\"]\n", gsop.aminoCount);
	fprintf(micaConnectFile, "set bonds {%d}\n", gsop.aminoCount + 1);
	fprintf(micaConnectFile, "$sel setbonds [list $bonds]\n$sel delete\n");
	fprintf(micaConnectFile, "set sel [atomselect top \"index %d\"]\n", gsop.aminoCount + 1);
	fprintf(micaConnectFile, "set bonds {%d}\n", gsop.aminoCount);
	fprintf(micaConnectFile, "$sel setbonds [list $bonds]\n$sel delete\n");
	for(i = 0; i < hc_indentation.surfaceSize; i++){
		for(j = 0; j < hc_indentation.surfaceSize; j++){
			int i1 = 2+i*hc_indentation.surfaceSize+j + gsop.aminoCount;
			int i2;
			fprintf(micaConnectFile, "set sel [atomselect top \"index %d\"]\n", i1);
			fprintf(micaConnectFile, "set bonds {");
			if(i != hc_indentation.surfaceSize - 1){
				i2 = 2+(i+1)*hc_indentation.surfaceSize+j + gsop.aminoCount;
				fprintf(micaConnectFile, "%d", i2);
			}
			if(j != hc_indentation.surfaceSize - 1){
				i2 = 2+i*hc_indentation.surfaceSize+j+1 + gsop.aminoCount;
				fprintf(micaConnectFile, " %d", i2);
			}
			fprintf(micaConnectFile, "}\n");
			fprintf(micaConnectFile, "$sel setbonds [list $bonds]\n$sel delete\n");

		}
	}
	fclose(micaConnectFile);
}

/*
 * Execute computations on GPU
 */
void IndentationPotential::compute(){
	if(this->discreteSurf){
	    indentationDiscreteSurf_kernel<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
	} else {
	    indentation_kernel<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
	}
}

IndentationTipUpdater::IndentationTipUpdater(IndentationPotential *indentation) {
	this->name = "Indentation";
    this->frequency = nav;
    this->indentation = indentation;
}

/*
 * Move cantilever
 */
void IndentationTipUpdater::update(){
	// Change moving direction to switch to retraction
	if(step == hc_indentation.retractionStep){
		printf("Starting tip retraction at step %ld.\n", hc_indentation.retractionStep);
		/*hc_indentation.direction.x = -hc_indentation.direction.x;
		hc_indentation.direction.y = -hc_indentation.direction.y;
		hc_indentation.direction.z = -hc_indentation.direction.z;*/
		hc_indentation.V = -hc_indentation.V;
	}
	if(step % this->frequency == 0){
		// Copy forces acting onto tip to host memory
		cudaMemcpy(hc_indentation.h_tipForces, hc_indentation.d_tipForces, gsop.aminoCount*sizeof(float4), cudaMemcpyDeviceToHost);
		int i;
		float4 f = make_float4(0.0, 0.0, 0.0, 0.0);
		// Sum up all force and set 0 value to forces
		for(i = 0; i < gsop.aminoCount; i++){
			f.x += hc_indentation.h_tipForces[i].x;
			f.y += hc_indentation.h_tipForces[i].y;
			f.z += hc_indentation.h_tipForces[i].z;
			hc_indentation.h_tipForces[i].x = 0.0f;
			hc_indentation.h_tipForces[i].y = 0.0f;
			hc_indentation.h_tipForces[i].z = 0.0f;
		}
		cudaMemcpy(hc_indentation.d_tipForces, hc_indentation.h_tipForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
		// Compute average force and its absolute value
		f.x /= (float)this->frequency;
		f.y /= (float)this->frequency;
		f.z /= (float)this->frequency;
		f.w = sqrtf(f.x*f.x + f.y*f.y + f.z*f.z);
		hc_indentation.tipForce = f;
		// Move chip or surface
		hc_indentation.dx += hc_indentation.V;//(step/nav)*hc_indentation.chipV;
		if(hc_indentation.moveSurface == 1){
			if(this->indentation->discreteSurf){
				for(i = 0; i < hc_indentation.surfaceBeadsCount; i++){
					hc_indentation.h_surfacePointsCoord[i].x = hc_indentation.surfacePointsR0[i].x + hc_indentation.direction.x*hc_indentation.dx;
					hc_indentation.h_surfacePointsCoord[i].y = hc_indentation.surfacePointsR0[i].y + hc_indentation.direction.y*hc_indentation.dx;
					hc_indentation.h_surfacePointsCoord[i].z = hc_indentation.surfacePointsR0[i].z + hc_indentation.direction.z*hc_indentation.dx;
				}
				cudaMemcpy(hc_indentation.d_surfacePointsCoord, hc_indentation.h_surfacePointsCoord,
						hc_indentation.surfaceBeadsCount*sizeof(float4), cudaMemcpyHostToDevice);
				generateMicaList_kernel<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
			} else {
				hc_indentation.micaR.x = hc_indentation.micaR0.x + hc_indentation.direction.x*hc_indentation.dx;
				hc_indentation.micaR.y = hc_indentation.micaR0.y + hc_indentation.direction.y*hc_indentation.dx;
				hc_indentation.micaR.z = hc_indentation.micaR0.z + hc_indentation.direction.z*hc_indentation.dx;
			}
		} else {
			hc_indentation.chipCoord.x = hc_indentation.chipCoord0.x + hc_indentation.direction.x*hc_indentation.dx;
			hc_indentation.chipCoord.y = hc_indentation.chipCoord0.y + hc_indentation.direction.y*hc_indentation.dx;
			hc_indentation.chipCoord.z = hc_indentation.chipCoord0.z + hc_indentation.direction.z*hc_indentation.dx;
			if(this->indentation->discreteSurf){
				generateMicaList_kernel<<<gsop.aminoCount/gsop.blockSize + 1, gsop.blockSize>>>();
			}
		}
		float4 dr;
		dr.x = hc_indentation.chipCoord.x - hc_indentation.tipCoord.x;
		dr.y = hc_indentation.chipCoord.y - hc_indentation.tipCoord.y;
		dr.z = hc_indentation.chipCoord.z - hc_indentation.tipCoord.z;
		dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		// Integrating equation of motion for cantilever tip (Langevin equation)
		float4 df;
		df.x = -f.x + dr.x*hc_indentation.cantileverKs;
		df.y = -f.y + dr.y*hc_indentation.cantileverKs;
		df.z = -f.z + dr.z*hc_indentation.cantileverKs;
		df.w = sqrtf(df.x*df.x + df.y*df.y + df.z*df.z);

		float mult = ((integrator->h*((float)nav))/(hc_indentation.tipZeta));

		if(hc_indentation.fixTransversal == 1){
			float acostheta =
					mult*df.x*hc_indentation.direction.x +
					mult*df.y*hc_indentation.direction.y +
					mult*df.z*hc_indentation.direction.z;
			hc_indentation.tipCoord.x += acostheta*hc_indentation.direction.x;
			hc_indentation.tipCoord.y += acostheta*hc_indentation.direction.y;
			hc_indentation.tipCoord.z += acostheta*hc_indentation.direction.z;
		} else {
			hc_indentation.tipCoord.x += mult*df.x;//hc_indentationChipCoordX;// + dx*hc_indentation.direction.x;
			hc_indentation.tipCoord.y += mult*df.y;//hc_indentationChipCoordY;// + dx*hc_indentation.direction.y;
			hc_indentation.tipCoord.z += mult*df.z;//hc_indentationChipCoordZ;// + dx*hc_indentation.direction.z;
		}
		// Moving cantilever dummy beads for representation purposes

		hc_indentation.chipCoordAv.x += hc_indentation.chipCoord.x;
		hc_indentation.chipCoordAv.y += hc_indentation.chipCoord.y;
		hc_indentation.chipCoordAv.z += hc_indentation.chipCoord.z;
		hc_indentation.tipCoordAv.x += hc_indentation.tipCoord.x;
		hc_indentation.tipCoordAv.y += hc_indentation.tipCoord.y;
		hc_indentation.tipCoordAv.z += hc_indentation.tipCoord.z;
		hc_indentation.fav.x += f.x;
		hc_indentation.fav.y += f.y;
		hc_indentation.fav.z += f.z;
		hc_indentation.fav.w += f.w;
		hc_indentation.kDeltaXAv += dr.w*hc_indentation.cantileverKs;
		// Output the resulting data on the screen and into the file
		if(step % hc_indentation.outputFreq == 0){
			int factor = hc_indentation.outputFreq/this->frequency;
			hc_indentation.chipCoordAv.x /= factor;
			hc_indentation.chipCoordAv.y /= factor;
			hc_indentation.chipCoordAv.z /= factor;
			hc_indentation.tipCoordAv.x /= factor;
			hc_indentation.tipCoordAv.y /= factor;
			hc_indentation.tipCoordAv.z /= factor;
			hc_indentation.fav.x /= factor;
			hc_indentation.fav.y /= factor;
			hc_indentation.fav.z /= factor;
			hc_indentation.fav.w /= factor;
			hc_indentation.kDeltaXAv /= factor;
			float fmod = sqrtf(hc_indentation.fav.x*hc_indentation.fav.x + hc_indentation.fav.y*hc_indentation.fav.y + hc_indentation.fav.z*hc_indentation.fav.z);
			float dxav = hc_indentation.dx - 0.5*factor*hc_indentation.V;
			printf("Indentation: %ld\t%5.2f\t%5.2f\t%5.2f\n", step, dxav, fmod, hc_indentation.kDeltaXAv);
			fprintf(hc_indentation.out, "%ld\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\n",
							step, hc_indentation.dx,
							hc_indentation.fav.w, fmod, abs(dr.w*hc_indentation.cantileverKs), hc_indentation.kDeltaXAv, hc_indentation.fav.x, hc_indentation.fav.y, hc_indentation.fav.z,
							hc_indentation.tipCoord.x, hc_indentation.tipCoord.y, hc_indentation.tipCoord.z,
							hc_indentation.chipCoord.x, hc_indentation.chipCoord.y, hc_indentation.chipCoord.z);
			fflush(hc_indentation.out);
			hc_indentation.chipCoordAv.x = 0.0f;
			hc_indentation.chipCoordAv.y = 0.0f;
			hc_indentation.chipCoordAv.z = 0.0f;
			hc_indentation.tipCoordAv.x = 0.0f;
			hc_indentation.tipCoordAv.y = 0.0f;
			hc_indentation.tipCoordAv.z = 0.0f;
			hc_indentation.fav.x = 0.0f;
			hc_indentation.fav.y = 0.0f;
			hc_indentation.fav.z = 0.0f;
			hc_indentation.fav.w = 0.0f;
			hc_indentation.kDeltaXAv = 0.0f;
		}
		// Update device data
		cudaMemcpyToSymbol(c_indentation, &hc_indentation, sizeof(IndentationConstant), 0, cudaMemcpyHostToDevice);
	}
}

void IndentationAminoUpdater::update(){
	if(step % this->frequency == 0){
		if(showTipMica){
			int i;
			sop.additionalAminos[1].x = hc_indentation.tipCoord.x;
			sop.additionalAminos[1].y = hc_indentation.tipCoord.y;
			sop.additionalAminos[1].z = hc_indentation.tipCoord.z;
			/*if(hc_indentation.moveSurface == 1){
					float displ = hc_indentation.surfaceStep;
					for(i = 0; i < hc_indentation.surfaceSize; i++){
						for(j = 0; j < hc_indentation.surfaceSize; j++){
							int i1 = i-hc_indentation.surfaceSize/2;
							int i2 = j-hc_indentation.surfaceSize/2;
							sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].x =
									hc_indentation.micaR.x + i1*displ*r1.x + i2*displ*r2.x;
							sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].y =
									hc_indentation.micaR.y + i1*displ*r1.y + i2*displ*r2.y;
							sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].z =
									hc_indentation.micaR.z + i1*displ*r1.z + i2*displ*r2.z;
						}
					}
				}*/


			if(hc_indentation.moveSurface == 1){
				for(i = 0; i < hc_indentation.surfaceBeadsCount; i++){
						/*float mult = hc_indentation.V*(this->frequency/hc_indentationUpdater.frequency);
						sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].x += mult*hc_indentation.direction.x;
						sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].y += mult*hc_indentation.direction.y;
						sop.additionalAminos[2+i*hc_indentation.surfaceSize+j].z += mult*hc_indentation.direction.z;*/

						sop.additionalAminos[2+i].x =
								hc_indentation.surfacePointsR0[i].x +
								hc_indentation.dx*hc_indentation.direction.x;
						sop.additionalAminos[2+i].y =
								hc_indentation.surfacePointsR0[i].y +
								hc_indentation.dx*hc_indentation.direction.y;
						sop.additionalAminos[2+i].z =
								hc_indentation.surfacePointsR0[i].z +
								hc_indentation.dx*hc_indentation.direction.z;

				}
			} else {
				sop.additionalAminos[0].x = hc_indentation.chipCoord.x + hc_indentation.cantileverVector.x;
				sop.additionalAminos[0].y = hc_indentation.chipCoord.y + hc_indentation.cantileverVector.y;
				sop.additionalAminos[0].z = hc_indentation.chipCoord.z + hc_indentation.cantileverVector.z;
			}
		}
	}
}

