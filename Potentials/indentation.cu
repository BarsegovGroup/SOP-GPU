/*
 * indentation.cu
 *
 *  Created on: Apr 9, 2010
 *      Author: zhmurov
 */
#include "../gsop.cuh"
#include "indentation.cuh"
#include "../Updaters/output_manager.cuh"
#include "indentation_kernel.cu"

int showTipMica;
int discreteSurf;

void addTipMicaParticles();

/*
 * Add potential computation and update functions to lists,
 * create timer.
 */

void createIndentationPotential(){
	if(indentationOn == 1 || getYesNoParameter(INDENTATION_ON_STRING, 0, 1) == 1){
		if(Ntr > 1){
			printf("Many-runs-per-GPU mode is not supported by indentation potential.\n");
			exit(0);
		}
		indentationOn = 1;
		sprintf(indentationPotential.name, "Indentation");
		discreteSurf = getYesNoParameter(INDENTATION_DISCRETE_SURF_STRING, 0, 1);
		if(discreteSurf){
			indentationPotential.compute = &computeIndentationDiscreteSurface;
		} else {
			indentationPotential.compute = &computeIndentation;
		}
		indentationPotential.computeEnergy = &computeIndentationEnergy;
		potentials[potentialsCount] = &indentationPotential;
		if(gsop.deviceProp.major == 2){
			cudaFuncSetCacheConfig(indentation_kernel, cudaFuncCachePreferL1);
		}
		potentialsCount++;
		sprintf(indentationUpdater.name, "Indentation");
		indentationUpdater.update = &updateTipPosition;
		indentationUpdater.frequency = nav;
		updaters[updatersCount] = &indentationUpdater;
		updatersCount++;
		initIndentation();
	}
}

/*
 * Initialize potential:
 * Tip and surface are both represented by a repulsive LJ-potential.
 * Tip is a sphere, connected to cantilever base by a spring.
 * Motion of a tip is described by Langevin equation but not integrated every step,
 * since tip is much heavier then particles.
 */
void initIndentation(){
	printf("Initializing indentation protocol...\n");
	// Read poteintial from configuration file
	getVectorParameter(INDENTATION_CHIP_POSITION_STRING, &indentation.chipCoord.x, &indentation.chipCoord.y, &indentation.chipCoord.z, 0.0, 0.0, 0.0, 0);
	indentation.chipCoord0 = indentation.chipCoord;
	getVectorParameter(INDENTATION_TIP_POSITION_STRING, &indentation.tipCoord.x, &indentation.tipCoord.y, &indentation.tipCoord.z,
			indentation.chipCoord.x, indentation.chipCoord.y, indentation.chipCoord.z, 1);
	getVectorParameter(INDENTATION_DIRECTION_STRING, &indentation.direction.x, &indentation.direction.y, &indentation.direction.z, 0.0, 0.0, 0.0, 0);
	getVectorParameter(INDENTATION_SURFACE_R0_STRING, &indentation.micaR0.x, &indentation.micaR0.y, &indentation.micaR0.z, 0.0, 0.0, 0.0, 0);
	indentation.micaR = indentation.micaR0;
	getVectorParameter(INDENTATION_SURFACE_N_STRING, &indentation.micaN.x, &indentation.micaN.y, &indentation.micaN.z, 0.0, 0.0, 0.0, 0);
	indentation.tipRadius = getFloatParameter(INDENTATION_TIP_RADIUS_STRING, 0.0, 0);
	indentation.cantileverKs = getFloatParameter(INDENTATION_KS_STRING, 0.0, 0);
	indentation.V = getFloatParameter(INDENTATION_DELTAX_STRING, 0.0, 0);
	indentation.dx = 0;
	indentation.moveSurface = getYesNoParameter(INDENTATION_MOVE_SURFACE, 0, 1);
	indentation.fixTransversal = getYesNoParameter(INDENTATION_FIX_TRANSVERSAL,
			DEFAULT_INDENTATION_FIX_TRANSVERSAL, 1);
	float tempVar;
	tempVar = getFloatParameter(INDENTATION_SIGMA_STRING, 1.0, 1);
	indentation.tipa6 = powf(getFloatParameter(INDENTATION_TIP_SIGMA_STRING, tempVar, 1), 6.0);
	indentation.surfa6 = powf(getFloatParameter(INDENTATION_SURF_SIGMA_STRING, tempVar, 1), 6.0);
	tempVar = getFloatParameter(INDENTATION_EL_STRING, 1.0, 1);
	indentation.tipel = getFloatParameter(INDENTATION_TIP_EL_STRING, tempVar, 1);
	indentation.surfel = getFloatParameter(INDENTATION_SURF_EL_STRING, tempVar, 1);

	float tipA, tipB, surfA, surfB;
	tipA = getFloatParameter(INDENTATION_TIP_A, DEFAULT_INDENTATION_TIP_A, 1);
	tipB = getFloatParameter(INDENTATION_TIP_B, DEFAULT_INDENTATION_TIP_B, 1);
	surfA = getFloatParameter(INDENTATION_SURFACE_A, DEFAULT_INDENTATION_SURFACE_A, 1);
	surfB = getFloatParameter(INDENTATION_SURFACE_B, DEFAULT_INDENTATION_SURFACE_B, 1);
	indentation.tipAprime = 12.0f*indentation.tipel*indentation.tipa6*indentation.tipa6*tipA;
	indentation.tipBprime = 6.0f*indentation.tipel*indentation.tipa6*tipB;
	indentation.surfAprime = 12.0f*indentation.surfel*indentation.surfa6*indentation.surfa6*surfA;
	indentation.surfBprime = 6.0f*indentation.surfel*indentation.surfa6*surfB;

	/*float x = 201;
	FILE* tempFile = fopen("temp.dat", "w");
	while(x < 300){
		float xprime = pow(x - indentation.tipRadius, 6.0);
		float y = indentation.el*((tipA*indentation.a6*indentation.a6)/(xprime*xprime) + (tipB*indentation.a6)/(xprime));
		fprintf(tempFile, "%f\t%f\n", x, y);
		x += 0.1;
	}
	fclose(tempFile);
	exit(0);*/

	indentation.tipZeta = getFloatParameter(INDENTATION_TIP_ZETA_STRING, 5000.0, 1);
	char indentation_filename[100];
	getMaskedParameter(indentation_filename, INDENTATION_OUTPUT_FILENAME_STRING, DEFAULT_INDENTATION_FILENAME, 1);

	printf("Initial cantilever base coordinates: (%5.2f, %5.2f, %5.2f).\n", indentation.chipCoord.x, indentation.chipCoord.y, indentation.chipCoord.z);

	// Normilize vectors
	float norm = sqrt(indentation.direction.x*indentation.direction.x + indentation.direction.y*indentation.direction.y + indentation.direction.z*indentation.direction.z);
	indentation.direction.x /= norm;
	indentation.direction.y /= norm;
	indentation.direction.z /= norm;
	norm = sqrt(indentation.micaN.x*indentation.micaN.x + indentation.micaN.y*indentation.micaN.y + indentation.micaN.z*indentation.micaN.z);
	indentation.micaN.x /= norm;
	indentation.micaN.y /= norm;
	indentation.micaN.z /= norm;

	// Allocatae memory, copy data to the device
	indentation.h_tipForces = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMalloc((void**)&indentation.d_tipForces, gsop.aminoCount*sizeof(float4));
	cudaMemcpy(indentation.d_tipForces, indentation.h_tipForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
	indentation.out = fopen(indentation_filename, "w");
	indentation.retractionStep = getLongIntegerParameter(INDENTATION_RETRACTION_STEP_STRING, -1, 1);
	checkCUDAError();

	// Add dummy particles to show surface and cantilever in dcd
	showTipMica = getYesNoParameter(INDENTATION_SHOW_TIP_SURFACE_STRING, 0, 1);
	if(!discreteSurf){
		printf("Discrete mica representation requires saving mica beads into dcd. Forcing \"showTipSurf\".\n");
		showTipMica = 1;
	}

	if(showTipMica){
		sprintf(additionalAminosUpdater.name, "AdditionalAminos");
		additionalAminosUpdater.update = &updateAdditionalAminos;
		additionalAminosUpdater.frequency = getIntegerParameter(OUTPUT_FREQUENCY_STRING, DEFAULT_OUTPUT_FREQUENCY, 1);
		updaters[updatersCount] = &additionalAminosUpdater;
		updatersCount++;
		addTipMicaParticles();
	} else {
		sop.additionalAminosCount = 0;
	}
	if(discreteSurf){

		indentation.pairsCutoff2 = getFloatParameter(INDENTATION_PAIRS_CUTOFF_STRING,
				DEFAULT_INDENTATION_PAIRS_CUTOFF, 1);
		indentation.pairsCutoff2 = indentation.pairsCutoff2*indentation.pairsCutoff2;
		indentation.surfaceBeadsCount = indentation.surfaceSize*indentation.surfaceSize;

		indentation.h_surfacePointsCoord = (float4*)calloc(indentation.surfaceSize*indentation.surfaceSize,
					sizeof(float4));
		int i, j;
		for(i = 0; i < indentation.surfaceSize; i++){
			for(j = 0; j < indentation.surfaceSize; j++){
				indentation.h_surfacePointsCoord[i*indentation.surfaceSize + j].x =
						indentation.surfacePointsR0[i*indentation.surfaceSize + j].x;
				indentation.h_surfacePointsCoord[i*indentation.surfaceSize + j].y =
						indentation.surfacePointsR0[i*indentation.surfaceSize + j].y;
				indentation.h_surfacePointsCoord[i*indentation.surfaceSize + j].z =
						indentation.surfacePointsR0[i*indentation.surfaceSize + j].z;
			}
		}


		printf("Computing maximum pairs in mica pairslist..\n");
		cudaMalloc((void**)&indentation.d_surfacePointsCoord, indentation.surfaceSize*indentation.surfaceSize*sizeof(float4));

		int maxPairs = 0;
		cudaMallocHost((void**)&indentation.h_micaListCounts, gsop.aminoCount*sizeof(int));
		for(i = 0; i < gsop.aminoCount; i++){
			indentation.h_micaListCounts[i] = 0;
		}
		float4 r1, r2;
		for(i = 0; i < gsop.aminoCount; i++){
			for(j = 0; j < indentation.surfaceBeadsCount; j++){
				r1 = gsop.h_coord[i];
				r2 = indentation.h_surfacePointsCoord[j];
				r1.x -= r2.x;
				r1.y -= r2.y;
				r1.z -= r2.z;
				r1.w = r1.x*r1.x + r1.y*r1.y + r1.z*r1.z;
				if(r1.w < indentation.pairsCutoff2){
					indentation.h_micaListCounts[i] ++;
				}
			}
		}
		for(i = 0; i < gsop.aminoCount; i++){
			if(indentation.h_micaListCounts[i] > maxPairs){
				maxPairs = indentation.h_micaListCounts[i];
			}
		}
		printf("Done. Maximum pairs found is %d.\n", maxPairs);
		maxPairs += 256;
		if(maxPairs > indentation.surfaceBeadsCount){
			maxPairs = indentation.surfaceBeadsCount;
		}


		cudaMallocHost((void**)&indentation.h_micaListCounts, gsop.aminoCount*sizeof(int));
		cudaMallocHost((void**)&indentation.h_micaList,	maxPairs*gsop.width*sizeof(int));
		cudaMalloc((void**)&indentation.d_micaListCounts, gsop.aminoCount*sizeof(int));
		cudaMalloc((void**)&indentation.d_micaList, maxPairs*gsop.width*sizeof(int));

		cudaMemcpy(indentation.d_surfacePointsCoord, indentation.h_surfacePointsCoord,
				indentation.surfaceSize*indentation.surfaceSize*sizeof(float4), cudaMemcpyHostToDevice);
	}
	indentation.outputFreq = getIntegerParameter(INDENTATION_OUTPUT_FREQ_STRING, DEFAULT_INDENTATION_OUTPUT_FREQ, 1);
	cudaMemcpyToSymbol(c_indentation, &indentation, sizeof(Indentation), 0, cudaMemcpyHostToDevice);
}

/*
 * Create dummy particles, showing grid as a surface and a couple of beads for cantilever
 */
void addTipMicaParticles(){
	indentation.surfaceSize = getIntegerParameter(INDENTATION_MICA_SIZE_STRING, 51, 1);
	indentation.surfaceStep = getFloatParameter(INDENTATION_MICA_STEP_STRING, 1.4, 1);
	sop.additionalAminosCount = 2 + indentation.surfaceSize*indentation.surfaceSize;
	sop.additionalAminos = (Atom*)calloc(sop.additionalAminosCount, sizeof(Atom));
	sop.additionalAminos[0].id = sop.aminoCount;
	strcpy(sop.additionalAminos[0].resName, "TIP");
	sop.additionalAminos[0].chain = 'T';
	sop.additionalAminos[0].x = indentation.chipCoord.x;
	sop.additionalAminos[0].y = indentation.chipCoord.y;
	sop.additionalAminos[0].z = indentation.chipCoord.z;
	sop.additionalAminos[0].resid = 1;
	sop.additionalAminos[1].id = sop.aminoCount + 1;
	strcpy(sop.additionalAminos[1].resName, "TIP");
	sop.additionalAminos[1].chain = 'T';
	sop.additionalAminos[1].x = indentation.tipCoord.x;
	sop.additionalAminos[1].y = indentation.tipCoord.y;
	sop.additionalAminos[1].z = indentation.tipCoord.z;
	sop.additionalAminos[1].resid = 2;
	float3 r1, r2, n, r0;
	float3 r;
	float len;
	r0 = indentation.direction;
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
	indentation.cantileverVector = r;
	sop.additionalAminos[0].x += r.x;
	sop.additionalAminos[0].y += r.y;
	sop.additionalAminos[0].z += r.z;

	indentation.surfacePointsR0 = (float3*)calloc(indentation.surfaceSize*indentation.surfaceSize,
			sizeof(float3));
	n = indentation.micaN;
	r0 = indentation.micaR0;
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
	float displ = indentation.surfaceStep;
	for(i = 0; i < indentation.surfaceSize; i++){
		for(j = 0; j < indentation.surfaceSize; j++){
			int i1 = i-indentation.surfaceSize/2;
			int i2 = j-indentation.surfaceSize/2;
			r.x = r0.x + i1*displ*r1.x + i2*displ*r2.x;
			r.y = r0.y + i1*displ*r1.y + i2*displ*r2.y;
			r.z = r0.z + i1*displ*r1.z + i2*displ*r2.z;
			sop.additionalAminos[2+i*indentation.surfaceSize+j].id =
					sop.aminoCount + 2 + i*indentation.surfaceSize+j;
			strcpy(sop.additionalAminos[2+i*indentation.surfaceSize+j].resName, "MIC");
			sop.additionalAminos[2+i*indentation.surfaceSize+j].chain = 'M';
			sop.additionalAminos[2+i*indentation.surfaceSize+j].x = r.x;
			sop.additionalAminos[2+i*indentation.surfaceSize+j].y = r.y;
			sop.additionalAminos[2+i*indentation.surfaceSize+j].z = r.z;
			indentation.surfacePointsR0[i*indentation.surfaceSize + j] = r;
		}
	}
	indentation.chipCoordAv.x = 0.0f;
	indentation.chipCoordAv.y = 0.0f;
	indentation.chipCoordAv.z = 0.0f;
	indentation.tipCoordAv.x = 0.0f;
	indentation.tipCoordAv.y = 0.0f;
	indentation.tipCoordAv.z = 0.0f;
	indentation.fav.x = 0.0f;
	indentation.fav.y = 0.0f;
	indentation.fav.z = 0.0f;
	indentation.fav.w = 0.0f;
	char tempFilename[100];
	getMaskedParameter(tempFilename, INDENTATION_VMD_CONNECT_SCRIPT, "connect_mica.vmd", 1);
	FILE* micaConnectFile = fopen(tempFilename, "w");
	printf("Dumping VMD script to connect surface elements into '%s'.\n", tempFilename);
	fprintf(micaConnectFile, "set sel [atomselect top \"index %d\"]\n", gsop.aminoCount);
	fprintf(micaConnectFile, "set bonds {%d}\n", gsop.aminoCount + 1);
	fprintf(micaConnectFile, "$sel setbonds [list $bonds]\n$sel delete\n");
	fprintf(micaConnectFile, "set sel [atomselect top \"index %d\"]\n", gsop.aminoCount + 1);
	fprintf(micaConnectFile, "set bonds {%d}\n", gsop.aminoCount);
	fprintf(micaConnectFile, "$sel setbonds [list $bonds]\n$sel delete\n");
	for(i = 0; i < indentation.surfaceSize; i++){
		for(j = 0; j < indentation.surfaceSize; j++){
			int i1 = 2+i*indentation.surfaceSize+j + gsop.aminoCount;
			int i2;
			fprintf(micaConnectFile, "set sel [atomselect top \"index %d\"]\n", i1);
			fprintf(micaConnectFile, "set bonds {");
			if(i != indentation.surfaceSize - 1){
				i2 = 2+(i+1)*indentation.surfaceSize+j + gsop.aminoCount;
				fprintf(micaConnectFile, "%d", i2);
			}
			if(j != indentation.surfaceSize - 1){
				i2 = 2+i*indentation.surfaceSize+j+1 + gsop.aminoCount;
				fprintf(micaConnectFile, " %d", i2);
			}
			fprintf(micaConnectFile, "}\n", i2);
			fprintf(micaConnectFile, "$sel setbonds [list $bonds]\n$sel delete\n");

		}
	}
	fclose(micaConnectFile);
}

/*
 * Execute computations on GPU
 */
inline void computeIndentation(){
	indentation_kernel<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
}

inline void computeIndentationDiscreteSurface(){
	indentationDiscreteSurf_kernel<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
}

/*
 * Stub for energy evaluation
 */
inline void computeIndentationEnergy(){
	// TODO (if needed)
}


/*
 * Move cantilever
 */
inline void updateTipPosition(){
	// Change moving direction to switch to retraction
	if(step == indentation.retractionStep){
		printf("Starting tip retraction at step %ld.\n", indentation.retractionStep);
		/*indentation.direction.x = -indentation.direction.x;
		indentation.direction.y = -indentation.direction.y;
		indentation.direction.z = -indentation.direction.z;*/
		indentation.V = -indentation.V;
	}
	if(step % indentationUpdater.frequency == 0){
		// Copy forces acting onto tip to host memory
		cudaMemcpy(indentation.h_tipForces, indentation.d_tipForces, gsop.aminoCount*sizeof(float4), cudaMemcpyDeviceToHost);
		int i;
		float4 f = make_float4(0.0, 0.0, 0.0, 0.0);
		// Sum up all force and set 0 value to forces
		for(i = 0; i < gsop.aminoCount; i++){
			f.x += indentation.h_tipForces[i].x;
			f.y += indentation.h_tipForces[i].y;
			f.z += indentation.h_tipForces[i].z;
			indentation.h_tipForces[i].x = 0.0f;
			indentation.h_tipForces[i].y = 0.0f;
			indentation.h_tipForces[i].z = 0.0f;
		}
		cudaMemcpy(indentation.d_tipForces, indentation.h_tipForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);
		// Compute average force and its absolute value
		f.x /= (float)indentationUpdater.frequency;
		f.y /= (float)indentationUpdater.frequency;
		f.z /= (float)indentationUpdater.frequency;
		f.w = sqrtf(f.x*f.x + f.y*f.y + f.z*f.z);
		indentation.tipForce = f;
		// Move chip or surface
		indentation.dx += indentation.V;//(step/nav)*indentation.chipV;
		if(indentation.moveSurface == 1){
			if(discreteSurf){
				for(i = 0; i < indentation.surfaceBeadsCount; i++){
					indentation.h_surfacePointsCoord[i].x = indentation.surfacePointsR0[i].x + indentation.direction.x*indentation.dx;
					indentation.h_surfacePointsCoord[i].y = indentation.surfacePointsR0[i].y + indentation.direction.y*indentation.dx;
					indentation.h_surfacePointsCoord[i].z = indentation.surfacePointsR0[i].z + indentation.direction.z*indentation.dx;
				}
				cudaMemcpy(indentation.d_surfacePointsCoord, indentation.h_surfacePointsCoord,
						indentation.surfaceBeadsCount*sizeof(float4), cudaMemcpyHostToDevice);
				generateMicaList_kernel<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
			} else {
				indentation.micaR.x = indentation.micaR0.x + indentation.direction.x*indentation.dx;
				indentation.micaR.y = indentation.micaR0.y + indentation.direction.y*indentation.dx;
				indentation.micaR.z = indentation.micaR0.z + indentation.direction.z*indentation.dx;
			}
		} else {
			indentation.chipCoord.x = indentation.chipCoord0.x + indentation.direction.x*indentation.dx;
			indentation.chipCoord.y = indentation.chipCoord0.y + indentation.direction.y*indentation.dx;
			indentation.chipCoord.z = indentation.chipCoord0.z + indentation.direction.z*indentation.dx;
			if(discreteSurf){
				generateMicaList_kernel<<<gsop.aminoCount/BLOCK_SIZE + 1, BLOCK_SIZE>>>();
			}
		}
		float4 dr;
		dr.x = indentation.chipCoord.x - indentation.tipCoord.x;
		dr.y = indentation.chipCoord.y - indentation.tipCoord.y;
		dr.z = indentation.chipCoord.z - indentation.tipCoord.z;
		dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		// Integrating equation of motion for cantilever tip (Langevin equation)
		float4 df;
		df.x = -f.x + dr.x*indentation.cantileverKs;
		df.y = -f.y + dr.y*indentation.cantileverKs;
		df.z = -f.z + dr.z*indentation.cantileverKs;
		df.w = sqrtf(df.x*df.x + df.y*df.y + df.z*df.z);

		float mult = ((integrator->h*((float)nav))/(indentation.tipZeta));

		if(indentation.fixTransversal == 1){
			float acostheta =
					mult*df.x*indentation.direction.x +
					mult*df.y*indentation.direction.y +
					mult*df.z*indentation.direction.z;
			indentation.tipCoord.x += acostheta*indentation.direction.x;
			indentation.tipCoord.y += acostheta*indentation.direction.y;
			indentation.tipCoord.z += acostheta*indentation.direction.z;
		} else {
			indentation.tipCoord.x += mult*df.x;//indentationChipCoordX;// + dx*indentation.direction.x;
			indentation.tipCoord.y += mult*df.y;//indentationChipCoordY;// + dx*indentation.direction.y;
			indentation.tipCoord.z += mult*df.z;//indentationChipCoordZ;// + dx*indentation.direction.z;
		}
		// Moving cantilever dummy beads for representation purposes

		indentation.chipCoordAv.x += indentation.chipCoord.x;
		indentation.chipCoordAv.y += indentation.chipCoord.y;
		indentation.chipCoordAv.z += indentation.chipCoord.z;
		indentation.tipCoordAv.x += indentation.tipCoord.x;
		indentation.tipCoordAv.y += indentation.tipCoord.y;
		indentation.tipCoordAv.z += indentation.tipCoord.z;
		indentation.fav.x += f.x;
		indentation.fav.y += f.y;
		indentation.fav.z += f.z;
		indentation.fav.w += f.w;
		indentation.kDeltaXAv += dr.w*indentation.cantileverKs;
		// Output the resulting data on the screen and into the file
		if(step % indentation.outputFreq == 0){
			int factor = indentation.outputFreq/indentationUpdater.frequency;
			indentation.chipCoordAv.x /= factor;
			indentation.chipCoordAv.y /= factor;
			indentation.chipCoordAv.z /= factor;
			indentation.tipCoordAv.x /= factor;
			indentation.tipCoordAv.y /= factor;
			indentation.tipCoordAv.z /= factor;
			indentation.fav.x /= factor;
			indentation.fav.y /= factor;
			indentation.fav.z /= factor;
			indentation.fav.w /= factor;
			indentation.kDeltaXAv /= factor;
			float fmod = sqrtf(indentation.fav.x*indentation.fav.x + indentation.fav.y*indentation.fav.y + indentation.fav.z*indentation.fav.z);
			float dxav = indentation.dx - 0.5*factor*indentation.V;
			printf("Indentation: %ld\t%5.2f\t%5.2f\t%5.2f\n", step, dxav, fmod, indentation.kDeltaXAv);
			fprintf(indentation.out, "%ld\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\n",
							step, indentation.dx,
							indentation.fav.w, fmod, abs(dr.w*indentation.cantileverKs), indentation.kDeltaXAv, indentation.fav.x, indentation.fav.y, indentation.fav.z,
							indentation.tipCoord.x, indentation.tipCoord.y, indentation.tipCoord.z,
							indentation.chipCoord.x, indentation.chipCoord.y, indentation.chipCoord.z);
			fflush(indentation.out);
			indentation.chipCoordAv.x = 0.0f;
			indentation.chipCoordAv.y = 0.0f;
			indentation.chipCoordAv.z = 0.0f;
			indentation.tipCoordAv.x = 0.0f;
			indentation.tipCoordAv.y = 0.0f;
			indentation.tipCoordAv.z = 0.0f;
			indentation.fav.x = 0.0f;
			indentation.fav.y = 0.0f;
			indentation.fav.z = 0.0f;
			indentation.fav.w = 0.0f;
			indentation.kDeltaXAv = 0.0f;
		}
		// Update device data
		cudaMemcpyToSymbol(c_indentation, &indentation, sizeof(Indentation), 0, cudaMemcpyHostToDevice);
	}
}

void updateAdditionalAminos(){
	if(step % additionalAminosUpdater.frequency == 0){
		if(showTipMica){
			int i;
			sop.additionalAminos[1].x = indentation.tipCoord.x;
			sop.additionalAminos[1].y = indentation.tipCoord.y;
			sop.additionalAminos[1].z = indentation.tipCoord.z;
			/*if(indentation.moveSurface == 1){
					float displ = indentation.surfaceStep;
					for(i = 0; i < indentation.surfaceSize; i++){
						for(j = 0; j < indentation.surfaceSize; j++){
							int i1 = i-indentation.surfaceSize/2;
							int i2 = j-indentation.surfaceSize/2;
							sop.additionalAminos[2+i*indentation.surfaceSize+j].x =
									indentation.micaR.x + i1*displ*r1.x + i2*displ*r2.x;
							sop.additionalAminos[2+i*indentation.surfaceSize+j].y =
									indentation.micaR.y + i1*displ*r1.y + i2*displ*r2.y;
							sop.additionalAminos[2+i*indentation.surfaceSize+j].z =
									indentation.micaR.z + i1*displ*r1.z + i2*displ*r2.z;
						}
					}
				}*/


			if(indentation.moveSurface == 1){
				for(i = 0; i < indentation.surfaceBeadsCount; i++){
						/*float mult = indentation.V*(additionalAminosUpdater.frequency/indentationUpdater.frequency);
						sop.additionalAminos[2+i*indentation.surfaceSize+j].x += mult*indentation.direction.x;
						sop.additionalAminos[2+i*indentation.surfaceSize+j].y += mult*indentation.direction.y;
						sop.additionalAminos[2+i*indentation.surfaceSize+j].z += mult*indentation.direction.z;*/

						sop.additionalAminos[2+i].x =
								indentation.surfacePointsR0[i].x +
								indentation.dx*indentation.direction.x;
						sop.additionalAminos[2+i].y =
								indentation.surfacePointsR0[i].y +
								indentation.dx*indentation.direction.y;
						sop.additionalAminos[2+i].z =
								indentation.surfacePointsR0[i].z +
								indentation.dx*indentation.direction.z;

				}
			} else {
				sop.additionalAminos[0].x = indentation.chipCoord.x + indentation.cantileverVector.x;
				sop.additionalAminos[0].y = indentation.chipCoord.y + indentation.cantileverVector.y;
				sop.additionalAminos[0].z = indentation.chipCoord.z + indentation.cantileverVector.z;
			}
		}
	}
}
