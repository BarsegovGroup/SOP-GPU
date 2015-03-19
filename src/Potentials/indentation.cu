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

#include "indentation_kernel.cu"

/*
 * Add potential computation and update functions to lists,
 * create timer.
 */

void createIndentationPotential(){
	if(gsop.indentationOn == 1 || parameters::indentation.get()){
		if(gsop.Ntr > 1){
			DIE("Many-runs-per-GPU mode is not supported by indentation potential.\n");
		}
		gsop.indentationOn = 1;

		int showTipMica;
		int discreteSurf;
		showTipMica = parameters::indentationShowTipSurf.get();
		discreteSurf = parameters::indentationDiscreteSurf.get();
		if(discreteSurf && !showTipMica){
			printf("Discrete mica representation requires saving mica beads into dcd. Forcing '%s'.\n", parameters::indentationShowTipSurf.name().c_str());
			showTipMica = 1;
		}

		/*
		 * Creating updaters and potential
		 */

		IndentationTipUpdater* tipUpdater;
		IndentationAminoUpdater* aminoUpdater;
		IndentationPotential* pot;

		if(showTipMica){
			aminoUpdater = new IndentationAminoUpdater();
		}
		pot = new IndentationPotential();
		tipUpdater =  new IndentationTipUpdater(pot);

		/*
		 * Adding the created updaters and potential
		 */

		updaters[updatersCount] = tipUpdater;
		updatersCount++;

		if(showTipMica){
			updaters[updatersCount] = aminoUpdater;
			updatersCount++;
		} else {
			sop.additionalAminos.resize(0);
		}

		potentials[potentialsCount] = pot;
		potentialsCount++;

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

	printf("Initializing indentation protocol...\n");
	// Read parameters from configuration file

	// Cantilever parameters
	hc_indentation.chipCoord = parameters::indentationChip.get();
	hc_indentation.chipCoord0 = hc_indentation.chipCoord;
	hc_indentation.tipCoord = parameters::indentationTip.get();
    hc_indentation.direction = parameters::indentationDirection.get();
	hc_indentation.tipRadius = parameters::indentationTipR.get();
	hc_indentation.cantileverKs = parameters::indentationKs.get();
	hc_indentation.V = parameters::indentationDeltaX.get();
	hc_indentation.dx = 0;
	hc_indentation.tipZeta = parameters::indentationTipZeta.get();
	printf("Initial cantilever base coordinates: (%5.2f, %5.2f, %5.2f).\n", hc_indentation.chipCoord.x, hc_indentation.chipCoord.y, hc_indentation.chipCoord.z);

	// Surface parameters
	this->discreteSurf = parameters::indentationDiscreteSurf.get();
	hc_indentation.moveSurface = parameters::indentationMoveSurface.get();

	hc_indentation.micaR0 = parameters::indentationSurfaceR0.get();
	hc_indentation.micaR = hc_indentation.micaR0;
	hc_indentation.micaN = parameters::indentationSurfaceN.get();
	hc_indentation.fixTransversal = parameters::indentationFixTrans.get();

	// Tip-amino and surface-amino interactions
	float tipSigma, surfSigma, tipel, surfel;

	tipSigma = parameters::indentationTipSigma.get();
	surfSigma = parameters::indentationSurfSigma.get();

	tipel = parameters::indentationTipEl.get();
	surfel = parameters::indentationSurfEl.get();

	float tipA, tipB, surfA, surfB;
	tipA = parameters::indentationTipA.get();
	tipB = parameters::indentationTipB.get();
	surfA = parameters::indentationSurfA.get();
	surfB = parameters::indentationSurfB.get();

	hc_indentation.tipAprime  = 12.0f*tipel *powf( tipSigma, 12.0f)* tipA;
	hc_indentation.tipBprime  =  6.0f*tipel *powf( tipSigma,  6.0f)* tipB;
	hc_indentation.surfAprime = 12.0f*surfel*powf(surfSigma, 12.0f)*surfA;
	hc_indentation.surfBprime =  6.0f*surfel*powf(surfSigma,  6.0f)*surfB;

	// Normalize vectors
	float norm = sqrt(hc_indentation.direction.x*hc_indentation.direction.x + hc_indentation.direction.y*hc_indentation.direction.y + hc_indentation.direction.z*hc_indentation.direction.z);
	hc_indentation.direction.x /= norm;
	hc_indentation.direction.y /= norm;
	hc_indentation.direction.z /= norm;
	norm = sqrt(hc_indentation.micaN.x*hc_indentation.micaN.x + hc_indentation.micaN.y*hc_indentation.micaN.y + hc_indentation.micaN.z*hc_indentation.micaN.z);
	hc_indentation.micaN.x /= norm;
	hc_indentation.micaN.y /= norm;
	hc_indentation.micaN.z /= norm;

	// Allocate memory, copy data to the device
	hc_indentation.h_tipForces = (float4*)calloc(gsop.aminoCount, sizeof(float4));
	cudaMalloc((void**)&hc_indentation.d_tipForces, gsop.aminoCount*sizeof(float4));
	cudaMemcpy(hc_indentation.d_tipForces, hc_indentation.h_tipForces, gsop.aminoCount*sizeof(float4), cudaMemcpyHostToDevice);

	checkCUDAError();

	if(this->discreteSurf){

		hc_indentation.pairsCutoff2 = parameters::indentationPairsCutoff.get();
		hc_indentation.pairsCutoff2 = hc_indentation.pairsCutoff2*hc_indentation.pairsCutoff2;
		hc_indentation.surfaceBeadsCount = hc_indentation.surfaceSize*hc_indentation.surfaceSize;

		hc_indentation.h_surfacePointsCoord = (float4*)calloc(hc_indentation.surfaceSize*hc_indentation.surfaceSize,
					sizeof(float4));
		int i, j;
		for(i = 0; i < hc_indentation.surfaceSize; i++){
			for(j = 0; j < hc_indentation.surfaceSize; j++){
				hc_indentation.h_surfacePointsCoord[i*hc_indentation.surfaceSize + j].x =
						hc_indentation.surfacePointsR0[i*hc_indentation.surfaceSize + j].x; // These are defined in IndentationAminoUpdater initialization
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
	cudaMemcpyToSymbol(c_indentation, &hc_indentation, sizeof(hc_indentation), 0, cudaMemcpyHostToDevice);

	if(deviceProp.major == 2){ // TODO: >= 2
		cudaFuncSetCacheConfig(indentation_kernel, cudaFuncCachePreferL1);
	}
}

/*
 * An updater that shifts the cantilever tip or surface
 */
IndentationTipUpdater::IndentationTipUpdater(IndentationPotential *indentation) {
	this->name = "Indentation";
    this->frequency = gsop.nav;
    this->indentation = indentation;
    outputFreq = parameters::indentationOutputFreq.get();
    chipCoordAv.x = 0.0f;
	chipCoordAv.y = 0.0f;
	chipCoordAv.z = 0.0f;
	tipCoordAv.x = 0.0f;
	tipCoordAv.y = 0.0f;
	tipCoordAv.z = 0.0f;
	fav.x = 0.0f;
	fav.y = 0.0f;
	fav.z = 0.0f;
	fav.w = 0.0f;
	kDeltaXAv = 0.0f;
	// Output parameters
	this->outputFile = safe_fopen(parameters::indentationOutput.get().c_str(), "w");
	this->retractionStep = parameters::indentationRetractionStep.get();
}

/*
 * Create dummy particles, showing grid as a surface and a couple of beads for cantilever
 * Since this updater saves the information of the tip and surface position it should update after IndentationTipUpdater
 * Since it creates the surface beads, it should be initialized before IndentationPotential
 */
IndentationAminoUpdater::IndentationAminoUpdater(){
	this->name = "AdditionalAminos";
	this->frequency = parameters::outputtiming.get();

	/*
	 * We do re-read these parameters in indentation potential, though it is used only for the first frame and only
	 * if IndentationAminoUpdater is not updated before dcd output.
	 */
	float3 chipCoord;
	float3 tipCoord;
	chipCoord = parameters::indentationChip.get();
    tipCoord = parameters::indentationTip.get();
    hc_indentation.micaR0 = parameters::indentationSurfaceR0.get();
    hc_indentation.micaN = parameters::indentationSurfaceN.get();
    hc_indentation.direction = parameters::indentationDirection.get();

	// Generating additional amino-acids to represent cantilever and surface
	hc_indentation.surfaceSize = parameters::indentationSurfaceSize.get();
	hc_indentation.surfaceStep = parameters::indentationSurfaceStep.get();
	sop.additionalAminos.resize( 2 + hc_indentation.surfaceSize*hc_indentation.surfaceSize );

	// Cantilever base amino acid
	sop.additionalAminos[0].id = sop.aminos.size();
	strcpy(sop.additionalAminos[0].resName, "TIP");
	strcpy(sop.additionalAminos[0].name, "CA");
	sop.additionalAminos[0].chain = 'T';
	sop.additionalAminos[0].x = chipCoord.x;
	sop.additionalAminos[0].y = chipCoord.y;
	sop.additionalAminos[0].z = chipCoord.z;
	sop.additionalAminos[0].resid = 1;
	sop.additionalAminos[0].beta = sop.additionalAminos[0].occupancy = 0.0f;
	sop.additionalAminos[0].altLoc = ' ';
	// Cantilever tip amino acid
	sop.additionalAminos[1].id = sop.aminos.size() + 1;
	strcpy(sop.additionalAminos[1].resName, "TIP");
	strcpy(sop.additionalAminos[1].name, "CA");
	sop.additionalAminos[1].chain = 'T';
	sop.additionalAminos[1].x = tipCoord.x;
	sop.additionalAminos[1].y = tipCoord.y;
	sop.additionalAminos[1].z = tipCoord.z;
	sop.additionalAminos[1].resid = 2;
	sop.additionalAminos[1].beta = sop.additionalAminos[1].occupancy = 0.0f;
	sop.additionalAminos[1].altLoc = ' ';

	// Shifting cantilever base amino for representation purpose
	// Finding vector, perpendicular to indentation direction
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
	// Shifting cantilever base amino along the found vector
	len = parameters::indentationCantLength.get();
	r.x *= len;
	r.y *= len;
	r.z *= len;
	printf("Cantilever representation vector: (%5.3f, %5.3f, %5.3f)\n", r.x, r.y, r.z);
	cantileverVector = r;
	sop.additionalAminos[0].x += r.x;
	sop.additionalAminos[0].y += r.y;
	sop.additionalAminos[0].z += r.z;

	// Generating amino-acids to represent the surface
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
			int amino_id = 2+i*hc_indentation.surfaceSize+j;
			sop.additionalAminos[amino_id].id =
					sop.aminos.size() + 2 + i*hc_indentation.surfaceSize+j;
			strcpy(sop.additionalAminos[amino_id].resName, "MIC");
			strcpy(sop.additionalAminos[amino_id].name, "CA");
			sop.additionalAminos[amino_id].chain = 'M';
			sop.additionalAminos[amino_id].resid = 1;
			sop.additionalAminos[amino_id].x = r.x;
			sop.additionalAminos[amino_id].y = r.y;
			sop.additionalAminos[amino_id].z = r.z;
			sop.additionalAminos[amino_id].beta = sop.additionalAminos[amino_id].occupancy = 0.0f;
			sop.additionalAminos[amino_id].altLoc = ' ';
			hc_indentation.surfacePointsR0[i*hc_indentation.surfaceSize + j] = r; // These will use in case discrete surface representation is selected
		}
	}

	/*
	 *  VMD script to connect surface beads (for representation purposes)
	 */
    std::string vmdFilename = parameters::indentationSurfConnectFile.get();
	FILE* micaConnectFile = safe_fopen(vmdFilename.c_str(), "w");
	printf("Dumping VMD script to connect surface elements into '%s'.\n", vmdFilename.c_str());
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

/*
 * Move cantilever or surface
 */
void IndentationTipUpdater::update(){
	// Change moving direction to switch to retraction
	if(gsop.step == this->retractionStep){
		printf("Starting tip retraction at step %ld.\n", this->retractionStep);
		/*hc_indentation.direction.x = -hc_indentation.direction.x;
		hc_indentation.direction.y = -hc_indentation.direction.y;
		hc_indentation.direction.z = -hc_indentation.direction.z;*/
		hc_indentation.V = -hc_indentation.V;
	}
	if(gsop.step % this->frequency == 0){
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

		float mult = ((integrator->h*((float)gsop.nav))/(hc_indentation.tipZeta));

		if(hc_indentation.fixTransversal == 1){
			float costheta =
					mult*df.x*hc_indentation.direction.x +
					mult*df.y*hc_indentation.direction.y +
					mult*df.z*hc_indentation.direction.z;
			hc_indentation.tipCoord.x += costheta*hc_indentation.direction.x;
			hc_indentation.tipCoord.y += costheta*hc_indentation.direction.y;
			hc_indentation.tipCoord.z += costheta*hc_indentation.direction.z;
		} else {
			hc_indentation.tipCoord.x += mult*df.x;//hc_indentationChipCoordX;// + dx*hc_indentation.direction.x;
			hc_indentation.tipCoord.y += mult*df.y;//hc_indentationChipCoordY;// + dx*hc_indentation.direction.y;
			hc_indentation.tipCoord.z += mult*df.z;//hc_indentationChipCoordZ;// + dx*hc_indentation.direction.z;
		}
		// Moving cantilever dummy beads for representation purposes

		chipCoordAv.x += hc_indentation.chipCoord.x;
		chipCoordAv.y += hc_indentation.chipCoord.y;
		chipCoordAv.z += hc_indentation.chipCoord.z;
		tipCoordAv.x += hc_indentation.tipCoord.x;
		tipCoordAv.y += hc_indentation.tipCoord.y;
		tipCoordAv.z += hc_indentation.tipCoord.z;
		fav.x += f.x;
		fav.y += f.y;
		fav.z += f.z;
		fav.w += f.w;
		kDeltaXAv += dr.w*hc_indentation.cantileverKs;
		// Output the resulting data on the screen and into the file
		if(gsop.step % outputFreq == 0){
			int factor = outputFreq/this->frequency;
			chipCoordAv.x /= factor;
			chipCoordAv.y /= factor;
			chipCoordAv.z /= factor;
			tipCoordAv.x /= factor;
			tipCoordAv.y /= factor;
			tipCoordAv.z /= factor;
			fav.x /= factor;
			fav.y /= factor;
			fav.z /= factor;
			fav.w /= factor;
			kDeltaXAv /= factor;
			float fmod = sqrtf(fav.x*fav.x + fav.y*fav.y + fav.z*fav.z);
			float dxav = hc_indentation.dx - 0.5*factor*hc_indentation.V;
			printf("Indentation: %ld\t%5.2f\t%5.2f\t%5.2f\n", gsop.step, dxav, fmod, kDeltaXAv);
			fprintf(this->outputFile, "%ld\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\t"
					"%8.5f\t%8.5f\t%8.5f\n",
							gsop.step, hc_indentation.dx,
							fav.w, fmod, abs(dr.w*hc_indentation.cantileverKs), kDeltaXAv, fav.x, fav.y, fav.z,
							hc_indentation.tipCoord.x, hc_indentation.tipCoord.y, hc_indentation.tipCoord.z,
							hc_indentation.chipCoord.x, hc_indentation.chipCoord.y, hc_indentation.chipCoord.z);
			fflush(this->outputFile);
			chipCoordAv.x = 0.0f;
			chipCoordAv.y = 0.0f;
			chipCoordAv.z = 0.0f;
			tipCoordAv.x = 0.0f;
			tipCoordAv.y = 0.0f;
			tipCoordAv.z = 0.0f;
			fav.x = 0.0f;
			fav.y = 0.0f;
			fav.z = 0.0f;
			fav.w = 0.0f;
			kDeltaXAv = 0.0f;
		}
		// Update device data
		cudaMemcpyToSymbol(c_indentation, &hc_indentation, sizeof(IndentationConstant), 0, cudaMemcpyHostToDevice);
	}
}

/*
 * Update cantilever tip/surface position
 */
void IndentationAminoUpdater::update(){
	if(gsop.step % this->frequency == 0){
		int i;
		sop.additionalAminos[1].x = hc_indentation.tipCoord.x;
		sop.additionalAminos[1].y = hc_indentation.tipCoord.y;
		sop.additionalAminos[1].z = hc_indentation.tipCoord.z;

		if(hc_indentation.moveSurface == 1){
			for(i = 0; i < hc_indentation.surfaceBeadsCount; i++){
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
			sop.additionalAminos[0].x = hc_indentation.chipCoord.x + cantileverVector.x;
			sop.additionalAminos[0].y = hc_indentation.chipCoord.y + cantileverVector.y;
			sop.additionalAminos[0].z = hc_indentation.chipCoord.z + cantileverVector.z;
		}
	}
}

