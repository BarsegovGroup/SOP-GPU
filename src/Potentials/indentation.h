/*
 * indentation.cuh
 *
 *  Created on: Apr 9, 2010
 *      Author: zhmurov
 */

#pragma once

#define INDENTATION_ON_STRING				"indentation"

#define INDENTATION_CHIP_POSITION_STRING	"indentationChip"
#define INDENTATION_TIP_POSITION_STRING		"indentationTip"
#define INDENTATION_DIRECTION_STRING		"indentationDirection"
#define INDENTATION_SURFACE_R0_STRING		"indentationSurfaceR0"
#define INDENTATION_SURFACE_N_STRING		"indentationSurfaceN"
#define INDENTATION_TIP_RADIUS_STRING		"indentationTipR"
#define INDENTATION_KS_STRING				"indentationKs"
#define INDENTATION_DELTAX_STRING			"indentationDeltaX"
#define INDENTATION_MOVE_SURFACE			"indentationMoveSurface"
#define INDENTATION_FIX_TRANSVERSAL			"indentationFixTrans"
#define INDENTATION_SIGMA_STRING			"indentationSigma"
#define INDENTATION_EL_STRING				"indentationEl"
#define INDENTATION_TIP_SIGMA_STRING		"indentationTipSigma"
#define INDENTATION_TIP_EL_STRING			"indentationTipEl"
#define INDENTATION_SURF_SIGMA_STRING		"indentationSurfSigma"
#define INDENTATION_SURF_EL_STRING			"indentationSurfEl"
#define INDENTATION_TIP_A					"indentationTipA"
#define INDENTATION_TIP_B					"indentationTipB"
#define INDENTATION_SURFACE_A				"indentationSurfA"
#define INDENTATION_SURFACE_B				"indentationSurfB"

#define INDENTATION_TIP_ZETA_STRING			"indentationTipZeta"
#define INDENTATION_OUTPUT_FILENAME_STRING	"indentationOutput"
#define INDENTATION_RETRACTION_STEP_STRING	"indentationRetractionStep"
#define INDENTATION_OUTPUT_FREQ_STRING		"indentationOutputFreq"

#define INDENTATION_DISCRETE_SURF_STRING	"indentationDiscreteSurf"
#define INDENTATION_PAIRS_CUTOFF_STRING		"indentationPairsCutoff"

#define INDENTATION_SHOW_TIP_SURFACE_STRING "indentationShowTipSurf"
#define INDENTATION_VMD_CONNECT_SCRIPT		"indentationSurfConnectFile"
#define INDENTATION_MICA_SIZE_STRING		"indentationSurfaceSize"
#define INDENTATION_MICA_STEP_STRING		"indentationSurfaceStep"
#define INDENTATION_CANTILEVER_LENGTH		"indentationCantLength"

#define DEFAULT_INDENTATION_FIX_TRANSVERSAL		0
#define DEFAULT_INDENTATION_FILENAME			"indentation.<name>_<author><run>_<stage>.dat"
#define DEFAULT_INDENTATION_TIP_A				0.0f
#define DEFAULT_INDENTATION_TIP_B				1.0f
#define DEFAULT_INDENTATION_SURFACE_A			0.0f
#define DEFAULT_INDENTATION_SURFACE_B			1.0f
#define DEFAULT_INDENTATION_SHOW_TIP_SURFACE 	0
#define DEFAULT_INDENTATION_MICA_SIZE			51
#define DEFAULT_INDENTATION_MICA_STEP			10
#define DEFAULT_INDENTATION_CANTILEVER_LENGTH	500.0f
#define DEFAULT_INDENTATION_OUTPUT_FREQ			1000
#define DEFAULT_INDENTATION_PAIRS_CUTOFF		40.0f

typedef struct {
	float3 tipCoord;
	float3 chipCoord;
	float3 chipCoord0;
	float3 direction;
	float V;
	float dx;
	int moveSurface;
	int fixTransversal;
	float tipRadius;
	float cantileverKs;
	float tipa6;
	float tipel;
	float surfa6;
	float surfel;
	float tipAprime;
	float tipBprime;
	float surfAprime;
	float surfBprime;
	float tipZeta;
	float3 micaN;
	float3 micaR0;
	float3 micaR;
	float4* h_tipForces;
	float4* d_tipForces;
	float4 tipForce;
	FILE* out;
	long int retractionStep;
	int showTipSurface;
	int surfaceSize;
	int surfaceBeadsCount;
	float surfaceStep;
	float3* surfacePointsR0;
	float4* h_surfacePointsCoord;
	float4* d_surfacePointsCoord;

	float pairsCutoff2;
	int* h_micaListCounts;
	int* d_micaListCounts;
	int* h_micaList;
	int* d_micaList;

	int outputFreq;
	float3 cantileverVector;
	float4 fav;
	float3 tipCoordAv;
	float3 chipCoordAv;
	float kDeltaXAv;

} Indentation;

void createIndentationPotential();
void initIndentation();
inline void computeIndentation();
inline void computeIndentationDiscreteSurface();
inline void computeIndentationEnergy();
inline void updateTipPosition();
void updateAdditionalAminos();

