/*
 * pulling.cuh
 *
 *  Created on: May 25, 2010
 *      Author: zhmurov
 */

#pragma once

#define PULLING_ON_STRING			"pulling"

#define PULLING_DELTAX_STRING		"deltax"
#define PULLING_KS_STRING			"k_trans"
#define PULLING_FCONST_STRING		"fconst"
#define PULLING_FIXED_COUNT_STRING	"fixed_beads"
#define PULLING_FIXED_STRING		"fixed"
#define PULLING_PULLED_COUNT_STRING	"pulled_beads"
#define PULLING_PULLED_STRING		"pulled"
#define PULLING_DIRECTION_STRING	"pullDirection"
#define PULLING_VECTOR_STRING		"pullVector"
#define PULLING_FIXED_END_STRING	"fixedEnd"
#define PULLING_PULLED_END_STRING	"pulledEnd"
#define PULLING_FILENAME			"pullOutput"
#define PULLING_FREQ				"pullFreq"

#define PULLING_DIRECTION_ENDTOEND_STRING	"endToEnd"
#define PULLING_DIRECTION_VECTOR_STRING		"vector"

#define DEFAULT_PULLING_KS			0.05f
#define DEFAULT_PULLING_DIRECTION	PULLING_DIRECTION_ENDTOEND_STRING
#define DEFAULT_PULLING_FILENAME	"pull.<name>_<author><run>_<stage>.dat"

struct Pulling {
	float Ks;
	float deltax;
	float fconst;
	int fixedEnd;
	int pulledEnd;
	int fixedCount;
	int pulledCount;
	int* fixed;
	int* pulled;
	float3* pullVector;
	float3* extForce;
	float4* h_extForces;
	float4* d_extForces;
	float3* cantCoord0;
	float3* cantCoord;
};

void createPullingPotential();
void initPulling();
inline void computePulling();
inline void computePullingEnergy();
inline void updatePulling();

