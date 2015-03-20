/*
 * indentation.cuh
 *
 *  Created on: Apr 9, 2010
 *      Author: zhmurov
 */

#pragma once

#include "../Util/parameters.h"

PARAMETER(indentation, bool, false, "true/false", "...")

PARAMETER_MANDATORY(indentationChip, float3, "?", "...")
PARAMETER(indentationTip, float3, parameters::indentationChip, "?", "...")
PARAMETER_MANDATORY(indentationDirection, float3, "?", "...")
PARAMETER_MANDATORY(indentationSurfaceR0, float3, "?", "...")
PARAMETER_MANDATORY(indentationSurfaceN, float3, "?", "...")
PARAMETER_MANDATORY(indentationTipR, float, "?", "...")
PARAMETER_MANDATORY(indentationKs, float, "?", "...")
PARAMETER_MANDATORY(indentationDeltaX, float, "?", "...")
PARAMETER(indentationMoveSurface, bool, false, "true/false", "...")
PARAMETER(indentationFixTrans, bool, false, "true/false", "...")
PARAMETER(indentationSigma, float, 1.0f, "?", "...")
PARAMETER(indentationTipSigma, float, parameters::indentationSigma, "?", "...")
PARAMETER(indentationSurfSigma, float, parameters::indentationSigma, "?", "...")
PARAMETER(indentationEl, float, 1.0f, "?", "...")
PARAMETER(indentationTipEl, float, parameters::indentationEl, "?", "...")
PARAMETER(indentationSurfEl, float, parameters::indentationEl, "?", "...")
PARAMETER(indentationTipA, float, 0.0f, "?", "...")
PARAMETER(indentationTipB, float, 1.0f, "?", "...")
PARAMETER(indentationSurfA, float, 0.0f, "?", "...")
PARAMETER(indentationSurfB, float, 1.0f, "?", "...")

PARAMETER(indentationTipZeta, float, 5000.0f, "?", "...")
PARAMETER(indentationOutput, std::string, "indentation.<name>_<run>.dat", "path", "...")
PARAMETER(indentationRetractionStep, long, -1, "?", "...")
PARAMETER(indentationOutputFreq, int, 1000, "steps", "...")

PARAMETER(indentationDiscreteSurf, bool, false, "true/false", "...")
PARAMETER(indentationPairsCutoff, float, 40.0f, "?", "...")

PARAMETER(indentationShowTipSurf, bool, false, "true/false", "...")
PARAMETER(indentationSurfConnectFile, std::string, "connect_mica.vmd", "path", "...");
PARAMETER(indentationSurfaceSize, int, 51, "?", "...")
PARAMETER(indentationSurfaceStep, float, 1.4, "?", "...")
PARAMETER(indentationCantLength, float, 500.0f, "?", "...")

void createIndentationPotential();

class IndentationPotential : public SOPPotential{
public:
    IndentationPotential();
    virtual ~IndentationPotential() { }

    virtual void compute();
    bool discreteSurf;
};

class IndentationTipUpdater : public SOPUpdater{
public:
    IndentationTipUpdater(IndentationPotential *indentation);
    virtual void update();
private:
    IndentationPotential *indentation;
	int outputFreq;
	float4 fav;
	float3 tipCoordAv;
	float3 chipCoordAv;
	float kDeltaXAv;
	FILE* outputFile; // .dat output file
	long int retractionStep; // Step at which the direction of the cantilever base movement will be reversed
};

class IndentationAminoUpdater : public SOPUpdater{
public:
    IndentationAminoUpdater();
    virtual void update();
private:
    float3 cantileverVector;
};

