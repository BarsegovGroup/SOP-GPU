/*
 * topio.h
 *
 *  Created on: May 24, 2009
 *      Author: zhmurov
 */

#pragma once

#include "pdbio.h"

struct PDBAtom;

struct CovalentBond{
	int i;
	int j;
	float r0;
};

struct NativeContact{
	int i;
	int j;
	float r0;
	float eh;
};

struct PossiblePair{
	int i;
	int j;
};

struct TOPPair{
	int i;
	int j;
	int func;
	float c0;
	float c1;
	float c2;
	float c3;

    void save(FILE* topFile) const;
    void read(FILE* topFile);
};

struct SOP {
	int aminoCount;
	int bondCount;
	int nativeCount;
	int pairCount;
	PDBAtom* aminos;
	CovalentBond* bonds;
	NativeContact* natives;
	PossiblePair* pairs;

	int additionalAminosCount;
	PDBAtom* additionalAminos;

    void load(const char* filename);
    void save(const char* filename) const;
};

