/*
 * topio.h
 *
 *  Created on: May 24, 2009
 *      Author: zhmurov
 */

#pragma once

#include <vector>
#include "pdbio.h"

struct PDBAtom;

struct TOPPair{
	int i;
	int j;
	int func;
	float c0;
	float c1;
	float c2;
	float c3;

    void save(FILE* topFile) const;
    void parse(char* buffer);
    void read(FILE* topFile);
};

struct CovalentBond{
    CovalentBond(const TOPPair& pair) : i(pair.i), j(pair.j), r0(pair.c0) {  }
    CovalentBond() {  }
	int i;
	int j;
	float r0;
};

struct NativeContact{
    NativeContact(const TOPPair& pair) : i(pair.i), j(pair.j), r0(pair.c0), eh(pair.c1) { }
    NativeContact() { }
	int i;
	int j;
	float r0;
	float eh;
};

struct PossiblePair{
    PossiblePair(const TOPPair& pair) : i(pair.i), j(pair.j) { }
    PossiblePair() { }
	int i;
	int j;
};

struct SOP {
    std::vector<PDBAtom> aminos;
    std::vector<CovalentBond> bonds;
    std::vector<NativeContact> natives;
    std::vector<PossiblePair> pairs;

    std::vector<PDBAtom> additionalAminos;

    void load(const char* filename);
    void save(const char* filename) const;
};

