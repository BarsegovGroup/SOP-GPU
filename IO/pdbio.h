#pragma once

#include <stdio.h>
#include "topio.h"

struct SOP;

/*
 * pdbio.h
 *
 *  Created on: Jul 5, 2009
 *      Author: zhmurov
 */

/*
 * Structures
 */

struct PDBAtom {
    int    id;
    char   name[5], chain, resName[4], altLoc;
    int    resid;
    double x, y, z;
    double occupancy;
    double beta;

    void print() const;
    void fprint(FILE* file) const;
    void parse(const char* line);
    void fromTOP(FILE* topFile);
};

struct PDBSSBond {
	int resid1, resid2;
	char chain1, chain2;
    
    void parse(const char* line);
};

struct PDBConnect {
	int* connectCount;
	int* connectMap;
};

struct PDB {
	int atomCount;
	PDBAtom* atoms;
	int ssCount;
	PDBSSBond* ssbonds;
	PDBConnect connections;

    void read(const char *filename);
    void write(const char *filename, bool printConnections = false) const;
    void fromSOP(const SOP& sop);
    void toSOP(SOP& sop) const;
};

void savePDB(const char* filename, const SOP& sop);
void readCoord(const char* filename, SOP& sop);

