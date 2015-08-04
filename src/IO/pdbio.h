#pragma once

#include <vector>
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

    PDBAtom(const char* line) { this->parse(line); }
    PDBAtom() { }
    void fprint(FILE* file, int custom_atomid = -1, bool atomid_in_hex = false) const;
    void parse(const char* line);
    void fromTOP(FILE* topFile);
};

struct PDBSSBond {
	int resid1, resid2;
	char chain1, chain2;

    PDBSSBond(const char *line) { this->parse(line); }
    PDBSSBond() { }
    void fprint(FILE* file, int bondid) const;
    void parse(const char* line);
};

struct PDBConnect {
	int* connectCount;
	int* connectMap;
};

struct PDB {
    std::vector<PDBAtom> atoms;
    std::vector<PDBSSBond> ssbonds;
	PDBConnect connections;

    void read(const char *filename);
    void readXYZ(const char *filename);
    void write(const char *filename, bool printConnections = false) const;
    void fromSOP(const SOP& sop);
    void toSOP(SOP& sop) const;
    bool is_xyz;
};

void savePDB(const char* filename, const SOP& sop);
void readCoord(const char* filename, SOP& sop);

