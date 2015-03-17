/*
 * aatocg.h
 *
 *  Created on: Nov 28, 2012
 *      Author: kias
 */

#ifndef AATOCG_H_
#define AATOCG_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

#define NAME_LENGTH 5
#define TYPE_LENGTH 5
#define RESNAME_LENGTH 5

typedef struct {
	char name[NAME_LENGTH];
} Atom;

typedef struct {
	char resname[NAME_LENGTH];
	char name[NAME_LENGTH];
	char type[TYPE_LENGTH];
	double mass;
	double charge;
	std::vector<Atom> atomsRepresents;
	std::vector<Atom> atomsCM;
	std::vector<Atom> connectedTo;
} Bead;

typedef struct {
	char resname[NAME_LENGTH];
	std::vector<Bead> beads;
} Residue;

typedef struct {
	char name;
	double mass;
} Mass;

typedef struct {
	std::vector<Residue> residues;
	std::vector<Mass> masses;
} CGConfig;

extern void readCGConfiguration(const char* filename, CGConfig *config);

#endif /* AATOCG_H_ */
