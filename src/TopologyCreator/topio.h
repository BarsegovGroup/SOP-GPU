/*
 * topio.h
 *
 *  Created on: May 24, 2009
 *      Author: zhmurov
 */

#pragma once

#define TOP_ATOMTYPE_LENGTH 16
#define TOP_ATOMNAME_LENGTH 8
#define TOP_RESNAME_LENGTH 8

#define TOP_SECTION_ATOMS		0
#define TOP_SECTION_PAIRS		1
#define TOP_SECTION_ANGLES		2
#define TOP_SECTION_DIHEDRALS	3

typedef struct {
	int id;
	char type[TOP_ATOMTYPE_LENGTH];
	int resid;
	char resName[TOP_RESNAME_LENGTH];
	char name[TOP_ATOMNAME_LENGTH];
	char chain;
	float charge;
	float mass;
} TOPAtom;

typedef struct {
	int i;
	int j;
	int func;
	float c0;
	float c1;
	float c2;
	float c3;
} TOPPair;

typedef struct {
	int i;
	int j;
	int k;
	int func;
	float c0;
	float c1;
	float c2;
	float c3;
} TOPAngle;

typedef struct {
	int i;
	int j;
	int k;
	int l;
	int func;
	float c0;
	float c1;
	float c2;
	float c3;
} TOPDihedral;

typedef struct {
	int atomCount;
	int bondCount;
	int pairsCount;
	int nativesCount;
	int angleCount;
	int dihedralCount;
	TOPAtom* atoms;
	TOPPair* bonds;
	TOPPair* natives;
	TOPPair* pairs;
	TOPAngle* angles;
	TOPDihedral* dihedrals;
} TOPData;

extern int readTOP(const char* filename, TOPData* topData);
extern void writeTOP(const char* filename, TOPData* topData);
