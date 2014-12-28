/*
 * sop.c
 *
 *  Created on: May 26, 2009
 *      Author: zhmurov
 */

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "IO/pdbio.h"
#include "IO/topio.h"

SOP sop;
PDB pdbdata;
SOP tandem;

int* pdbRef;
int* aminoRef;
int* firstAtomInResid;
int* lastAtomInResid;

float getDistanceAtoms(int i1, int i2);
float getDistanceBeads(int i, int j);
float getOccupancy(int i);
float getBeta(int i);

float getX(int i);
float getY(int i);
float getZ(int i);

int checkCovalent(int i, int j);
int checkNative(int i, int j);
int checkPairs(int i, int j);
int checkPossiblePairs(int i, int j);

float getR0(int i, int j);
float getEh(int i, int j);

char pdb_filename[FILENAME_LENGTH];
char top_filename[FILENAME_LENGTH];
char coord_filename[FILENAME_LENGTH];

float a;
float eh;
float R_limit_bond;
float SC_limit_bond;
int covalentLJ;
float pairs_threshold;
float covalent_cutoff;
float SS_cutoff;

int createTandem;
int linkerLength;
int monomerCount;
int tandemDirection = 1;
float tandemVectorX;
float tandemVectorY;
float tandemVectorZ;

int firstResid;
int lastResid;

typedef struct {
	int resid1;
	int resid2;
	char residName1[5];
	char residName2[5];
} CovalentLinker;

int covLinkersCount;
CovalentLinker* covLinkers;
float covLinkerCutoff;

int checkCovalent(int i, int j){
	int i1 = pdbRef[i];
	int i2 = pdbRef[j];
	// For backbone
	if((pdbdata.atoms[i1].resid == pdbdata.atoms[i2].resid + 1 && pdbdata.atoms[i1].chain == pdbdata.atoms[i2].chain)
			|| (pdbdata.atoms[i2].resid == pdbdata.atoms[i1].resid + 1 && pdbdata.atoms[i2].chain == pdbdata.atoms[i1].chain)){
		if(getDistanceAtoms(i1, i2) < covalent_cutoff && abs(i-j) == 1){
			return 1;
		}
	}
	// For SS-bonds
	int k;
	for(k = 0; k < pdbdata.ssCount; k++){
		if((pdbdata.atoms[i1].resid == pdbdata.ssbonds[k].resid1) && (pdbdata.atoms[i1].chain == pdbdata.ssbonds[k].chain1) &&
				(pdbdata.atoms[i2].resid == pdbdata.ssbonds[k].resid2) && (pdbdata.atoms[i2].chain == pdbdata.ssbonds[k].chain2)){
			if(getDistanceAtoms(i1, i2) < SS_cutoff){
				return 1;
			}
		}
		if((pdbdata.atoms[i1].resid == pdbdata.ssbonds[k].resid2) && (pdbdata.atoms[i1].chain == pdbdata.ssbonds[k].chain2) &&
						(pdbdata.atoms[i2].resid == pdbdata.ssbonds[k].resid1) && (pdbdata.atoms[i2].chain == pdbdata.ssbonds[k].chain1)){
			if(getDistanceAtoms(i1, i2) < SS_cutoff){
				return 1;
			}
		}
	}
	// For Linkers
	for(k = 0; k < covLinkersCount; k++){
		if((pdbdata.atoms[i1].resid == covLinkers[k].resid1) &&	(strncmp(pdbdata.atoms[i1].resName, covLinkers[k].residName1, 3) == 0) &&
				(pdbdata.atoms[i2].resid == covLinkers[k].resid2) && (strncmp(pdbdata.atoms[i2].resName, covLinkers[k].residName2, 3) == 0)){
			if(getDistanceAtoms(i1, i2) < covLinkerCutoff){
				//printf("\n1. Linking residues %d and %d. R = %f\n", i, j, getDistanceAtoms(i1, i2));
				return 1;
			} else {
				return 0;
			}
		}
		if((pdbdata.atoms[i1].resid == covLinkers[k].resid2) &&	(strncmp(pdbdata.atoms[i1].resName, covLinkers[k].residName2, 3) == 0) &&
				(pdbdata.atoms[i2].resid == covLinkers[k].resid1) && (strncmp(pdbdata.atoms[i2].resName, covLinkers[k].residName1, 3) == 0)){
			if(getDistanceAtoms(i1, i2) < covLinkerCutoff){
				//printf("\n2. Linking residues %d and %d. R = %f\n", i, j, getDistanceAtoms(i1, i2));
				return 1;
			} else {
				return 0;
			}
		}
	}
	return 0;
}

int checkNative(int i, int j){
	if(!checkCovalent(i, j)){
		int i1 = pdbRef[i];
		int i2 = pdbRef[j];
		if(!(abs(pdbdata.atoms[i1].resid - pdbdata.atoms[i2].resid) <= 2
				&& pdbdata.atoms[i1].chain == pdbdata.atoms[i2].chain)){
			if(strcmp(pdbdata.atoms[i1].name, "CA") == 0 && strcmp(pdbdata.atoms[i2].name, "CA") == 0){
				if(getDistanceAtoms(i1, i2) < R_limit_bond){
	#ifdef DEBUG1
				printf("CA-CA (%f):\n", getDistanceAtoms(i1, i2));
				printAtom(pdbdata.atoms[i1]);
				printAtom(pdbdata.atoms[i2]);
	#endif
					return 1;
				}
			}
		}
		if(SC_limit_bond > 0.0f){
			if(!(abs(pdbdata.atoms[i1].resid - pdbdata.atoms[i2].resid) <= 2
								&& pdbdata.atoms[i1].chain == pdbdata.atoms[i2].chain)){
				for(i1 = firstAtomInResid[i]; i1 <= lastAtomInResid[i]; i1++){
					for(i2 = firstAtomInResid[j]; i2 <= lastAtomInResid[j]; i2++){
						if(strcmp(pdbdata.atoms[i1].name, "CA") != 0 && strcmp(pdbdata.atoms[i2].name, "CA") != 0
								&& strcmp(pdbdata.atoms[i1].name, "N") != 0 && strcmp(pdbdata.atoms[i2].name, "N") != 0
								&& strcmp(pdbdata.atoms[i1].name, "C") != 0 && strcmp(pdbdata.atoms[i2].name, "C") != 0
								&& strcmp(pdbdata.atoms[i1].name, "O") != 0 && strcmp(pdbdata.atoms[i2].name, "O") != 0 ){
							if(getDistanceAtoms(i1, i2) < SC_limit_bond){
								return 1;
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

int checkPairs(int i, int j){
	int i1 = pdbRef[i];
	int i2 = pdbRef[j];
	float dr = getDistanceAtoms(i1 ,i2);
	if(!checkPossiblePairs(i, j)){
		return 0;
	}
	if((dr < pairs_threshold)){
		return 1;
	}

	return 0;

}

int checkPossiblePairs(int i, int j){
	/*int i1 = pdbRef[i];
	int i2 = pdbRef[j];
	float dr = getDistanceAtoms(i1, i2);*/
	if(pairs_threshold > 0.0f){
		int i1 = pdbRef[i];
		int i2 = pdbRef[j];
		float dr = getDistanceAtoms(i1, i2);
		if(checkCovalent(i,j)){
			if(covalentLJ){
				if(dr < pairs_threshold){
					return 1;
				} else {
					return 0;
				}
			} else {
				return 0;
			}
		}
		if(checkNative(i,j)){
			return 0;
		}
		if(i != j){
		//if((dr < R_limit_bond && abs(i-j) <= 2 && i != j) || dr > R_limit_bond){
			if(dr < pairs_threshold){
				return 1;
			} else {
				return 0;
			}
		}
	}
	return 0;
}

float getR0(int i, int j){
	float dr = getDistanceBeads(i, j);
	if(dr < 1.0 && i != j){
		printf("WARNING: Residues %d and %d are closer than 1.0A!\n", i, j);
	}
	//if(dr < R_limit_bond && i != j+1 && j != i+1 && i != j){
	//if((abs(i-j) >= 3 && checkNative(i, j)) || (abs(i-j) <= 2 && dr <= R_limit_bond)){
	if(checkNative(i, j) || checkCovalent(i,j)){
		return dr;
	}
	return a;
}

float getEh(int i, int j){
	if(eh == -1.0){
		return sqrtf(getOccupancy(i)*getOccupancy(j));
	} else if(eh == -2.0){
		return sqrtf(getBeta(i)*getBeta(j));
	} else {
		return eh;
	}
}

float getDistanceAtoms(int i1, int i2){
	float dr;
	float x = pdbdata.atoms[i1].x - pdbdata.atoms[i2].x;
	float y = pdbdata.atoms[i1].y - pdbdata.atoms[i2].y;
	float z = pdbdata.atoms[i1].z - pdbdata.atoms[i2].z;
	dr = sqrtf(x*x+y*y+z*z);
	return dr;
}

float getDistanceBeads(int i, int j){
	int i1 = pdbRef[i];
	int i2 = pdbRef[j];
	return getDistanceAtoms(i1, i2);
}


float getOccupancy(int i){
	return pdbdata.atoms[pdbRef[i]].occupancy;
}

float getBeta(int i){
	return pdbdata.atoms[pdbRef[i]].beta;
}

float getX(int i){
	return sop.aminos[i].x;
}

float getY(int i){
	return sop.aminos[i].y;
}

float getZ(int i){
	return sop.aminos[i].z;
}

