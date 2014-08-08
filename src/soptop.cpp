/*
 * gsop-top.c
 *
 *  Created on: Oct 21, 2010
 *      Author: zhmurov
 */

#include "sop.cpp"
#include "IO/configreader.h"
#include "IO/topio.h"
#include "structure.h"
#include "IO/pdbio.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define version "trunk"
#define FILENAME_LENGTH 100

/*float R_limit_bond;
float SC_limit_bond;
float a;
int covalentLJ;
float pairs_cutoff;
float pairs_threshold;


int createTandem;
int linkerLength;
int monomerCount;
int tandemDirection = 1;
float tandemVectorX;
float tandemVectorY;
float tandemVectorZ;
int firstResid;
int lastResid;*/

SOP sop;
PDB pdbdata;
SOP tandem;

//char pdb_filename[FILENAME_LENGTH];
char top_filename[FILENAME_LENGTH];
char coord_filename[FILENAME_LENGTH];

void createModel();

int main(int argc, char *argv[]){

	printf("==========================\n");
	printf("gSOP Topology creator version %s\n", version);
	printf("==========================\n");

	if(argc < 1){
		printf("ERROR: Configuration file should be specified.\n");
		exit(-1);
	}

	parseParametersFile(argv[1]); //, argc, argv);

	createModel();
}


void createModel(){

	char ehString[30];
	getParameter(ehString, "eh", "1.5", 1);
	if(strcmp(ehString, "O") == 0){
		printf("Taking eh's values from occupancy column.\n");
		eh = -1.0;
	} else if(strcmp(ehString, "B") == 0){
		printf("Taking eh's values from beta column.\n");
		eh = -2.0;
	} else {
		eh = atof(ehString);
		if(eh == 0){
			printf("WARNING: Value of eh in a configuration file is not valid. Should be:\n "
					" - Positive integer, or\n"
					" - 'O' to take values from occupancy column, or\n"
					" - 'B' to take from beta column.\n"
					"eh is set to 0.\n");
		}
	}
	createTandem = getYesNoParameter("createTandem", 0, 1);
	if(createTandem){
		linkerLength = getIntegerParameter("linkerLength", 0, 0);
		monomerCount = getIntegerParameter("monomerCount", 0, 0);
		char tandemDirectionString[100];
		getMaskedParameter(tandemDirectionString, "tandemDirection", "endToEnd", 1);
		if(strcmp(tandemDirectionString, "endToEnd") == 0){
			tandemDirection = 1;
			printf("Direction of tandem linker is a monomer end-to-end direction.\n");
		} else
		if(strcmp(tandemDirectionString, "vector") == 0){
			tandemDirection = 2;
			tandemVectorX = getFloatParameter("tandemVectorX", 0, 0);
			tandemVectorY = getFloatParameter("tandemVectorY", 0, 0);
			tandemVectorZ = getFloatParameter("tandemVectorZ", 0, 0);
			printf("Tandem linker direction is (%3.2f, %3.2f, %3.2f) vector.\n", tandemVectorX, tandemVectorY, tandemVectorZ);
		} else {
			tandemDirection = 1;
			printf("WARNING: Direction of tandem linker is not specified or wrong. Assuming monomer end-to-end direction.\n");
		}
		firstResid = getIntegerParameter("fixedEnd", 0, 0);
		lastResid = getIntegerParameter("pulledEnd", 0, 0);
	}
	getMaskedParameter(pdb_filename, "structure", "", 0);
	getMaskedParameter(top_filename, "topology", "", 0);
	getMaskedParameter(coord_filename, "coordinates", "", 0);


	R_limit_bond = getFloatParameter("R_limit_bond", 8.0f, 1);
	SC_limit_bond = getFloatParameter("SC_limit_bond", 0.0f, 1);
	covalentLJ = getYesNoParameter("covalentRepulsion", 0, 1);
	a = getFloatParameter("a", 3.8f, 1);
	pairs_threshold = getFloatParameter("pairs_threshold", 1000.0f, 1);
	covalent_cutoff = getFloatParameter("covalent_cutoff", 10.0f, 1);
	SS_cutoff = getFloatParameter("SS_cutoff", 10.0f, 1);

	int i, j, k, i1, i2;
	covLinkersCount = getIntegerParameter("covLinkersCount", 0, 1);
	if(covLinkersCount > 0){
		covLinkers = (CovalentLinker*)calloc(covLinkersCount, sizeof(CovalentLinker));
		covLinkerCutoff = getFloatParameter("covLinkersCutoff", 5.0f, 1);
		char covLinker[50];
		char covLinkerName[50];
		char* pch;
		for(i = 0; i < covLinkersCount; i++){
			sprintf(covLinkerName, "covLinker%d", i+1);
			getMaskedParameter(covLinker, covLinkerName, "", 0);
			pch = strtok(covLinker, "-\t");
			strncpy(covLinkers[i].residName1, pch, 3);
			covLinkers[i].resid1 = atoi(&pch[3]);
			pch = strtok(NULL, "\n\r\t ");
			strncpy(covLinkers[i].residName2, pch, 3);
			covLinkers[i].resid2 = atoi(&pch[3]);
			printf("Covalent linker between residues %s%d and %s%d will be added.\n",
					covLinkers[i].residName1, covLinkers[i].resid1,
					covLinkers[i].residName2, covLinkers[i].resid2);
		}
	}

	pdbdata.read(pdb_filename);
	printf("Creating SOP-model...\n");

	sop.aminoCount = 0;
	//Counting number of Calphas
	for(i = 0; i < pdbdata.atomCount; i++){
		if(strcmp(pdbdata.atoms[i].name, "CA") == 0){
			sop.aminoCount ++;
		}
	}
	printf("Found %d residues.\n", sop.aminoCount);

	pdbRef = (int*)malloc(sop.aminoCount*sizeof(int));
	firstAtomInResid = (int*)malloc(sop.aminoCount*sizeof(int));
	lastAtomInResid = (int*)malloc(sop.aminoCount*sizeof(int));
	aminoRef = (int*)malloc(pdbdata.atomCount*sizeof(int));
	sop.aminos = (PDBAtom*)malloc(sop.aminoCount*sizeof(PDBAtom));
	/*nativeContacts = (char**)malloc(sop.aminoCount*sop.aminoCount*sizeof(char));
	for(i = 0; i < sop.aminoCount; i++){
		nativeContacts[i] = (char*)malloc(sop.aminoCount*sizeof(char));
	}*/

	int index = 0;
	for(i = 0; i < pdbdata.atomCount; i++){
		if(strcmp(pdbdata.atoms[i].name, "CA") == 0){
			pdbRef[index] = i;
			sop.aminos[index] = pdbdata.atoms[i];
			index ++;
		}

	}

	index = -1;
	int lastAminoResid = -1;// pdbdata.atoms[0].resid;
	char lastAminoChain = 0;//pdbdata.atoms[0].chain;
	for(i = 0; i < pdbdata.atomCount; i++){
		if(lastAminoResid != pdbdata.atoms[i].resid || lastAminoChain != pdbdata.atoms[i].chain){
			lastAminoResid = pdbdata.atoms[i].resid;
			lastAminoChain = pdbdata.atoms[i].chain;
			index ++;
		}
		aminoRef[i] = index;
#ifdef DEBUG1
	printf("%d - %d\n", i, index);
#endif
	}

	index = 0;
	firstAtomInResid[0] = 0;
	lastAtomInResid[sop.aminoCount - 1] = pdbdata.atomCount - 1;
	for(i = 0; i < sop.aminoCount - 1; i++){
		while(aminoRef[index] == aminoRef[index+1]){
			index ++;
		}
		index ++;
		lastAtomInResid[i] = index;
		firstAtomInResid[i+1] = index+1;
	}

	/*for(i = 0; i < sop.aminoCount; i++){
		printf("Residue %d : %d - %d\n", i, firstAtomInResid[i], lastAtomInResid[i]);
	}*/


/*	for(i = 0; i < sop.aminoCount; i++){
		for(j = 0; j < sop.aminoCount; j++){
			nativeContacts[i][j] = 0;
		}
	}*/

/*	int sc_count = 0, ca_count = 0;
	for(i1 = 0; i1 < pdbdata.atomCount; i1++){
		for(i2 = 0; i2 < pdbdata.atomCount; i2++){
			i = aminoRef[i1];
			j = aminoRef[i2];
			if(!(abs(pdbdata.atoms[i1].resid - pdbdata.atoms[i2].resid) <= 2
					&& pdbdata.atoms[i1].chain == pdbdata.atoms[i2].chain)){
				if(strcmp(pdbdata.atoms[i1].name, "CA") == 0 && strcmp(pdbdata.atoms[i2].name, "CA") == 0){
					if(getDistanceAtoms(i1, i2) < R_limit_bond && nativeContacts[i][j] == 0){
#ifdef DEBUG1
					printf("CA-CA (%f):\n", getDistanceAtoms(i1, i2));
					printAtom(pdbdata.atoms[i1]);
					printAtom(pdbdata.atoms[i2]);
#endif
						nativeContacts[i][j] = 1;
						ca_count ++;
					}
				}
			}
		}
	}
	if(SC_limit_bond > 0.0f){
		for(i1 = 0; i1 < pdbdata.atomCount; i1++){
				for(i2 = 0; i2 < pdbdata.atomCount; i2++){
					i = aminoRef[i1];
					j = aminoRef[i2];
					if(!(abs(pdbdata.atoms[i1].resid - pdbdata.atoms[i2].resid) <= 2
										&& pdbdata.atoms[i1].chain == pdbdata.atoms[i2].chain)){
						if(strcmp(pdbdata.atoms[i1].name, "CA") != 0 && strcmp(pdbdata.atoms[i2].name, "CA") != 0
								&& strcmp(pdbdata.atoms[i1].name, "N") != 0 && strcmp(pdbdata.atoms[i2].name, "N") != 0
								&& strcmp(pdbdata.atoms[i1].name, "C") != 0 && strcmp(pdbdata.atoms[i2].name, "C") != 0
								&& strcmp(pdbdata.atoms[i1].name, "O") != 0 && strcmp(pdbdata.atoms[i2].name, "O") != 0 ){
							if(getDistanceAtoms(i1, i2) < SC_limit_bond && nativeContacts[i][j] == 0){
	#ifdef DEBUG1
							printf("SC-SC (%f):\n", getDistanceAtoms(i1, i2));
							printAtom(pdbdata.atoms[i1]);
							printAtom(pdbdata.atoms[i2]);
	#endif
								nativeContacts[i][j] = 2;
								sc_count ++;
							}
						}
					}
				}
		}
	}
	printf("CA native contacts: %d. \nSC native contacts: %d. \nTotal:%d(%d)\n", ca_count, sc_count, ca_count + sc_count, (ca_count + sc_count)/2);*/


	sop.bondCount = 0;
	sop.nativeCount = 0;
	sop.pairCount = 0;
	//long int computed = 0;
	for(i = 0; i < sop.aminoCount; i++){
		for(j = i + 1; j < sop.aminoCount; j++){
			if(checkCovalent(i, j)){
				sop.bondCount ++;
			}
			if(checkNative(i, j)){
				sop.nativeCount ++;
			}
			if(checkPossiblePairs(i, j)){
				sop.pairCount ++;
			}
			//computed++;
		}
		printf("Amino: %d. Bonds: %d. Native contacts: %d. Pairs: %d. (%4.2f%% completed)\r",
				i,
				sop.bondCount,
				sop.nativeCount,
				sop.pairCount,
				((float)i+1)*100.0f/((float)(sop.aminoCount)));
	}
	printf("\n");

	printf("Total per monomer:\n%d amino acids\n%d covalent bonds\n%d native contacts\n%d pairs\n=======\n",
			sop.aminoCount, sop.bondCount, sop.nativeCount, sop.pairCount);


	sop.aminos = (PDBAtom*)calloc(sop.aminoCount, sizeof(PDBAtom));
	sop.bonds = (CovalentBond*)calloc(sop.bondCount, sizeof(CovalentBond));
	sop.natives = (NativeContact*)calloc(sop.nativeCount, sizeof(NativeContact));
	sop.pairs = (PossiblePair*)calloc(sop.pairCount, sizeof(PossiblePair));



	index = 0;
	for(i = 0; i < pdbdata.atomCount; i++){
		if(strcmp(pdbdata.atoms[i].name, "CA") == 0){
			memcpy(&sop.aminos[index], &pdbdata.atoms[i], sizeof(PDBAtom));
			index ++;
		}
	}
	sop.bondCount = 0;
	sop.nativeCount = 0;
	sop.pairCount = 0;

	for(i = 0; i < sop.aminoCount; i++){
		for(j = i + 1; j < sop.aminoCount; j++){
			if(checkCovalent(i, j)){
				sop.bonds[sop.bondCount].i = i;
				sop.bonds[sop.bondCount].j = j;
				sop.bonds[sop.bondCount].r0 = getR0(i, j);
				sop.bondCount ++;
			}
			if(checkNative(i, j)){
				sop.natives[sop.nativeCount].i = i;
				sop.natives[sop.nativeCount].j = j;
				sop.natives[sop.nativeCount].r0 = getR0(i, j);
				sop.natives[sop.nativeCount].eh = getEh(i, j);
				sop.nativeCount ++;
			}
			if(checkPossiblePairs(i, j)){
				sop.pairs[sop.pairCount].i = i;
				sop.pairs[sop.pairCount].j = j;
				sop.pairCount ++;
			}
		}
		printf("Creating model ... %5.2f%% completed\r", ((float)i+1)*100.0f/((float)(sop.aminoCount)));
	}
	printf("\n");

	if(createTandem){


		printf("Building a tandem of %d\n", monomerCount);

		tandem.aminoCount = monomerCount*sop.aminoCount + (monomerCount + 1)*linkerLength;
		tandem.bondCount = monomerCount*sop.bondCount + (monomerCount + 1)*(linkerLength + 1) - 2;
		tandem.nativeCount = monomerCount*sop.nativeCount;
		tandem.pairCount = (tandem.aminoCount * tandem.aminoCount - tandem.aminoCount) / 2 - tandem.nativeCount;
		if(covalentLJ == 0){
			tandem.pairCount -= tandem.bondCount;
		}

		printf("Total per %d-mer:\n%d amino acids\n%d covalent bonds\n%d native contacts\n%d pairs\n=======\n",
				monomerCount, tandem.aminoCount, tandem.bondCount, tandem.nativeCount, tandem.pairCount);

		tandem.aminos = (PDBAtom*)calloc(tandem.aminoCount, sizeof(PDBAtom));
		tandem.bonds = (CovalentBond*)calloc(tandem.bondCount, sizeof(CovalentBond));
		tandem.natives = (NativeContact*)calloc(tandem.nativeCount, sizeof(NativeContact));
		tandem.pairs = (PossiblePair*)calloc(tandem.pairCount, sizeof(PossiblePair));

		index = 0;
		PDBAtom linkerAA;
		strcpy(linkerAA.resName, "GLY");
		strcpy(linkerAA.name, "CA");
		linkerAA.chain = 'l';
		linkerAA.occupancy = sop.aminos[0].occupancy;
		linkerAA.beta = sop.aminos[0].beta;
		linkerAA.resid = 1;
		linkerAA.id = 1;

		float rx, ry, rz, r;
		float rxn, ryn, rzn;
		rx = sop.aminos[lastResid].x - sop.aminos[firstResid].x;
		ry = sop.aminos[lastResid].y - sop.aminos[firstResid].y;
		rz = sop.aminos[lastResid].z - sop.aminos[firstResid].z;
		if(tandemDirection == 1){
			r = sqrt(rx*rx + ry*ry + rz*rz);
			rxn = rx/r;
			ryn = ry/r;
			rzn = rz/r;
		} else {
			r = sqrt(tandemVectorX*tandemVectorX + tandemVectorY*tandemVectorY + tandemVectorZ*tandemVectorZ);
			rxn = tandemVectorX/r;
			ryn = tandemVectorY/r;
			rzn = tandemVectorZ/r;
		}

		// Copy monomer amino acids
		for(i = 0; i < monomerCount; i++){
			for(j = 0; j < sop.aminoCount; j++){
				tandem.aminos[i*(sop.aminoCount + linkerLength) + j + linkerLength] = sop.aminos[j];
				tandem.aminos[i*(sop.aminoCount + linkerLength) + j + linkerLength].x += i*(rx + rxn*a*(linkerLength + 1));
				tandem.aminos[i*(sop.aminoCount + linkerLength) + j + linkerLength].y += i*(ry + ryn*a*(linkerLength + 1));
				tandem.aminos[i*(sop.aminoCount + linkerLength) + j + linkerLength].z += i*(rz + rzn*a*(linkerLength + 1));
			}
		}

		// Add linker amino acids
		for(k = 0; k < linkerLength; k++){
			linkerAA.x = sop.aminos[firstResid].x - (k+1)*a*rxn;
			linkerAA.y = sop.aminos[firstResid].y - (k+1)*a*ryn;
			linkerAA.z = sop.aminos[firstResid].z - (k+1)*a*rzn;
			tandem.aminos[linkerLength - (k + 1)] = linkerAA;
		}
		for(i = 0; i < monomerCount; i++){
			for(k = 0; k < linkerLength; k++){
				linkerAA.x = sop.aminos[lastResid].x + (rx + (linkerLength + 1)*a*rxn)*i + (k+1)*a*rxn;
				linkerAA.y = sop.aminos[lastResid].y + (ry + (linkerLength + 1)*a*ryn)*i + (k+1)*a*ryn;
				linkerAA.z = sop.aminos[lastResid].z + (rz + (linkerLength + 1)*a*rzn)*i + (k+1)*a*rzn;
				tandem.aminos[(i+1)*(sop.aminoCount) + (i+1)*linkerLength + k] = linkerAA;
			}
		}

		// Copy monomer covalent bonds
		for(i = 0; i < monomerCount; i++){
			for(j = 0; j < sop.bondCount; j++){
				tandem.bonds[i*(sop.bondCount + linkerLength + 1) + linkerLength + j].i = sop.bonds[j].i + i*(sop.aminoCount + linkerLength) + linkerLength;
				tandem.bonds[i*(sop.bondCount + linkerLength + 1) + linkerLength + j].j = sop.bonds[j].j + i*(sop.aminoCount + linkerLength) + linkerLength;
				tandem.bonds[i*(sop.bondCount + linkerLength + 1) + linkerLength + j].r0 = sop.bonds[j].r0;
			}
		}
		// Adding covalent bonds for linker
		for(k = 0; k < linkerLength; k++){
			CovalentBond bond;
			bond.i = k;
			bond.j = k + 1;
			bond.r0 = a;
			tandem.bonds[k] = bond;
		}
		tandem.bonds[linkerLength - 1].j = firstResid + linkerLength;
		int kmax = linkerLength + 1;
		for(i = 0; i < monomerCount; i++){
			if(i == monomerCount - 1){
				kmax = linkerLength;
			}
			for(k = 0; k < kmax; k++){
				CovalentBond bond;
				bond.i = (i + 1)*(sop.aminoCount + linkerLength) + k - 1;
				bond.j = (i + 1)*(sop.aminoCount + linkerLength) + k;
				bond.r0 = a;
				tandem.bonds[(i + 1)*(sop.bondCount + linkerLength) + i + k] = bond;
				//printf("%d: %d-%d\n", (i + 1)*(sop.bondCount + linkerLength) + i + k, (i + 1)*(sop.aminoCount + linkerLength) + k -1, (i + 1)*(sop.aminoCount + linkerLength) + k);
			}
			tandem.bonds[(i + 1)*(sop.bondCount + linkerLength) + i].i = i*(sop.aminoCount + linkerLength) + linkerLength + lastResid;
			if(i != monomerCount - 1){
				tandem.bonds[(i + 1)*(sop.bondCount + linkerLength) + i + linkerLength].j = (i + 1)*(sop.aminoCount + linkerLength) + linkerLength + firstResid;
			}
		}

		// Copy monomer native contacts. There is no contacts between monomers and in the linker.
		for(i = 0; i < monomerCount; i++){
			for(j = 0; j < sop.nativeCount; j++){
				tandem.natives[i*sop.nativeCount + j].i = sop.natives[j].i + i*(sop.aminoCount + linkerLength) + linkerLength;
				tandem.natives[i*sop.nativeCount + j].j = sop.natives[j].j + i*(sop.aminoCount + linkerLength) + linkerLength;
				tandem.natives[i*sop.nativeCount + j].r0 = sop.natives[j].r0;
				tandem.natives[i*sop.nativeCount + j].eh = sop.natives[j].eh;
			}
		}

		// Buiding repulsion potential list. All pairs that are not in a native contact and don't form a covalent bond will be added to this list.
		int index = 0;
		int found = 0;
		for(i = 0; i < tandem.aminoCount; i++){
			for(j = i + 1; j < tandem.aminoCount; j++){
				if(covalentLJ == 0){
					k = 0;
					while(found == 0 && k < tandem.bondCount){
						if((tandem.bonds[k].i == i && tandem.bonds[k].j == j) ||
								tandem.bonds[k].j == i && tandem.bonds[k].i == j){
							found = 1;
						}
						k++;
					}
				}
				k = 0;
				while(found == 0 && k < tandem.nativeCount){
					if((tandem.natives[k].i == i && tandem.natives[k].j == j) ||
						tandem.natives[k].j == i && tandem.natives[k].i == j){
						found = 1;
					}
					k++;
				}
				if(found == 0){
					//printf("%d: %d - %d\n", index, i, j);
					tandem.pairs[index].i = i;
					tandem.pairs[index].j = j;
					index ++;
				}
				found = 0;
			}
		}

		sop = tandem;
	}
	sop.save(top_filename);
	FILE* test = fopen(coord_filename, "r");
	if(test != NULL){
		fclose(test);
		printf("Coordinates file '%s' exists. Overwrite (y/n)? ", coord_filename);
		char input;
		scanf("%c", &input);
		if(input == 'y'){
			savePDB(coord_filename, sop);
		}
	} else {
		savePDB(coord_filename, sop);
	}
}
