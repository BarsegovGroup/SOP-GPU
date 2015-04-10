/*
 * gsop-top.c
 *
 *  Created on: Oct 21, 2010
 *      Author: zhmurov
 */

#include "sop.hpp"
#include "IO/configreader.h"
#include "IO/topio.h"
#include "IO/pdbio.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define VERSION "trunk"

using namespace configreader; // << Temporary kludge to keep soptop working without getting deep into it

void createModel();

int main(int argc, char *argv[]){

	printf("==========================\n");
	printf("SOP-GPU Topology creator version %s\n", VERSION);
	printf("==========================\n");

	if(argc < 1){
		printf("ERROR: Configuration file should be specified.\n");
		exit(-1);
	}

	parseParametersFile(argv[1]); //, argc, argv);

	createModel();
}


void createModel(){

    std::string ehString = getParameterAs<std::string>("eh", "1.5");
	if(ehString == "O"){
		printf("Taking eh's values from occupancy column.\n");
		eh = -1.0;
	} else if(ehString == "B"){
		printf("Taking eh's values from beta column.\n");
		eh = -2.0;
	} else {
		eh = atof(ehString.c_str());
		if(eh == 0){
			printf("WARNING: Value of eh in a configuration file is not valid. Should be:\n "
					" - Positive integer, or\n"
					" - 'O' to take values from occupancy column, or\n"
					" - 'B' to take from beta column.\n"
					"eh is set to 0.\n");
		}
	}
	createTandem = getParameterAs<bool>("createTandem", false);
	if(createTandem){
		linkerLength = getParameterAs<int>("linkerLength");
		monomerCount = getParameterAs<int>("monomerCount");
        std::string tandemDirectionString = getParameterAs<std::string>("tandemDirection", "endToEnd");
		if(tandemDirectionString == "endToEnd"){
			tandemDirection = 1;
			printf("Direction of tandem linker is a monomer end-to-end direction.\n");
		} else
		if(tandemDirectionString == "vector"){
			tandemDirection = 2;
			tandemVectorX = getParameterAs<float>("tandemVectorX");
			tandemVectorY = getParameterAs<float>("tandemVectorY");
			tandemVectorZ = getParameterAs<float>("tandemVectorZ");
			printf("Tandem linker direction is (%3.2f, %3.2f, %3.2f) vector.\n", tandemVectorX, tandemVectorY, tandemVectorZ);
		} else {
			tandemDirection = 1;
			printf("WARNING: Direction of tandem linker is not specified or wrong. Assuming monomer end-to-end direction.\n");
		}
		firstResid = getParameterAs<int>("fixedEnd");
		lastResid = getParameterAs<int>("pulledEnd");
	}
	pdb_filename = getParameterAs<std::string>("structure");
	top_filename = getParameterAs<std::string>("topology");
	coord_filename = getParameterAs<std::string>("coordinates");


	R_limit_bond = getParameterAs<float>("R_limit_bond", 8.0f);
	SC_limit_bond = getParameterAs<float>("SC_limit_bond", 0.0f);
	covalentLJ = getParameterAs<bool>("covalentRepulsion", false);
	a = getParameterAs<float>("a", 3.8f);
	pairs_threshold = getParameterAs<float>("pairs_threshold", 1000.0f);
	covalent_cutoff = getParameterAs<float>("covalent_cutoff", 10.0f);
	SS_cutoff = getParameterAs<float>("SS_cutoff", 10.0f);

	int i, j, k, i1, i2;
	covLinkersCount = getParameterAs<int>("covLinkersCount", 0);
	if(covLinkersCount > 0){
		covLinkers = (CovalentLinker*)calloc(covLinkersCount, sizeof(CovalentLinker));
		covLinkerCutoff = getParameterAs<float>("covLinkersCutoff", 5.0f);
		char *covLinker;
		char covLinkerName[50];
		char* pch;
		for(i = 0; i < covLinkersCount; i++){
			sprintf(covLinkerName, "covLinker%d", i+1);
            const std::string tmp_str = getParameterAs<std::string>(covLinkerName);
            covLinker = strdup(tmp_str.c_str());
			pch = strtok(covLinker, "-\t");
			strncpy(covLinkers[i].residName1, pch, 3);
			covLinkers[i].resid1 = atoi(&pch[3]);
			pch = strtok(NULL, "\n\r\t ");
			strncpy(covLinkers[i].residName2, pch, 3);
			covLinkers[i].resid2 = atoi(&pch[3]);
            free(covLinker);
			printf("Covalent linker between residues %s%d and %s%d will be added.\n",
					covLinkers[i].residName1, covLinkers[i].resid1,
					covLinkers[i].residName2, covLinkers[i].resid2);
		}
	}

	pdbdata.read(pdb_filename.c_str());
	printf("Creating SOP-model...\n");

	int sop_aminoCount = 0;
	//Counting number of Calphas
	for(i = 0; i < pdbdata.atoms.size(); i++){
		if(strcmp(pdbdata.atoms[i].name, "CA") == 0){
			sop_aminoCount ++;
		}
	}
	printf("Found %d residues.\n", sop_aminoCount);

	pdbRef = (int*)malloc(sop.aminos.size()*sizeof(int));
	firstAtomInResid = (int*)malloc(sop.aminos.size()*sizeof(int));
	lastAtomInResid = (int*)malloc(sop.aminos.size()*sizeof(int));
	aminoRef = (int*)malloc(pdbdata.atoms.size()*sizeof(int));
	sop.aminos.resize(sop_aminoCount);
	/*nativeContacts = (char**)malloc(sop.aminos.size()*sop.aminos.size()*sizeof(char));
	for(i = 0; i < sop.aminos.size(); i++){
		nativeContacts[i] = (char*)malloc(sop.aminos.size()*sizeof(char));
	}*/

	int index = 0;
	for(i = 0; i < pdbdata.atoms.size(); i++){
		if(strcmp(pdbdata.atoms[i].name, "CA") == 0){
			pdbRef[index] = i;
			sop.aminos[index] = pdbdata.atoms[i];
			index ++;
		}

	}

	index = -1;
	int lastAminoResid = -1;// pdbdata.atoms[0].resid;
	char lastAminoChain = 0;//pdbdata.atoms[0].chain;
	for(i = 0; i < pdbdata.atoms.size(); i++){
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
	lastAtomInResid[sop.aminos.size() - 1] = pdbdata.atoms.size() - 1;
	for(i = 0; i < sop.aminos.size() - 1; i++){
		while(aminoRef[index] == aminoRef[index+1]){
			index ++;
		}
		index ++;
		lastAtomInResid[i] = index;
		firstAtomInResid[i+1] = index+1;
	}

	/*for(i = 0; i < sop.aminos.size(); i++){
		printf("Residue %d : %d - %d\n", i, firstAtomInResid[i], lastAtomInResid[i]);
	}*/


/*	for(i = 0; i < sop.aminos.size(); i++){
		for(j = 0; j < sop.aminos.size(); j++){
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


	int sop_bondCount = 0;
	int sop_nativeCount = 0;
	int sop_pairCount = 0;
	for(i = 0; i < sop.aminos.size(); i++){
		for(j = i + 1; j < sop.aminos.size(); j++){
			if(checkCovalent(i, j)){
				sop_bondCount ++;
			}
			if(checkNative(i, j)){
				sop_nativeCount ++;
			}
			if(checkPossiblePairs(i, j)){
				sop_pairCount ++;
			}
		}
		printf("Amino: %d. Bonds: %d. Native contacts: %d. Pairs: %d. (%4.2f%% completed)\r",
				i,
				sop_bondCount,
				sop_nativeCount,
				sop_pairCount,
				((float)i+1)*100.0f/((float)(sop.aminos.size())));
	}
	printf("\n");

	printf("Total per monomer:\n%ld amino acids\n%d covalent bonds\n%d native contacts\n%d pairs\n=======\n",
			sop.aminos.size(), sop_bondCount, sop_nativeCount, sop_pairCount);


	sop.bonds.resize(sop_bondCount);
	sop.natives.resize(sop_nativeCount);
	sop.pairs.resize(sop_pairCount);



	index = 0;
	for(i = 0; i < pdbdata.atoms.size(); i++){
		if(strcmp(pdbdata.atoms[i].name, "CA") == 0){
			memcpy(&sop.aminos[index], &pdbdata.atoms[i], sizeof(PDBAtom));
			index ++;
		}
	}
	sop_bondCount = 0;
	sop_nativeCount = 0;
	sop_pairCount = 0;

	for(i = 0; i < sop.aminos.size(); i++){
		for(j = i + 1; j < sop.aminos.size(); j++){
			if(checkCovalent(i, j)){
				sop.bonds[sop_bondCount].i = i;
				sop.bonds[sop_bondCount].j = j;
				sop.bonds[sop_bondCount].r0 = getR0(i, j);
				sop_bondCount ++;
			}
			if(checkNative(i, j)){
				sop.natives[sop_nativeCount].i = i;
				sop.natives[sop_nativeCount].j = j;
				sop.natives[sop_nativeCount].r0 = getR0(i, j);
				sop.natives[sop_nativeCount].eh = getEh(i, j);
				sop_nativeCount ++;
			}
			if(checkPossiblePairs(i, j)){
				sop.pairs[sop_pairCount].i = i;
				sop.pairs[sop_pairCount].j = j;
				sop_pairCount ++;
			}
		}
		printf("Creating model ... %5.2f%% completed\r", ((float)i+1)*100.0f/((float)(sop.aminos.size())));
	}
	printf("\n");

	if(createTandem){


		printf("Building a tandem of %d\n", monomerCount);

		int tandem_aminoCount = monomerCount*sop.aminos.size() + (monomerCount + 1)*linkerLength;
		int tandem_bondCount = monomerCount*sop.bonds.size() + (monomerCount + 1)*(linkerLength + 1) - 2;
		int tandem_nativeCount = monomerCount*sop.natives.size();
		int tandem_pairCount = (tandem_aminoCount * tandem_aminoCount - tandem_aminoCount) / 2 - tandem_nativeCount;

		if(covalentLJ == 0){
			tandem_pairCount -= tandem_bondCount;
		}


		tandem.aminos.resize(tandem_aminoCount);
		tandem.bonds.resize(tandem_bondCount);
		tandem.natives.resize(tandem_nativeCount);
		tandem.pairs.resize(tandem_pairCount);

		printf("Total per %d-mer:\n%ld amino acids\n%ld covalent bonds\n%ld native contacts\n%ld pairs\n=======\n",
				monomerCount, tandem.aminos.size(), tandem.bonds.size(), tandem.natives.size(), tandem.pairs.size());

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
			for(j = 0; j < sop.aminos.size(); j++){
				tandem.aminos[i*(sop.aminos.size() + linkerLength) + j + linkerLength] = sop.aminos[j];
				tandem.aminos[i*(sop.aminos.size() + linkerLength) + j + linkerLength].x += i*(rx + rxn*a*(linkerLength + 1));
				tandem.aminos[i*(sop.aminos.size() + linkerLength) + j + linkerLength].y += i*(ry + ryn*a*(linkerLength + 1));
				tandem.aminos[i*(sop.aminos.size() + linkerLength) + j + linkerLength].z += i*(rz + rzn*a*(linkerLength + 1));
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
				tandem.aminos[(i+1)*(sop.aminos.size()) + (i+1)*linkerLength + k] = linkerAA;
			}
		}

		// Copy monomer covalent bonds
		for(i = 0; i < monomerCount; i++){
			for(j = 0; j < sop.bonds.size(); j++){
				tandem.bonds[i*(sop.bonds.size() + linkerLength + 1) + linkerLength + j].i = sop.bonds[j].i + i*(sop.aminos.size() + linkerLength) + linkerLength;
				tandem.bonds[i*(sop.bonds.size() + linkerLength + 1) + linkerLength + j].j = sop.bonds[j].j + i*(sop.aminos.size() + linkerLength) + linkerLength;
				tandem.bonds[i*(sop.bonds.size() + linkerLength + 1) + linkerLength + j].r0 = sop.bonds[j].r0;
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
				bond.i = (i + 1)*(sop.aminos.size() + linkerLength) + k - 1;
				bond.j = (i + 1)*(sop.aminos.size() + linkerLength) + k;
				bond.r0 = a;
				tandem.bonds[(i + 1)*(sop.bonds.size() + linkerLength) + i + k] = bond;
				//printf("%d: %d-%d\n", (i + 1)*(sop.bonds.size() + linkerLength) + i + k, (i + 1)*(sop.aminos.size() + linkerLength) + k -1, (i + 1)*(sop.aminos.size() + linkerLength) + k);
			}
			tandem.bonds[(i + 1)*(sop.bonds.size() + linkerLength) + i].i = i*(sop.aminos.size() + linkerLength) + linkerLength + lastResid;
			if(i != monomerCount - 1){
				tandem.bonds[(i + 1)*(sop.bonds.size() + linkerLength) + i + linkerLength].j = (i + 1)*(sop.aminos.size() + linkerLength) + linkerLength + firstResid;
			}
		}

		// Copy monomer native contacts. There is no contacts between monomers and in the linker.
		for(i = 0; i < monomerCount; i++){
			for(j = 0; j < sop.natives.size(); j++){
				tandem.natives[i*sop.natives.size() + j].i = sop.natives[j].i + i*(sop.aminos.size() + linkerLength) + linkerLength;
				tandem.natives[i*sop.natives.size() + j].j = sop.natives[j].j + i*(sop.aminos.size() + linkerLength) + linkerLength;
				tandem.natives[i*sop.natives.size() + j].r0 = sop.natives[j].r0;
				tandem.natives[i*sop.natives.size() + j].eh = sop.natives[j].eh;
			}
		}

		// Buiding repulsion potential list. All pairs that are not in a native contact and don't form a covalent bond will be added to this list.
		int index = 0;
		int found = 0;
		for(i = 0; i < tandem.aminos.size(); i++){
			for(j = i + 1; j < tandem.aminos.size(); j++){
				if(covalentLJ == 0){
					k = 0;
					while(found == 0 && k < tandem.bonds.size()){
						if((tandem.bonds[k].i == i && tandem.bonds[k].j == j) ||
								tandem.bonds[k].j == i && tandem.bonds[k].i == j){
							found = 1;
						}
						k++;
					}
				}
				k = 0;
				while(found == 0 && k < tandem.natives.size()){
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
	sop.save(top_filename.c_str());
	FILE* test = fopen(coord_filename.c_str(), "r");
	if(test!=NULL){
        fclose(test);
		printf("Coordinates file '%s' exists. Overwrite (y/n)? ", coord_filename.c_str());
		char input;
		while (scanf("%c", &input) != 1);
		if(input == 'y'){
			savePDB(coord_filename.c_str(), sop);
		}
	} else {
		savePDB(coord_filename.c_str(), sop);
	}
}
