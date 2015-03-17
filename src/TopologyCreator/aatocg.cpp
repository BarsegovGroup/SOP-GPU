/*
 * aatocg.cpp
 *
 *  Created on: Nov 28, 2012
 *      Author: kias
 */

#include "aatocg.h"

#define BUF_SIZE 1024


void readCGConfiguration(const char* filename, CGConfig *config){
	FILE* file = fopen(filename, "r");
	char buffer[BUF_SIZE];
	if(file != NULL){
		Residue resid;
		Bead bead;
		char* pch;
		while(fgets(buffer, BUF_SIZE, file) != NULL){
			if(strncmp(buffer, "MASS", 4) == 0){
				Mass mass;
				pch = strtok(buffer, " \t\n\r");
				pch = strtok(NULL, " \t\n\r");
				pch = strtok(NULL, " \t\n\r");
				mass.name = pch[0];
				pch = strtok(NULL, " \t\n\r");
				mass.mass = atof(pch);
				config->masses.push_back(mass);
			}
			if(strncmp(buffer, "RESI", 4) == 0){
				pch = strtok(buffer, " \t\n\r");
				pch = strtok(NULL, " \t\n\r");
				sprintf(resid.resname, "%s", pch);

			}
			if(strncmp(buffer, "ENDRESI", 4) == 0){
				config->residues.push_back(resid);
				resid.beads.clear();
			}
			if(strncmp(buffer, "BEAD", 4) == 0){
				pch = strtok(buffer, " \t\n\r");
				pch = strtok(NULL, " \t\n\r");
				sprintf(bead.name, "%s", pch);
				pch = strtok(NULL, " \t\n\r");
				sprintf(bead.type, "%s", pch);
			}
			if(strncmp(buffer, "ENDBEAD", 4) == 0){
				resid.beads.push_back(bead);
				bead.atomsCM.clear();
				bead.atomsRepresents.clear();
				bead.connectedTo.clear();
			}
			if(strncmp(buffer, "REPR", 4) == 0){
				pch = strtok(buffer, " \t\n\r");
				pch = strtok(NULL, " \t\n\r");
				while(pch != NULL){
					Atom atom;
					sprintf(atom.name, "%s", pch);
					bead.atomsRepresents.push_back(atom);
					pch = strtok(NULL, " \t\n\r");
				}
			}
			if(strncmp(buffer, "COOR", 4) == 0){
				pch = strtok(buffer, " \t\n\r");
				pch = strtok(NULL, " \t\n\r");
				while(pch != NULL){
					Atom atom;
					sprintf(atom.name, "%s", pch);
					bead.atomsCM.push_back(atom);
					pch = strtok(NULL, " \t\n\r");
				}
			}
			if(strncmp(buffer, "CONN", 4) == 0){
				pch = strtok(buffer, " \t\n\r");
				pch = strtok(NULL, " \t\n\r");
				while(pch != NULL){
					Atom atom;
					sprintf(atom.name, "%s", pch);
					bead.connectedTo.push_back(atom);
					pch = strtok(NULL, " \t\n\r");
				}
			}
			if(strncmp(buffer, "CHAR", 4) == 0){
				pch = strtok(buffer, " \t\n\r");
				pch = strtok(NULL, " \t\n\r");
				bead.charge = atof(pch);
			}
		}
	}
	unsigned int i, j, k, l;

	for(i = 0; i < config->residues.size(); i++){
		Residue resid = config->residues.at(i);
		for(j = 0; j < resid.beads.size(); j++){
			Bead bead = resid.beads.at(j);
			bead.mass = 0.0;
			for(k = 0; k < bead.atomsRepresents.size(); k++){
				for(l = 0; l < config->masses.size(); l++){
					if(config->masses.at(l).name == bead.atomsRepresents.at(k).name[0]){
						bead.mass += config->masses.at(l).mass;
					}
				}
			}
			config->residues.at(i).beads.at(j).mass = bead.mass;
		}
	}


	printf("Masses:\n");
	for(i = 0; i < config->masses.size(); i++){
		Mass mass = config->masses.at(i);
		printf("%c %f\n", mass.name, mass.mass);
	}
	for(i = 0; i < config->residues.size(); i++){
		Residue resid = config->residues.at(i);
		printf("\n\nResidue: %s\n", resid.resname);
		for(j = 0; j < resid.beads.size(); j++){
			Bead bead = resid.beads.at(j);
			printf("\nBead: %s-%s (q = %f, m = %f)\n", bead.name, bead.type, bead.charge, bead.mass);
			printf("Represents:");
			for(k = 0; k < bead.atomsRepresents.size(); k++){
				printf(" %s", bead.atomsRepresents.at(k).name);
			}
			printf("\n");
			printf("Coordinates:");
			for(k = 0; k < bead.atomsCM.size(); k++){
				printf(" %s", bead.atomsCM.at(k).name);
			}
			printf("\n");
			printf("Connected to:");
			for(k = 0; k < bead.connectedTo.size(); k++){
				printf(" %s", bead.connectedTo.at(k).name);
			}
			printf("\n");
		}
	}


	fclose(file);
}
