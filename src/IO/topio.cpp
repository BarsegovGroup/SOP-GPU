/*
 * topio.cpp
 *
 *  Created on: May 24, 2009
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Util/wrapper.h"
#include "topio.h"

#define BUF_SIZE 256

PDBAtom parseAtom(char* line);

enum TopologySection { NONE, ATOMS, BONDS, NATIVES, PAIRS };

TopologySection get_section(const char *line){
	if(strncmp(line, "[ atoms ]", 9) == 0){
        return ATOMS;
    }
	if(strncmp(line, "[ bonds ]", 9) == 0){
        return BONDS;
    }
    if(strncmp(line, "[ native ]", 10) == 0){
        return NATIVES;
    }
    if(strncmp(line, "[ pairs ]", 9) == 0){
        return PAIRS;
    }
    return NONE;
}

void SOP::load(const char* filename){
	int line;
	char buffer[BUF_SIZE];
	printf("Reading topology.\n");
	FILE* topFile = safe_fopen(filename, "r");

    TopologySection section = NONE;

    line = 0;
	while(safe_fgets(buffer, BUF_SIZE, topFile) != NULL){
        line++;
        if (buffer[0] == ';' || buffer[0] == '\n' || buffer[0] == '\0'){
            continue;
        }
        if (buffer[0] == '['){
            section = get_section(buffer);
            if (section == NONE){
                DIE("Malformed section header at %s:%d", filename, line);
            }
            continue;
        }
        switch(section) {
			TOPPair pair;
            case ATOMS:
                this->aminos.push_back(parseAtom(buffer));
                break;
            case BONDS:
                pair.parse(buffer);
                this->bonds.push_back(CovalentBond(pair));
                break;
            case NATIVES:
                pair.parse(buffer);
                this->natives.push_back(NativeContact(pair));
                break;
            case PAIRS:
                pair.parse(buffer);
                this->pairs.push_back(PossiblePair(pair));
                break;
            case NONE:
                DIE("Malformed topology: meaningful line at %s:%d outside any section", filename, line);
        }
    }
	fclose(topFile);
	printf("Done reading topology.\n");
}

void SOP::save(const char* filename) const {
	int i;
	FILE* topFile = safe_fopen(filename, "w");
	fprintf(topFile, "; Created by topio.cpp utility\n\n");

	fprintf(topFile, "[ atoms ]\n");
	fprintf(topFile, ";   nr       type  resnr residue  atom   cgnr     charge       mass\n");
	for(i = 0; i < this->aminos.size(); i++){
		fprintf(topFile, "%6d", i);
		fprintf(topFile, "%11s", this->aminos[i].name);
		fprintf(topFile, "%7d", this->aminos[i].resid);
		fprintf(topFile, "%7s", this->aminos[i].resName);
		fprintf(topFile, "%7s", this->aminos[i].name);
		fprintf(topFile, "%7c", this->aminos[i].chain);
		fprintf(topFile, "%11.2f", this->aminos[i].occupancy);
		fprintf(topFile, "%11.3f", this->aminos[i].beta);
		fprintf(topFile, "\n");
	}

	fprintf(topFile, "\n");

	fprintf(topFile, "[ bonds ]\n");
	fprintf(topFile, ";  ai    aj funct            c0            c1            c2            c3\n");
	for(i = 0; i < this->bonds.size(); i++){
		TOPPair pair;
		pair.i = this->bonds[i].i;
		pair.j = this->bonds[i].j;
		pair.func = 1;
		pair.c0 = this->bonds[i].r0;
		pair.c1 = 0.0f;
		pair.save(topFile);
	}

	fprintf(topFile, "\n");

	fprintf(topFile, "[ native ]\n");
	fprintf(topFile, ";  ai    aj funct            c0            c1            c2            c3\n");
	for(i = 0; i < this->natives.size(); i++){
		TOPPair pair;
		pair.i = this->natives[i].i;
		pair.j = this->natives[i].j;
		pair.func = 1;
		pair.c0 = this->natives[i].r0;
		pair.c1 = this->natives[i].eh;
		pair.c2 = 0.0f;
		pair.save(topFile);
	}

	fprintf(topFile, "\n");

	fprintf(topFile, "[ pairs ]\n");
	fprintf(topFile, ";  ai    aj funct            c0            c1            c2            c3\n");
	for(i = 0; i < this->pairs.size(); i++){
		TOPPair pair;
		pair.i = this->pairs[i].i;
		pair.j = this->pairs[i].j;
		pair.func = 1;
		pair.c0 = 0.0f;
		pair.save(topFile);
	}

	fclose(topFile);
}


void TOPPair::save(FILE* topFile) const {
	fprintf(topFile, "%5d", this->i);
	fprintf(topFile, " ");
	fprintf(topFile, "%5d", this->j);
	fprintf(topFile, " ");
	fprintf(topFile, "%5d", this->func);
	if(this->c0 != 0){
		fprintf(topFile, " ");
		fprintf(topFile, "%10.5f", this->c0);
		if(this->c1 != 0){
			fprintf(topFile, " ");
			fprintf(topFile, "%10.5f", this->c1);
			if(this->c2 != 0){
				fprintf(topFile, " ");
				fprintf(topFile, "%10.5f", this->c2);
				if(this->c3 != 0){
					fprintf(topFile, " ");
					fprintf(topFile, "%10.5f", this->c3);
				}
			}
		}
	}
	fprintf(topFile, "\n");
}

PDBAtom parseAtom(char* line){
    PDBAtom atom;
	char* pch;

	pch = strtok(line, " ");
	atom.id = atoi(pch);

	pch = strtok(NULL, " ");

	pch = strtok(NULL, " ");
	atom.resid = atoi(pch);

	pch = strtok(NULL, " ");
	strcpy(atom.resName, pch);

	pch = strtok(NULL, " ");
	strcpy(atom.name, pch);

	pch = strtok(NULL, " ");
	atom.chain = pch[0];

	pch = strtok(NULL, " ");
	atom.occupancy = atof(pch);

	pch = strtok(NULL, " ");
	atom.beta = atof(pch);
    
    atom.altLoc = ' ';

    return atom;
}

void TOPPair::parse(char *buffer){
    char *pch;

	pch = strtok(buffer, " ");
	this->i = atoi(pch);

	pch = strtok(NULL, " ");
	this->j = atoi(pch);

	pch = strtok(NULL, " ");
	this->func = atoi(pch);

	pch = strtok(NULL, " ");
	if(pch == NULL){
		return;
	}
	this->c0 = atof(pch);

	pch = strtok(NULL, " ");
	if(pch == NULL){
		return;
	}
	this->c1 = atof(pch);

	pch = strtok(NULL, " ");
	if(pch == NULL){
		return;
	}
	this->c2 = atof(pch);

	pch = strtok(NULL, " ");
	if(pch == NULL){
		return;
	}
	this->c3 = atof(pch);

}

void TOPPair::read(FILE* topFile){
	char buffer[BUF_SIZE];
	safe_fgets(buffer, BUF_SIZE, topFile);
    this->parse(buffer);
}

