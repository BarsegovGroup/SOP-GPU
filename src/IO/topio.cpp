/*
 * topio.cpp
 *
 *  Created on: May 24, 2009
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../def_param.h"
#include "../Util/wrapper.h"
#include "topio.h"

#define buf_size 256

int countRows(FILE* topFile);
void readAtom(PDBAtom* atom, FILE* topFile);

void SOP::load(const char* filename){
	int i;
	printf("Reading topology.\n");
	FILE* topFile = safe_fopen(filename, "r");
	char buffer[buf_size];
	if (topFile != NULL ){
		while(safe_fgets(buffer, buf_size, topFile) != NULL){
			if(strncmp(buffer, "[ atoms ]", 9) == 0){
				printf("Counting atoms...\n");
				this->aminoCount = countRows(topFile);
				printf("%d found.\n", this->aminoCount);
			}
			if(strncmp(buffer, "[ bonds ]", 9) == 0){
				printf("Counting covalent bonds...\n");
				this->bondCount = countRows(topFile);
				printf("%d found.\n", this->bondCount);
			}
			if(strncmp(buffer, "[ native ]", 10) == 0){
				printf("Counting native contacts...\n");
				this->nativeCount = countRows(topFile);
				printf("%d found.\n", this->nativeCount);
			}
			if(strncmp(buffer, "[ pairs ]", 9) == 0){
				printf("Counting pairs...\n");
				this->pairCount = countRows(topFile);
				printf("%d found.\n", this->pairCount);
			}
		}
	} else {
		printf("ERROR: can't find topology file '%s'.\n", filename);
		exit(-1);
	}


	this->aminos = (PDBAtom*)calloc(this->aminoCount, sizeof(PDBAtom));
	this->bonds = (CovalentBond*)calloc(this->bondCount, sizeof(CovalentBond));
	this->natives = (NativeContact*)calloc(this->nativeCount, sizeof(NativeContact));
	this->pairs = (PossiblePair*)calloc(this->pairCount, sizeof(PossiblePair));

	rewind(topFile);
	while(safe_fgets(buffer, buf_size, topFile) != NULL){
		if(strncmp(buffer, "[ atoms ]", 9) == 0){
			safe_fgets(buffer, buf_size, topFile);
			for(i = 0; i < this->aminoCount; i++){
				readAtom(&this->aminos[i], topFile);
			}
		}
		if(strncmp(buffer, "[ bonds ]", 9) == 0){
			safe_fgets(buffer, buf_size, topFile);
			for(i = 0; i < this->bondCount; i++){
				TOPPair pair;
                pair.read(topFile);
				this->bonds[i].i = pair.i;
				this->bonds[i].j = pair.j;
				this->bonds[i].r0 = pair.c0;
			}
		}
		if(strncmp(buffer, "[ native ]", 10) == 0){
			safe_fgets(buffer, buf_size, topFile);
			for(i = 0; i < this->nativeCount; i++){
				TOPPair pair; 
                pair.read(topFile);
				this->natives[i].i = pair.i;
				this->natives[i].j = pair.j;
				this->natives[i].r0 = pair.c0;
				this->natives[i].eh = pair.c1;

			}
		}
		if(strncmp(buffer, "[ pairs ]", 9) == 0){
			safe_fgets(buffer, buf_size, topFile);
			for(i = 0; i < this->pairCount; i++){
				TOPPair pair;
                pair.read(topFile);
				this->pairs[i].i = pair.i;
				this->pairs[i].j = pair.j;
			}
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
	for(i = 0; i < this->aminoCount; i++){
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
	for(i = 0; i < this->bondCount; i++){
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
	for(i = 0; i < this->nativeCount; i++){
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
	for(i = 0; i < this->pairCount; i++){
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


int countRows(FILE* topFile){
	char buffer[buf_size];
	char* pch;
	int result = 0;
	int skip;
	char* eof;
	do{
		eof = safe_fgets(buffer, buf_size, topFile);
		pch = strtok(buffer, " ");
		//printf("'%s'\n", pch);
		if(strcmp(pch, ";") != 0){
			result ++;
			skip = 0;
		} else {
			skip = 1;
		}
	} while((strcmp(pch,"0")==0 ||atoi(pch) != 0 || skip == 1) && eof != NULL);
	return result - 1;
}

void readAtom(PDBAtom* atom, FILE* topFile){
	char buffer[buf_size];
	char* pch;
	safe_fgets(buffer, buf_size, topFile);

	pch = strtok(buffer, " ");
	atom->id = atoi(pch);

	pch = strtok(NULL, " ");

	pch = strtok(NULL, " ");
	atom->resid = atoi(pch);

	pch = strtok(NULL, " ");
	strcpy(atom->resName, pch);

	pch = strtok(NULL, " ");
	strcpy(atom->name, pch);

	pch = strtok(NULL, " ");
	atom->chain = pch[0];

	pch = strtok(NULL, " ");
	atom->occupancy = atof(pch);

	pch = strtok(NULL, " ");
	atom->beta = atof(pch);

}

void TOPPair::read(FILE* topFile){
	char buffer[buf_size];
	char* pch;
	safe_fgets(buffer, buf_size, topFile);

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

