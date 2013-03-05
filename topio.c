/*
 * topio.cpp
 *
 *  Created on: May 24, 2009
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "def_param.h"
#include "topio.h"

#define buf_size 256

void loadTOP(char* filename);
void saveTOP(char* filename);

void savePair(FILE* topFile, TOPPair pair);

int countRows(FILE* topFile);
void readAtom(Atom* atom, FILE* topFile);
TOPPair readPair(FILE* topFile);

void loadTOP(char* filename){
	int i;
	printf("Reading topology.\n");
	FILE* topFile = fopen(filename, "r");
	char buffer[buf_size];
	if (topFile != NULL ){
		while(fgets(buffer, buf_size, topFile) != NULL){
			if(strncmp(buffer, "[ atoms ]", 9) == 0){
				printf("Counting atoms...\n");
				sop.aminoCount = countRows(topFile);
				printf("%d found.\n", sop.aminoCount);
			}
			if(strncmp(buffer, "[ bonds ]", 9) == 0){
				printf("Counting covalent bonds...\n");
				sop.bondCount = countRows(topFile);
				printf("%d found.\n", sop.bondCount);
			}
			if(strncmp(buffer, "[ native ]", 10) == 0){
				printf("Counting native contacts...\n");
				sop.nativeCount = countRows(topFile);
				printf("%d found.\n", sop.nativeCount);
			}
			if(strncmp(buffer, "[ pairs ]", 9) == 0){
				printf("Counting pairs...\n");
				sop.pairCount = countRows(topFile);
				printf("%d found.\n", sop.pairCount);
			}
		}
	} else {
		printf("ERROR: can't find topology file '%s'.\n", filename);
		exit(-1);
	}


	sop.aminos = (Atom*)calloc(sop.aminoCount, sizeof(Atom));
	sop.bonds = (CovalentBond*)calloc(sop.bondCount, sizeof(CovalentBond));
	sop.natives = (NativeContact*)calloc(sop.nativeCount, sizeof(NativeContact));
	sop.pairs = (PossiblePair*)calloc(sop.pairCount, sizeof(PossiblePair));

	rewind(topFile);
	while(fgets(buffer, buf_size, topFile) != NULL){
		if(strncmp(buffer, "[ atoms ]", 9) == 0){
			fgets(buffer, buf_size, topFile);
			for(i = 0; i < sop.aminoCount; i++){
				readAtom(&sop.aminos[i], topFile);
			}
		}
		if(strncmp(buffer, "[ bonds ]", 9) == 0){
			fgets(buffer, buf_size, topFile);
			for(i = 0; i < sop.bondCount; i++){
				TOPPair pair = readPair(topFile);
				sop.bonds[i].i = pair.i;
				sop.bonds[i].j = pair.j;
				sop.bonds[i].r0 = pair.c0;
			}
		}
		if(strncmp(buffer, "[ native ]", 10) == 0){
			fgets(buffer, buf_size, topFile);
			for(i = 0; i < sop.nativeCount; i++){
				TOPPair pair = readPair(topFile);
				sop.natives[i].i = pair.i;
				sop.natives[i].j = pair.j;
				sop.natives[i].r0 = pair.c0;
				sop.natives[i].eh = pair.c1;

			}
		}
		if(strncmp(buffer, "[ pairs ]", 9) == 0){
			fgets(buffer, buf_size, topFile);
			for(i = 0; i < sop.pairCount; i++){
				TOPPair pair = readPair(topFile);
				sop.pairs[i].i = pair.i;
				sop.pairs[i].j = pair.j;
			}
		}
	}

	fclose(topFile);
	printf("Done reading topology.\n");
}

void saveTOP(char* filename){
	int i;
	FILE* topFile = fopen(filename, "w");
	fprintf(topFile, "; Created by topio.c utility\n\n");
	fprintf(topFile, "[ atoms ]\n");
	fprintf(topFile, ";   nr       type  resnr residue  atom   cgnr     charge       mass\n");
	for(i = 0; i < sop.aminoCount; i++){
		fprintf(topFile, "%6d", i);
		fprintf(topFile, "%11s", sop.aminos[i].name);
		fprintf(topFile, "%7d", sop.aminos[i].resid);
		fprintf(topFile, "%7s", sop.aminos[i].resName);
		fprintf(topFile, "%7s", sop.aminos[i].name);
		fprintf(topFile, "%7c", sop.aminos[i].chain);
		fprintf(topFile, "%11.2f", sop.aminos[i].occupancy);
		fprintf(topFile, "%11.3f", sop.aminos[i].beta);
		fprintf(topFile, "\n");
	}

	fprintf(topFile, "\n");

	fprintf(topFile, "[ bonds ]\n");
	fprintf(topFile, ";  ai    aj funct            c0            c1            c2            c3\n");
	for(i = 0; i < sop.bondCount; i++){
		TOPPair pair;
		pair.i = sop.bonds[i].i;
		pair.j = sop.bonds[i].j;
		pair.func = 1;
		pair.c0 = sop.bonds[i].r0;
		pair.c1 = 0.0f;
		savePair(topFile, pair);
	}

	fprintf(topFile, "\n");

	fprintf(topFile, "[ native ]\n");
	fprintf(topFile, ";  ai    aj funct            c0            c1            c2            c3\n");
	for(i = 0; i < sop.nativeCount; i++){
		TOPPair pair;
		pair.i = sop.natives[i].i;
		pair.j = sop.natives[i].j;
		pair.func = 1;
		pair.c0 = sop.natives[i].r0;
		pair.c1 = sop.natives[i].eh;
		pair.c2 = 0.0f;
		savePair(topFile, pair);
	}

	fprintf(topFile, "\n");

	fprintf(topFile, "[ pairs ]\n");
	fprintf(topFile, ";  ai    aj funct            c0            c1            c2            c3\n");
	for(i = 0; i < sop.pairCount; i++){
		TOPPair pair;
		pair.i = sop.pairs[i].i;
		pair.j = sop.pairs[i].j;
		pair.func = 1;
		pair.c0 = 0.0f;
		savePair(topFile, pair);
	}


	fclose(topFile);

}


void savePair(FILE* topFile, TOPPair pair){
	fprintf(topFile, "%5d", pair.i);
	fprintf(topFile, " ");
	fprintf(topFile, "%5d", pair.j);
	fprintf(topFile, " ");
	fprintf(topFile, "%5d", pair.func);
	if(pair.c0 != 0){
		fprintf(topFile, " ");
		fprintf(topFile, "%10.5f", pair.c0);
		if(pair.c1 != 0){
			fprintf(topFile, " ");
			fprintf(topFile, "%10.5f", pair.c1);
			if(pair.c2 != 0){
				fprintf(topFile, " ");
				fprintf(topFile, "%10.5f", pair.c2);
				if(pair.c3 != 0){
					fprintf(topFile, " ");
					fprintf(topFile, "%10.5f", pair.c3);
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
		eof = fgets(buffer, buf_size, topFile);
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

void readAtom(Atom* atom, FILE* topFile){
	char buffer[buf_size];
	char* pch;
	fgets(buffer, buf_size, topFile);

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

TOPPair readPair(FILE* topFile){
	TOPPair pair;
	char buffer[buf_size];
	char* pch;
	fgets(buffer, buf_size, topFile);

	pch = strtok(buffer, " ");
	pair.i = atoi(pch);

	pch = strtok(NULL, " ");
	pair.j = atoi(pch);

	pch = strtok(NULL, " ");
	pair.func = atoi(pch);

	pch = strtok(NULL, " ");
	if(pch == NULL){
		return pair;
	}
	pair.c0 = atof(pch);

	pch = strtok(NULL, " ");
	if(pch == NULL){
		return pair;
	}
	pair.c1 = atof(pch);

	pch = strtok(NULL, " ");
	if(pch == NULL){
		return pair;
	}
	pair.c2 = atof(pch);

	pch = strtok(NULL, " ");
	if(pch == NULL){
		return pair;
	}
	pair.c3 = atof(pch);
	return pair;
}

