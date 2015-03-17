/*
 * topio.c
 *
 *  Created on: May 24, 2009
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define UNWRAP
#ifdef UNWRAP
# define safe_fopen fopen
# define safe_fgets fgets
# define safe_fread fread
# define DIE(format, ...) do{ printf(format, ##__VA_ARGS__); exit(-1); }while(0);
#else
# include "../Util/wrapper.h"
#endif

#include "topio.h"

#define BUF_SIZE 256

int countRowsInTOP(FILE* topFile);
TOPAtom readAtomLineFromTOP(FILE* topFile);
TOPPair readPairLineFromTOP(FILE* topFile);
void savePair(FILE* topFile, TOPPair pair);

int readTOP(const char* filename, TOPData* topData){
	printf("Reading topology from '%s'.\n", filename);
	FILE* topFile = safe_fopen(filename, "r");
	char buffer[BUF_SIZE];
	if (topFile != NULL ){
		while(safe_fgets(buffer, BUF_SIZE, topFile) != NULL){
			if(strstr(buffer, "[ atoms ]") != 0){
				printf("Counting atoms...\n");
				topData->atomCount = countRowsInTOP(topFile);
				printf("%d found.\n", topData->atomCount);
			}
			if(strstr(buffer, "[ bonds ]") != 0){
				printf("Counting bonds...\n");
				topData->bondCount = countRowsInTOP(topFile);
				printf("%d found.\n", topData->bondCount);
			}
			if(strstr(buffer, "[ pairs ]") != 0){
				printf("Counting pairs...\n");
				topData->pairsCount = countRowsInTOP(topFile);
				printf("%d found.\n", topData->pairsCount);
			}
			if(strstr(buffer, "[ angles ]") != 0){
				printf("Counting angles...\n");
				topData->angleCount = countRowsInTOP(topFile);
				printf("%d found.\n", topData->angleCount);
			}
			if(strstr(buffer, "[ dihedrals ]") != 0){
				printf("Counting dihedrals...\n");
				topData->dihedralCount = countRowsInTOP(topFile);
				printf("%d found.\n", topData->dihedralCount);
			}
		}
		topData->atoms = (TOPAtom*)calloc(topData->atomCount, sizeof(TOPAtom));
		topData->bonds = (TOPPair*)calloc(topData->bondCount, sizeof(TOPPair));
		topData->pairs = (TOPPair*)calloc(topData->pairsCount, sizeof(TOPPair));
		topData->angles = (TOPAngle*)calloc(topData->angleCount, sizeof(TOPAngle));
		topData->dihedrals = (TOPDihedral*)calloc(topData->dihedralCount, sizeof(TOPDihedral));
	} else {
		DIE("ERROR: cant find topology file '%s'.", filename);
	}

	rewind(topFile);
	int count = 0;
	while(safe_fgets(buffer, BUF_SIZE, topFile) != NULL){
		if(strstr(buffer, "[ atoms ]") != 0){
			count = 0;
			while(count < topData->atomCount){
				topData->atoms[count] = readAtomLineFromTOP(topFile);
				if(topData->atoms[count].id != -1){
					count++;
				}
			}
		}
		if(strstr(buffer, "[ bonds ]") != 0){
			count = 0;
			while(count < topData->bondCount){
				topData->bonds[count] = readPairLineFromTOP(topFile);
				if(topData->bonds[count].i != -1){
					count++;
				}
			}
		}
		if(strstr(buffer, "[ pairs ]") != 0){
			count = 0;
			while(count < topData->pairsCount){
				topData->pairs[count] = readPairLineFromTOP(topFile);
				if(topData->pairs[count].i != -1){
					count++;
				}
			}
		}
		if(strstr(buffer, "[ angles ]") != 0){
//TODO
		}
		if(strstr(buffer, "[ dihedrals ]") != 0){
//TODO
		}
	}

	fclose(topFile);
	printf("Done reading the topology section.\n");
	return count;
}

void writeTOP(const char* filename, TOPData* topData){
	int i;
	FILE* topFile = safe_fopen(filename, "w");
	fprintf(topFile, "; Created by topio.c utility\n\n");
	fprintf(topFile, "[ atoms ]\n");
	fprintf(topFile, ";   nr       type  resnr residue  atom   cgnr     charge       mass\n");
	for(i = 0; i < topData->atomCount; i++){
		fprintf(topFile, "%6d", i);
		fprintf(topFile, "%11s", topData->atoms[i].type);
		fprintf(topFile, "%7d", topData->atoms[i].resid);
		fprintf(topFile, "%7s", topData->atoms[i].resName);
		fprintf(topFile, "%7s", topData->atoms[i].name);
		fprintf(topFile, "%7c", topData->atoms[i].chain);
		fprintf(topFile, "%11.2f", topData->atoms[i].charge);
		fprintf(topFile, "%11.3f", topData->atoms[i].mass);
		fprintf(topFile, "\n");
	}

	fprintf(topFile, "\n");

	fprintf(topFile, "[ bonds ]\n");
	fprintf(topFile, ";  ai    aj funct            c0            c1            c2            c3\n");
	for(i = 0; i < topData->bondCount; i++){
		savePair(topFile, topData->bonds[i]);
	}

	fprintf(topFile, "\n");

	fprintf(topFile, "[ native ]\n");
	fprintf(topFile, ";  ai    aj funct            c0            c1            c2            c3\n");
	for(i = 0; i < topData->nativesCount; i++){
		savePair(topFile, topData->natives[i]);
	}

	fprintf(topFile, "\n");

	fprintf(topFile, "[ pairs ]\n");
	fprintf(topFile, ";  ai    aj funct            c0            c1            c2            c3\n");
	for(i = 0; i < topData->pairsCount; i++){
		savePair(topFile, topData->pairs[i]);
	}

	fprintf(topFile, "\n");

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


int countRowsInTOP(FILE* topFile){
	char buffer[BUF_SIZE];
	char* pch;
	int result = 0;
	int skip;
	char* eof;
	do{
		eof = safe_fgets(buffer, BUF_SIZE, topFile);
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

TOPAtom readAtomLineFromTOP(FILE* topFile){
	char buffer[BUF_SIZE];
	char* pch;
	TOPAtom atom;
	safe_fgets(buffer, BUF_SIZE, topFile);

	if(strncmp(buffer, ";", 1) != 0){
		pch = strtok(buffer, " ");
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
		atom.charge = atof(pch);

		pch = strtok(NULL, " ");
		atom.mass = atof(pch);
	} else {
		atom.id = -1;
	}

	return atom;
}

TOPPair readPairLineFromTOP(FILE* topFile){

	TOPPair pair;
	char buffer[BUF_SIZE];
	char* pch;
	safe_fgets(buffer, BUF_SIZE, topFile);
	if(strncmp(buffer, ";", 1) != 0){
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
	} else {
		pair.i = -1;
	}
	return pair;
}

