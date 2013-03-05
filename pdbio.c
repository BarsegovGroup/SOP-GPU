#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "def_param.h"

/*
 * pdbio.c
 *
 *  Created on: Nov 9, 2008
 *      Author: zhmurov
 */

#define buf_size 80

void readPDB(char* pdb_filename, PDB* pdbdata);
void printAtom(Atom atomdata);

void savePDB(char* filename);
void readCoord(char* filename);

void readPDB(char* pdb_filename, PDB* pdbdata){
	printf("Reading %s.\n", pdb_filename);
	int ss_count = 0, atoms_count = 0;
	char buffer[buf_size], *s;
	FILE* file = fopen(pdb_filename, "r");
	if ( file != NULL ){
		while(fgets(buffer, buf_size, file) != NULL){
			if(strncmp(buffer,"SSBOND",6) == 0){
				ss_count++;
			}
			if(strncmp(buffer, "ATOM", 4) == 0){
				atoms_count++;
			}
		}
		printf("Found %d atoms, %d ss-bounds.\n", atoms_count, ss_count);

		pdbdata->atomCount = atoms_count;
		pdbdata->ssCount = ss_count;
		pdbdata->atoms = (Atom*)malloc(atoms_count*sizeof(Atom));
		pdbdata->ssbonds = (SSBond*)malloc(ss_count*sizeof(SSBond));

		int current_atom = 0;
		int current_ss = 0;

		rewind(file);

		while(fgets(buffer, buf_size, file) != NULL){
			char* pch = strtok(buffer, " ");
			if(strcmp(pch, "SSBOND") == 0){
				char chain;
				char resid[4];
				strncpy(&chain, &buffer[15], 1);
				strncpy(resid, &buffer[17], 4);
				pdbdata->ssbonds[current_ss].chain1 = chain;
				pdbdata->ssbonds[current_ss].resid1 = atoi(resid);
				strncpy(&chain, &buffer[29], 1);
				strncpy(resid, &buffer[31], 4);
				pdbdata->ssbonds[current_ss].chain2 = chain;
				pdbdata->ssbonds[current_ss].resid2 = atoi(resid);

#ifdef DEBUG
				printf("SS #%d: %c%d - %c%d\n", current_ss, pdbdata->ssbonds[current_ss].chain1, pdbdata->ssbonds[current_ss].resid1,
						pdbdata->ssbonds[current_ss].chain2, pdbdata->ssbonds[current_ss].resid2);
#endif
				current_ss++;
			}
			if(strcmp(pch, "ATOM") == 0){
				char id[6];
				char atomName[5];
				char resName[4], chain;
				char resid[5];
				char x[9], y[9], z[9];
				char occupancy[7];
				char beta[7];
				strncpy(id, &buffer[6], 5);
				id[5] = '\0';
				strncpy(atomName, &buffer[12], 4);
				atomName[4] = '\0';
				strncpy(resName, &buffer[17], 3);
				resName[3] = '\0';
				strncpy(&chain, &buffer[21], 1);
				strncpy(resid, &buffer[22], 4);
				resid[4] = '\0';
				//printf("%s\n", resid);
				strncpy(x, &buffer[30], 8);
				x[8] = '\0';
				strncpy(y, &buffer[38], 8);
				y[8] = '\0';
				strncpy(z, &buffer[46], 8);
				z[8] = '\0';
				strncpy(occupancy, &buffer[54], 6);
				occupancy[6] = '\0';
				strncpy(beta, &buffer[60], 6);
				beta[6] = '\0';
				strcpy(pdbdata->atoms[current_atom].name, strtok(atomName, " "));
				pdbdata->atoms[current_atom].name[4] = 0;
				pdbdata->atoms[current_atom].chain = chain;
				pdbdata->atoms[current_atom].resid = atoi(resid);
				strcpy(pdbdata->atoms[current_atom].resName, resName);
				pdbdata->atoms[current_atom].resName[3] = 0;
				pdbdata->atoms[current_atom].id = atoi(id);
				pdbdata->atoms[current_atom].x = atof(x);
				pdbdata->atoms[current_atom].y = atof(y);
				pdbdata->atoms[current_atom].z = atof(z);
				pdbdata->atoms[current_atom].occupancy = atof(occupancy);
				pdbdata->atoms[current_atom].beta = atof(beta);
				pdbdata->atoms[current_atom].ssbond_aminoid = -1;
#ifdef DEBUG1
		printf("ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
				pdbdata->atoms[current_atom].id,
				pdbdata->atoms[current_atom].name,
				pdbdata->atoms[current_atom].resName,
				pdbdata->atoms[current_atom].chain,
				pdbdata->atoms[current_atom].resid,
				pdbdata->atoms[current_atom].x,
				pdbdata->atoms[current_atom].y,
				pdbdata->atoms[current_atom].z,
				pdbdata->atoms[current_atom].occupancy,
				pdbdata->atoms[current_atom].beta);
#endif
				current_atom ++;
			}

		}
		/*for(int i = 0; i < sop.aminoCount; i++){
			for(int j = 0; j < sop.aminoCount; j++){
				for(int k = 0; k < ss_count; k++){
					if((sop.aminos[i].resid == SSbonds[k].resid1) && (sop.aminos[i].chain == SSbonds[k].chain1) &&
							(sop.aminos[j].resid == SSbonds[k].resid2) && (sop.aminos[j].chain == SSbonds[k].chain2)){
						printf("Creating SS bound b/w resid %d and %d.\n", sop.aminos[i].id, sop.aminos[j].id);
						sop.aminos[i].ssbond_aminoid = sop.aminos[j].id;
						sop.aminos[j].ssbond_aminoid = sop.aminos[i].id;
					}
				}
			}
		}*/
	printf("Done reading '%s'.\n", pdb_filename);
	fclose(file);
	} else {
		perror(pdb_filename);
		exit(0);
	}
}

void readCoord(char* filename){
	//printf("Reading coordinates from '%s'.\n", filename);
	char buffer[buf_size];
	FILE* file = fopen(filename, "r");
	if(file == NULL){
		printf("ERROR: Coordinates file %s can not be found.\n", filename);
		exit(0);
	}
	int index = 0;
	while(fgets(buffer, buf_size, file) != NULL){
		char* pch = strtok(buffer, " ");
		if(strcmp(pch, "ATOM") == 0){
			char x[8], y[8], z[8];
			strncpy(x, &buffer[30], 8);
			strncpy(y, &buffer[38], 8);
			strncpy(z, &buffer[46], 8);
			sop.aminos[index].x = atof(x);
			sop.aminos[index].y = atof(y);
			sop.aminos[index].z = atof(z);
#ifdef DEBUG
			printf("%d: %f, %f, %f\n", index, sop.aminos[index].x, sop.aminos[index].y, sop.aminos[index].z);
#endif
			index++;
		}
	}
	if(index == 0){
		printf("ERROR: Can't read pdb file.\n");
		fclose(file);
		exit(0);
	}
	fclose(file);
	//printf("Done reading coordinates.\n");
}

void savePDB(char* filename){
	FILE* file = fopen(filename, "w");
	int i;
	for(i = 0; i < sop.aminoCount; i++){
		fprintf(file, "ATOM %6d  CA  %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
								i + 1,
								sop.aminos[i].resName,
								sop.aminos[i].chain,
								sop.aminos[i].resid,
								sop.aminos[i].x,
								sop.aminos[i].y,
								sop.aminos[i].z,
								sop.aminos[i].occupancy,
								sop.aminos[i].beta);
	}
	if(sop.additionalAminosCount != 0){
		for(i = 0; i < sop.additionalAminosCount; i++){
				fprintf(file, "ATOM %6d  CA  %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
										sop.aminoCount + i + 1,
										sop.additionalAminos[i].resName,
										sop.additionalAminos[i].chain,
										sop.additionalAminos[i].resid,
										sop.additionalAminos[i].x,
										sop.additionalAminos[i].y,
										sop.additionalAminos[i].z,
										sop.additionalAminos[i].occupancy,
										sop.additionalAminos[i].beta);
			}
	}
	int lastAmino = -1;
	for(int i = 0; i < sop.bondCount; i++){
		if(lastAmino != sop.bonds[i].i){
			if(lastAmino != -1){
				fprintf(file, "\n");
			}
			fprintf(file, "CONECT%5d", sop.bonds[i].i + 1);
			lastAmino = sop.bonds[i].i;
		}
		fprintf(file, "%5d", sop.bonds[i].j + 1);
	}
	fprintf(file, "\nEND");
	fclose(file);
}

void printAtom(Atom atomdata){
	printf("ATOM  %5d  CA  %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
			atomdata.id,
			atomdata.name,
			atomdata.chain,
			atomdata.resid,
			atomdata.x,
			atomdata.y,
			atomdata.z,
			atomdata.occupancy,
			atomdata.beta);
}
