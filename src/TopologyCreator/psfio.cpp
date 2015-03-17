/*
 * psfio.c
 *
 *	Utility to read/write PSF (Protein Structure File) files.
 *
 *  Created on: Mar 19, 2009
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#define UNWRAP
#ifdef UNWRAP
# define safe_fopen fopen
# define safe_fgets fgets
# define safe_fread fread
# define DIE(format, ...) do{ printf(format, ##__VA_ARGS__); exit(-1); }while(0);
#else
# include "../Util/wrapper.h"
#endif
#include "psfio.h"

#define BUF_SIZE 512
//#define DEBUGPSFIO

/*
 * Private methods
 */

void readAtoms(FILE* psfFile, PSF* psfData);
void readBonds(FILE* psfFile, PSF* psfData);
void readAngles(FILE* psfFile, PSF* psfData);
void readDihedrals(FILE* psfFile, PSF* psfData);
void readImpropers(FILE* psfFile, PSF* psfData);
void readCMAP(FILE* psfFile, PSF* psfData);
void readExclusions(FILE* psfFile, PSF* psfData);

void writeAtoms(FILE* psfFile, PSF* psfData);
void writeBonds(FILE* psfFile, PSF* psfData);
void writeAngles(FILE* psfFile, PSF* psfData);
void writeDihedrals(FILE* psfFile, PSF* psfData);
void writeImpropers(FILE* psfFile, PSF* psfData);
void writeCMAP(FILE* psfFile, PSF* psfData);
void writeExclusions(FILE* psfFile, PSF* psfData);

/*
 * Read the PSF file
 * Parameters:
 * 		filename: name of the file to read
 * 		psfData: pointer to a PSF object to write data into
 */
void readPSF(const char* filename, PSF* psfData){
	printf("Reading structure data from '%s'.\n", filename);
	FILE* psfFile = safe_fopen(filename, "r");
	if(psfFile != NULL){
		readAtoms(psfFile, psfData);
		readBonds(psfFile, psfData);
		readAngles(psfFile, psfData);
		readDihedrals(psfFile, psfData);
		readImpropers(psfFile, psfData);
		readCMAP(psfFile, psfData);
		readExclusions(psfFile, psfData);
		fclose(psfFile);
		printf("Done reading psf.\n");
	} else {
		DIE("Error opening '%s': %s", filename, strerror(errno));
	}
}

/*
 * Read atoms section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to read
 * 		psfData: pointer to the 'PSF' object to read data into
 */
void readAtoms(FILE* psfFile, PSF* psfData){
	int i,j;
	char * pch;
	char buffer[BUF_SIZE];
	do{
		safe_fgets(buffer, BUF_SIZE, psfFile);
	} while (strstr(buffer, "!NATOM") == NULL);
	pch = strtok(buffer, " ");
	psfData->natom = atoi(pch);
	printf("Found %d atoms.\n", psfData->natom);
	psfData->atoms = (PSFAtom*)malloc(psfData->natom*sizeof(PSFAtom));
	for(i = 0; i < psfData->natom; i++){
		safe_fgets(buffer, BUF_SIZE, psfFile);
		// Read ID
		pch = strtok(buffer, " ");
		psfData->atoms[i].id = atoi(pch);
		// Read SegmentID
		pch = strtok (NULL, " ");
		for(j = 0; j < PSF_SEGMENT_SIZE; j++){
			psfData->atoms[i].segment[j] = ' ';
		}
		psfData->atoms[i].segment[PSF_SEGMENT_SIZE] = '\0';
		strcpy(psfData->atoms[i].segment, pch);
		// Read ResNum
		pch = strtok (NULL, " ");
		psfData->atoms[i].resid = atoi(pch);
		// Read ResName
		pch = strtok (NULL, " ");
		for(j = 0; j < PSF_RESNAME_SIZE; j++){
			psfData->atoms[i].resName[j] = ' ';
		}
		psfData->atoms[i].resName[PSF_RESNAME_SIZE] = '\0';
		strcpy(psfData->atoms[i].resName, pch);
		// Read atom name
		pch = strtok (NULL, " ");
		for(j = 0; j < PSF_NAME_SIZE; j++){
			psfData->atoms[i].name[j] = ' ';
		}
		psfData->atoms[i].name[PSF_NAME_SIZE] = '\0';
		strcpy(psfData->atoms[i].name, pch);
		// Read atom type
		pch = strtok (NULL, " ");
		for(j = 0; j < PSF_TYPE_SIZE; j++){
			psfData->atoms[i].type[j] = ' ';
		}
		psfData->atoms[i].type[PSF_TYPE_SIZE] = '\0';
		strcpy(psfData->atoms[i].type, pch);
		// Read atom charge
		pch = strtok (NULL, " ");
		psfData->atoms[i].q = atof(pch);
		// Read atom mass
		pch = strtok (NULL, " ");
		psfData->atoms[i].m = atof(pch);
	}
	rewind(psfFile);
#ifdef DEBUGPSFIO
	for(i = 0; i < psfData->natom; i++){
		printf("%5d %s %5d %s %s %s %f %f\n",
				psfData->atoms[i].id, psfData->atoms[i].segment, psfData->atoms[i].resid,
				psfData->atoms[i].resName, psfData->atoms[i].name, psfData->atoms[i].type,
				psfData->atoms[i].q, psfData->atoms[i].m);
	}
#endif
	float q = 0.0f;
	float m = 0.0f;
	for(i = 0; i < psfData->natom; i++){
		m += psfData->atoms[i].m;
		q += psfData->atoms[i].q;
	}
	printf("Total mass: %f\n", m);
	printf("Total charge: %f\n", q);

}

/*
 * Read bonds section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to read
 * 		psfData: pointer to the 'PSF' object to read data into
 */
void readBonds(FILE* psfFile, PSF* psfData){
	int i;
	char * pch;
	char buffer[BUF_SIZE];
	do{
		safe_fgets(buffer, BUF_SIZE, psfFile);
	} while (strstr(buffer, "!NBOND") == NULL);
	pch = strtok(buffer, " ");
	psfData->nbond = atoi(pch);

	printf("Found %d bonds.\n", psfData->nbond);
	psfData->bonds = (PSFBond*)malloc(psfData->nbond*sizeof(PSFBond));
	i = 0;
	while(i < psfData->nbond){
		if(i % 4 == 0){
			safe_fgets(buffer, BUF_SIZE, psfFile);
			pch = strtok(buffer, " ");
		} else {
			pch = strtok(NULL, " ");
		}
		psfData->bonds[i].i = atoi(pch);
		pch = strtok(NULL, " ");
		psfData->bonds[i].j = atoi(pch);
#ifdef DEBUGPSFIO
			printf("%d - %d\n", psfData->bonds[i].i, psfData->bonds[i].j);
#endif
		i++;
	}
/*	for(i = 0; i < psfData->nbond; i++){
		safe_fgets(buffer, 17, psfFile);
		pch = strtok(buffer, " ");
		if(strcmp(pch, "\n") == 0){
			i--;
		} else {
			psfData->bonds[i].i = atoi(pch);
			pch = strtok(NULL, " ");
			psfData->bonds[i].j = atoi(pch);
#ifdef DEBUGPSFIO
			printf("%d - %d\n", psfData->bonds[i].i, psfData->bonds[i].j);
#endif
		}
	}*/
	rewind(psfFile);
}

/*
 * Read angles section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to read
 * 		psfData: pointer to the 'PSF' object to read data into
 */
void readAngles(FILE* psfFile, PSF* psfData){
	int i;
	char * pch;
	char buffer[BUF_SIZE];
	do{
		safe_fgets(buffer, BUF_SIZE, psfFile);
	} while (strstr(buffer, "!NTHETA") == NULL);
	pch = strtok(buffer, " ");
	psfData->ntheta = atoi(pch);

	printf("Found %d angles.\n", psfData->ntheta);
	psfData->angles = (PSFAngle*)malloc(psfData->ntheta*sizeof(PSFAngle));

	i = 0;
	while(i < psfData->ntheta){
		if(i % 3 == 0){
			safe_fgets(buffer, BUF_SIZE, psfFile);
			pch = strtok(buffer, " ");
		} else {
			pch = strtok(NULL, " ");
		}
		psfData->angles[i].i = atoi(pch);
		pch = strtok(NULL, " ");
		psfData->angles[i].j = atoi(pch);
		pch = strtok(NULL, " ");
		psfData->angles[i].k = atoi(pch);
#ifdef DEBUGPSFIO
			printf("%d - %d - %d\n", psfData->angles[i].i, psfData->angles[i].j, psfData->angles[i].k);
#endif
		i++;

	}

	/*for(i = 0; i < psfData->ntheta; i++){
		safe_fgets(buffer, 25, psfFile);
		pch = strtok(buffer, " ");
		if(strcmp(pch, "\n") == 0){
			i--;
		} else {
			psfData->angles[i].i = atoi(pch);
			pch = strtok(NULL, " ");
			psfData->angles[i].j = atoi(pch);
			pch = strtok(NULL, " ");
			psfData->angles[i].k = atoi(pch);
#ifdef DEBUGPSFIO
			printf("%d - %d - %d\n", psfData->angles[i].i, psfData->angles[i].j, psfData->angles[i].k);
#endif
		}
	}*/
	rewind(psfFile);
}


/*
 * Read dihedrals section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to read
 * 		psfData: pointer to the 'PSF' object to read data into
 */
void readDihedrals(FILE* psfFile, PSF* psfData){
	int i;
	char * pch;
	char buffer[BUF_SIZE];
	do{
		safe_fgets(buffer, BUF_SIZE, psfFile);
	} while (strstr(buffer, "!NPHI") == NULL);
	pch = strtok(buffer, " ");
	psfData->nphi = atoi(pch);

	printf("Found %d dihedrals.\n", psfData->nphi);
	psfData->dihedrals = (PSFDihedral*)malloc(psfData->nphi*sizeof(PSFDihedral));
	i = 0;
	while(i < psfData->nphi){
		if(i % 2 == 0){
			safe_fgets(buffer, BUF_SIZE, psfFile);
			pch = strtok(buffer, " ");
		} else {
			pch = strtok(NULL, " ");
		}
		psfData->dihedrals[i].i = atoi(pch);
		pch = strtok(NULL, " ");
		psfData->dihedrals[i].j = atoi(pch);
		pch = strtok(NULL, " ");
		psfData->dihedrals[i].k = atoi(pch);
		pch = strtok(NULL, " ");
		psfData->dihedrals[i].l = atoi(pch);

#ifdef DEBUGPSFIO
			printf("%d - %d - %d - %d\n", psfData->dihedrals[i].i, psfData->dihedrals[i].j, psfData->dihedrals[i].k, psfData->dihedrals[i].l);
#endif
		i++;
}
/*	for(i = 0; i < psfData->nphi; i++){
		safe_fgets(buffer, 33, psfFile);
		pch = strtok(buffer, " ");
		if(strcmp(pch, "\n") == 0){
			i--;
		} else {
			psfData->dihedrals[i].i = atoi(pch);
			pch = strtok(NULL, " ");
			psfData->dihedrals[i].j = atoi(pch);
			pch = strtok(NULL, " ");
			psfData->dihedrals[i].k = atoi(pch);
			pch = strtok(NULL, " ");
			psfData->dihedrals[i].l = atoi(pch);
#ifdef DEBUGPSFIO
			printf("%d - %d - %d - %d\n", psfData->dihedrals[i].i, psfData->dihedrals[i].j, psfData->dihedrals[i].k, psfData->dihedrals[i].l);
#endif
		}
	}*/
	rewind(psfFile);
}

/*
 * Read impropers section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to read
 * 		psfData: pointer to the 'PSF' object to read data into
 */
void readImpropers(FILE* psfFile, PSF* psfData){
	int i;
	char * pch;
	char buffer[BUF_SIZE];
	do{
		safe_fgets(buffer, BUF_SIZE, psfFile);
	} while (strstr(buffer, "!NIMPHI") == NULL);
	pch = strtok(buffer, " ");
	psfData->nimphi = atoi(pch);

	printf("Found %d impropers.\n", psfData->nimphi);
	psfData->impropers = (PSFImproper*)malloc(psfData->nimphi*sizeof(PSFImproper));
	i = 0;
	while(i < psfData->nimphi){
		if(i % 2 == 0){
			safe_fgets(buffer, BUF_SIZE, psfFile);
			pch = strtok(buffer, " ");
		} else {
			pch = strtok(NULL, " ");
		}
		psfData->impropers[i].i = atoi(pch);
		pch = strtok(NULL, " ");
		psfData->impropers[i].j = atoi(pch);
		pch = strtok(NULL, " ");
		psfData->impropers[i].k = atoi(pch);
		pch = strtok(NULL, " ");
		psfData->impropers[i].l = atoi(pch);
#ifdef DEBUGPSFIO
			printf("%d - %d - %d - %d\n", psfData->impropers[i].i, psfData->impropers[i].j, psfData->impropers[i].k, psfData->impropers[i].l);
#endif
		i++;
	}
	/*for(i = 0; i < psfData->nimphi; i++){
		safe_fgets(buffer, 33, psfFile);
		pch = strtok(buffer, " ");
		if(strcmp(pch, "\n") == 0){
			i--;
		} else {
			psfData->impropers[i].i = atoi(pch);
			pch = strtok(NULL, " ");
			psfData->impropers[i].j = atoi(pch);
			pch = strtok(NULL, " ");
			psfData->impropers[i].k = atoi(pch);
			pch = strtok(NULL, " ");
			psfData->impropers[i].l = atoi(pch);
#ifdef DEBUGPSFIO
			printf("%d - %d - %d - %d\n", psfData->impropers[i].i, psfData->impropers[i].j, psfData->impropers[i].k, psfData->impropers[i].l);
#endif
		}
	}*/
	rewind(psfFile);
}

/*
 * Read CMAP corrections section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to read
 * 		psfData: pointer to the 'PSF' object to read data into
 */
void readCMAP(FILE* psfFile, PSF* psfData){
	int i;
	char * pch;
	char buffer[BUF_SIZE];
	do{
		safe_fgets(buffer, BUF_SIZE, psfFile);
	} while (strstr(buffer, "!NCRTERM") == NULL && !feof(psfFile));
	if (feof(psfFile)) {
		psfData->ncmap = 0;
		printf("CMAP section not found.\n");
	} else {
		pch = strtok(buffer, " ");
		psfData->ncmap = atoi(pch);

		printf("Found %d CMAP corrections.\n", psfData->ncmap);
		psfData->cmaps = (PSFCMAP*)malloc(psfData->ncmap*sizeof(PSFCMAP));
		for(i = 0; i < psfData->ncmap; i++){
			safe_fgets(buffer, BUF_SIZE, psfFile);
			pch = strtok(buffer, " ");
			if(strcmp(pch, "\n") == 0){
				i--;
			} else {
				psfData->cmaps[i].i1 = atoi(pch);
				pch = strtok(NULL, " ");
				psfData->cmaps[i].j1 = atoi(pch);
				pch = strtok(NULL, " ");
				psfData->cmaps[i].k1 = atoi(pch);
				pch = strtok(NULL, " ");
				psfData->cmaps[i].l1 = atoi(pch);
				pch = strtok(NULL, " ");
				psfData->cmaps[i].i2 = atoi(pch);
				pch = strtok(NULL, " ");
				psfData->cmaps[i].j2 = atoi(pch);
				pch = strtok(NULL, " ");
				psfData->cmaps[i].k2 = atoi(pch);
				pch = strtok(NULL, " ");
				psfData->cmaps[i].l2 = atoi(pch);
#ifdef DEBUGPSFIO
				printf("%d - %d - %d - %d - %d - %d - %d - %d\n",
					psfData->cmaps[i].i1, psfData->cmaps[i].j1, psfData->cmaps[i].k1, psfData->cmaps[i].l1,
					psfData->cmaps[i].i2, psfData->cmaps[i].j2, psfData->cmaps[i].k2, psfData->cmaps[i].l2);
#endif
			}
		}
	}
	rewind(psfFile);
}

void readExclusions(FILE* psfFile, PSF* psfData){
	int i;
	char * pch;
	char buffer[BUF_SIZE];
	do{
		safe_fgets(buffer, BUF_SIZE, psfFile);
	} while (strstr(buffer, "!NNB") == NULL);
	pch = strtok(buffer, " ");
	psfData->nnb = atoi(pch);

	printf("Found %d explicit exclusions.\n", psfData->nnb);
	psfData->nbExclusions = (int*)malloc(psfData->nnb*sizeof(int));
	psfData->nbExclusionsCounts = (int*)malloc(psfData->natom*sizeof(int));
	if(psfData->nnb != 0){
		for(i = 0; i < psfData->nnb; i++){
			if(i % 8 == 0){
				safe_fgets(buffer, BUF_SIZE, psfFile);
				pch = strtok(buffer, " \n\r\t");
			} else {
				pch = strtok(NULL, " \n\r\t");
			}
			psfData->nbExclusions[i] = atoi(pch);
#ifdef DEBUGPSFIO
			printf("%d  ", psfData->nbExclusions[i]);
			if((i+1) % 8 == 0){
				printf("\n");
			}
#endif
		}
#ifdef DEBUGPSFIO
		printf("\nExclusion groups:\n");
#endif
		for(i = 0; i < psfData->natom; i++){
			if(i % 8 == 0){
				safe_fgets(buffer, BUF_SIZE, psfFile);
				pch = strtok(buffer, " \n\r\t");
			} else {
				pch = strtok(NULL, " \n\r\t");
			}
			psfData->nbExclusionsCounts[i] = atoi(pch);
#ifdef DEBUGPSFIO
			printf("%d  ", psfData->nbExclusionsCounts[i]);
			if((i+1) % 8 == 0){
				printf("\n");
			}
#endif
		}
#ifdef DEBUGPSFIO
		printf("\n");
#endif
	} else {
		for(i = 0; i < psfData->natom; i++){
			psfData->nbExclusionsCounts[i] = 0;
		}
	}
	rewind(psfFile);
}

/*
 * Write the PSF data into PSF file
 * Parameters:
 * 		filename: name of the file to write into (will be overwritten/created)
 * 		psfData: pointer to a PSF object to get data from
 */
void writePSF(const char* filename, PSF* psfData){
	FILE* psfFile = safe_fopen(filename, "w");

	fprintf(psfFile, "PSF\n");
	fprintf(psfFile, "\n");
	fprintf(psfFile, "       1 !NTITLE\n");
	fprintf(psfFile, " REMARKS generated by psfio.c utility\n");
	fprintf(psfFile, "\n");
	writeAtoms(psfFile, psfData);
	writeBonds(psfFile, psfData);
	writeAngles(psfFile, psfData);
	writeDihedrals(psfFile, psfData);
	writeImpropers(psfFile, psfData);
	writeCMAP(psfFile, psfData);
	fclose(psfFile);
}

/*
 * Writes atoms section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to write into
 * 		psfData: pointer to the 'PSF' object to get data from
 */
void writeAtoms(FILE* psfFile, PSF* psfData){

	int i;
	fprintf(psfFile, "%8d !NATOM\n", psfData->natom);
	for(i = 0; i < psfData->natom; i++){
		fprintf(psfFile, "%8d ", psfData->atoms[i].id);
		fprintf(psfFile, "%-5s", psfData->atoms[i].segment);
		fprintf(psfFile, "%-4d ", psfData->atoms[i].resid);
		fprintf(psfFile, "%-4s", psfData->atoms[i].resName);
		fprintf(psfFile, " ");
		fprintf(psfFile, "%-5s", psfData->atoms[i].name);
		fprintf(psfFile, "%-5s ", psfData->atoms[i].type);
		if(psfData->atoms[i].q >= 0){
			fprintf(psfFile, " %7.6f", psfData->atoms[i].q);
		} else {
			fprintf(psfFile, "%8.6f", psfData->atoms[i].q);
		}
		fprintf(psfFile, "      ");
		fprintf(psfFile, "%8.4f", psfData->atoms[i].m);
		fprintf(psfFile, "\n");
	}
	fprintf(psfFile, "\n");
}

/*
 * Writes bonds section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to write into
 * 		psfData: pointer to the 'PSF' object to get data from
 */
void writeBonds(FILE* psfFile, PSF* psfData){
	int i;
	fprintf(psfFile, "%8d !NBOND: bonds\n", psfData->nbond);
	for(i = 0; i < psfData->nbond; i++){
		fprintf(psfFile, "%8d%8d",
				psfData->bonds[i].i,
				psfData->bonds[i].j);
		if(i != 0 && i%4 == 3){
			fprintf(psfFile, "\n");
		}
	}
	fprintf(psfFile, "\n\n");
}

/*
 * Writes angles section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to write into
 * 		psfData: pointer to the 'PSF' object to get data from
 */
void writeAngles(FILE* psfFile, PSF* psfData){
	int i;
	fprintf(psfFile, "%8d !NTHETA: angles\n", psfData->ntheta);
	for(i = 0; i < psfData->ntheta; i++){
		fprintf(psfFile, "%8d%8d%8d",
				psfData->angles[i].i,
				psfData->angles[i].j,
				psfData->angles[i].k);
		if(i != 0 && i%3 == 2){
			fprintf(psfFile, "\n");
		}
	}
	fprintf(psfFile, "\n\n");
}

/*
 * Writes dihedrals section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to write into
 * 		psfData: pointer to the 'PSF' object to get data from
 */
void writeDihedrals(FILE* psfFile, PSF* psfData){
	int i;
	fprintf(psfFile, "%8d !NPHI: dihedrals\n", psfData->nphi);
	for(i = 0; i < psfData->nphi; i++){
		fprintf(psfFile, "%8d%8d%8d%8d",
				psfData->dihedrals[i].i,
				psfData->dihedrals[i].j,
				psfData->dihedrals[i].k,
				psfData->dihedrals[i].l);
		if(i != 0 && i%2 == 1){
			fprintf(psfFile, "\n");
		}
	}

	fprintf(psfFile, "\n\n");
}

/*
 * Writes impropers section of a PSF file
 * Parameters:
 * 		psfFile: pointer to the file to write into
 * 		psfData: pointer to the 'PSF' object to get data from
 */
void writeImpropers(FILE* psfFile, PSF* psfData){
	int i;
	fprintf(psfFile, "%8d !IMNPHI: impropers\n", psfData->nimphi);
	for(i = 0; i < psfData->nimphi; i++){
		fprintf(psfFile, "%8d%8d%8d%8d",
				psfData->impropers[i].i,
				psfData->impropers[i].j,
				psfData->impropers[i].k,
				psfData->impropers[i].l);
		if(i != 0 && i%2 == 1){
			fprintf(psfFile, "\n");
		}
	}
	fprintf(psfFile, "\n\n");
}

void writeCMAP(FILE* psfFile, PSF* psfData){
	int i;
	fprintf(psfFile, "%8d !NCRTERM: cross-terms\n", psfData->ncmap);
	for(i = 0; i < psfData->ncmap; i++){
		fprintf(psfFile, "%8d%8d%8d%8d%8d%8d%8d%8d\n",
				psfData->cmaps[i].i1, psfData->cmaps[i].j1, psfData->cmaps[i].k1, psfData->cmaps[i].l1,
				psfData->cmaps[i].i2, psfData->cmaps[i].j2, psfData->cmaps[i].k2, psfData->cmaps[i].l2);
	}
	fprintf(psfFile, "\n\n");
}

void writeExclusions(FILE* psfFile, PSF* psfData){

}
