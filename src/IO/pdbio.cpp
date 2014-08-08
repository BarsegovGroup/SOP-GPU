#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "pdbio.h"
#ifdef UNWRAP
# define safe_fopen fopen
# define safe_fgets fgets
# define safe_fread fread
# define DIE(format, ...) do{ printf(format, ##__VA_ARGS__); exit(-1); }while(0);
#else
# include "../Util/wrapper.h"
#endif

/*
 * pdbio.c
 *
 *  Created on: Nov 9, 2008
 *      Author: zhmurov
 */

#define BUF_SIZE 80

/*
 * Parses data from PDB (Protein Data Bank) file format into PDB object.
 * Only ATOM and SSBOND entries are considered, some fields are ignored
 * (marked befor corresponding method).
 *
 * For detailed information about PDB file format, please visit
 * http://www.wwpdb.org/documentation/format23/sect6.html
 *
 * Parameters:
 * 		filename: name of a pdb file to parse
 */

void PDB::read(const char* filename){
	printf("Reading %s.\n", filename);
	int ss_count = 0, atoms_count = 0;
	char buffer[BUF_SIZE];
	FILE* file = safe_fopen(filename, "r");
	if ( file != NULL ){
		while(fgets(buffer, BUF_SIZE, file) != NULL){
			if(strncmp(buffer,"SSBOND",6) == 0){
				ss_count++;
			}
			if(strncmp(buffer, "ATOM", 4) == 0){
				atoms_count++;
			}
		}
		printf("Found %d atoms.\n", atoms_count);

		this->atomCount = atoms_count;
		this->ssCount = ss_count;
		this->atoms = (PDBAtom*)malloc(atoms_count*sizeof(PDBAtom));
		this->ssbonds = (PDBSSBond*)malloc(ss_count*sizeof(PDBSSBond));

		int current_atom = 0;
		int current_ss = 0;

		rewind(file);

		while(fgets(buffer, BUF_SIZE, file) != NULL){
			char* pch = strtok(buffer, " ");
			if(strcmp(pch, "SSBOND") == 0){
				this->ssbonds[current_ss].parse(buffer);
				current_ss++;
			}
			if(strcmp(pch, "ATOM") == 0){
                this->atoms[current_atom].parse(buffer);
				current_atom ++;
			}

		}
	printf("Done reading '%s'.\n", filename);
	fclose(file);
	} else {
        DIE("Error opening '%s': %s", filename, strerror(errno));
	}
}

/*
 * Parses single line of 'ATOM' entry from pdb file.
 * ATOM entry format in PDB:
 *
 * COLUMNS      DATA TYPE        FIELD      DEFINITION
 * ------------------------------------------------------
 *  1 -  6      Record name      "ATOM    "
 *  7 - 11      Integer          id		    Atom serial number.
 * 13 - 16      Atom             name       Atom name.
 * 17           Character        altLoc     Alternate location indicator.
 * 18 - 20      Residue name     resName    Residue name.
 * 22           Character        chainID    Chain identifier.
 * 23 - 26      Integer          resSeq     Residue sequence number.
 * 27           AChar            iCode      Code for insertion of residues. (Ignored)
 * 31 - 38      Real(8.3)        x          Orthogonal coordinates for X in Angstroms
 * 39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in Angstroms
 * 47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in Angstroms
 * 55 - 60      Real(6.2)        occupancy  Occupancy.
 * 61 - 66      Real(6.2)        tempFactor Temperature factor.
 * 77 - 78      LString(2)       element    Element symbol, right-justified. (Ignored)
 * 79 - 80      LString(2)       charge     Charge on the atom. (Ignored)
 *
 * For more details visit http://www.wwpdb.org/documentation/format23/sect9.html
 *
 *
 * Method parameters:
 * 		pdbData:	pointer to pdbData object to read into
 * 		line: line from the file to parse
 * 		currentSSBond: indicates location in array of ssbonds
 * 			where new bond should be saved
 */

void PDBAtom::parse(const char* line){
	char id[6];
	char atomName[5];
	char resName[4], chain, altLoc;
	char resid[5];
	char x[9], y[9], z[9];
	char occupancy[7];
	char beta[7];
	strncpy(id, &line[6], 5);
	id[5] = '\0';
	strncpy(atomName, &line[12], 4);
	atomName[4] = '\0';
	altLoc = line[16];
	strncpy(resName, &line[17], 3);
	resName[3] = '\0';
	chain = line[21];
	strncpy(resid, &line[22], 4);
	resid[4] = '\0';
	strncpy(x, &line[30], 8);
	x[8] = '\0';
	strncpy(y, &line[38], 8);
	y[8] = '\0';
	strncpy(z, &line[46], 8);
	z[8] = '\0';
	strncpy(occupancy, &line[54], 6);
	occupancy[6] = '\0';
	strncpy(beta, &line[60], 6);
	beta[6] = '\0';
	strcpy(this->name, strtok(atomName, " "));
	this->altLoc = altLoc;
	this->name[4] = 0;
	this->chain = chain;
	this->resid = atoi(resid);
	strcpy(this->resName, resName);
	this->resName[3] = 0;
	this->id = atoi(id);
	this->x = atof(x);
	this->y = atof(y);
	this->z = atof(z);
	this->occupancy = atof(occupancy);
	this->beta = atof(beta);
#ifdef DEBUGPDBIO
	this->print();
#endif
}

/*
 * Parses single line of 'SSBOND' entry from pdb file.
 * SSBOND entry format in PDB:
 *
 * COLUMNS        DATA TYPE       FIELD         DEFINITION
 * -------------------------------------------------------------------
 *  1 -  6        Record name     "SSBOND"
 *  8 - 10        Integer         serNum       Serial number.
 * 12 - 14        LString(3)      "CYS"        Residue name.
 * 16             Character       chainID1     Chain identifier.
 * 18 - 21        Integer         seqNum1      Residue sequence number.
 * 22             AChar           icode1       Insertion code. (Ignored)
 * 26 - 28        LString(3)      "CYS"        Residue name.
 * 30             Character       chainID2     Chain identifier.
 * 32 - 35        Integer         seqNum2      Residue sequence number.
 * 36             AChar           icode2       Insertion code. (Ignored)
 * 60 - 65        SymOP           sym1         Symmetry oper for 1st resid (Ignored)
 * 67 - 72        SymOP           sym2         Symmetry oper for 2nd resid (Ignored)
 *
 * For more details visit http://www.wwpdb.org/documentation/format23/sect6.html
 *
 * Method parameters:
 * 		line: line from the file to parse
 */
void PDBSSBond::parse(const char* line){
	char chain;
	char resid[5];
	chain = line[15];
    resid[4] = '\0';
	strncpy(resid, &line[17], 4);
	this->chain1 = chain;
	this->resid1 = atoi(resid);
	chain = line[29];
	strncpy(resid, &line[31], 4);
	this->chain2 = chain;
	this->resid2 = atoi(resid);
}


void savePDB(const char* filename, const SOP& sop){
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



/*
 * Saves data from PDB object into a PDB file format
 *
 * For detailed information about PDB file format, please visit
 * http://www.wwpdb.org/documentation/format23/sect6.html
 *
 * Parameters:
 * 		filename: name of a file to save into (will be overwritten/created)
 */

void PDB::write(const char* filename, bool printConnections) const {
	printf("Saving PDB '%s'...\n", filename);
	FILE* file = safe_fopen(filename, "w");
	int i, j;
	for(i = 0; i < this->ssCount; i++){
		fprintf(file, "SSBOND %3d CYS %c %4d    CYS %c %4d\n",
								i + 1,
								this->ssbonds[i].chain1,
								this->ssbonds[i].resid1,
								this->ssbonds[i].chain2,
								this->ssbonds[i].resid2);
	}
	if(this->atomCount < 100000){
		for(i = 0; i < this->atomCount; i++){
			fprintf(file, "ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
									i + 1,
									this->atoms[i].name,
									this->atoms[i].altLoc,
									this->atoms[i].resName,
									this->atoms[i].chain,
									this->atoms[i].resid,
									this->atoms[i].x,
									this->atoms[i].y,
									this->atoms[i].z,
									this->atoms[i].occupancy,
									this->atoms[i].beta);
		}
	} else {
		for(i = 0; i < this->atomCount; i++){
			fprintf(file, "ATOM  %5x %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
									i + 1,
									this->atoms[i].name,
									this->atoms[i].altLoc,
									this->atoms[i].resName,
									this->atoms[i].chain,
									this->atoms[i].resid,
									this->atoms[i].x,
									this->atoms[i].y,
									this->atoms[i].z,
									this->atoms[i].occupancy,
									this->atoms[i].beta);
		}
	}

	if(printConnections){
		for(i = 0; i < this->atomCount; i++){
			fprintf(file, "CONECT%5d", i+1);
			for(j = 0; j < this->connections.connectCount[i]; j++){
				fprintf(file, "%5d", this->connections.connectMap[j*this->atomCount + i]+1);
			}
			fprintf(file, "\n");
		}
	}
	fprintf(file, "END");
	fclose(file);
	printf("Done saving PDB.\n");
}

/*
 * Prints PDBAtom object in a PDB format.
 * Parameters:
 * 		atomData: PDBAtom object to print
 */
void PDBAtom::print() const {
    // TODO: use PDBAtom::fprint
	printf("ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
			id,
			name,
			altLoc,
			resName,
			chain,
			resid,
			x,
			y,
			z,
			occupancy,
			beta);
}

/*
 * Prints PDBAtom object in a PDB format into a file.
 * Parameters:
 * 		atomData: PDBAtom object to print
 */
void PDBAtom::fprint(FILE* file) const {
	fprintf(file, "ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %c\n",
			id,
			name,
			altLoc,
			resName,
			' ',
			resid,
			x,
			y,
			z,
			occupancy,
			beta,
			chain);
}

void readCoord(const char* filename, SOP& sop){
	//printf("Reading coordinates from '%s'.\n", filename);
	char buffer[BUF_SIZE+1];
	FILE* file = fopen(filename, "r");
	if(file == NULL){
		printf("ERROR: Coordinates file %s can not be found.\n", filename);
		exit(-1);
	}
	int index = 0;
	while(fgets(buffer, BUF_SIZE, file) != NULL){
		//char* pch = strtok(buffer, " ");
		if(strncmp(buffer, "ATOM", 4) == 0){
			char x[9], y[9], z[9];
			strncpy(x, &buffer[30], 9);
			strncpy(y, &buffer[38], 9);
			strncpy(z, &buffer[46], 9);
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
		exit(-1);
	}
	fclose(file);
	//printf("Done reading coordinates.\n");
}

