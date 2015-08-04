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
	this->is_xyz = (!strcmp(filename + strlen(filename)-4, ".xyz"));
	if (this->is_xyz) {
		this->readXYZ(filename);
		return;
	}
	int ss_count = 0, atoms_count = 0;

	char buffer[BUF_SIZE];
	FILE* file = safe_fopen(filename, "r");

	while(fgets(buffer, BUF_SIZE, file) != NULL){
		if(strncmp(buffer,"SSBOND",6) == 0){
            this->ssbonds.push_back(PDBSSBond(buffer));
		}
		if(strncmp(buffer, "ATOM", 4) == 0){
            this->atoms.push_back(PDBAtom(buffer));
		}
	}

	printf("Found %zu atoms and %zu SS-bonds.\n", this->atoms.size(), this->ssbonds.size());
	printf("Done reading '%s'.\n", filename);
	fclose(file);
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
	this->id = atoi(id);

	strncpy(atomName, &line[12], 4);
	atomName[4] = '\0';
	strcpy(this->name, strtok(atomName, " "));
	this->name[4] = '\0';

	altLoc = line[16];
	this->altLoc = altLoc;

	strncpy(resName, &line[17], 3);
	resName[3] = '\0';
	strcpy(this->resName, resName);
	this->resName[3] = '\0';

	chain = line[21];
	this->chain = chain;

	strncpy(resid, &line[22], 4);
	resid[4] = '\0';
	this->resid = atoi(resid);

	strncpy(x, &line[30], 8);
	x[8] = '\0';
	this->x = atof(x);

	strncpy(y, &line[38], 8);
	y[8] = '\0';
	this->y = atof(y);

	strncpy(z, &line[46], 8);
	z[8] = '\0';
	this->z = atof(z);

	strncpy(occupancy, &line[54], 6);
	occupancy[6] = '\0';
	this->occupancy = atof(occupancy);

	strncpy(beta, &line[60], 6);
	beta[6] = '\0';
	this->beta = atof(beta);
#ifdef DEBUGPDBIO
	this->fprint(stdout);
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
	this->chain1 = chain;

	strncpy(resid, &line[17], 4);
    resid[4] = '\0';
	this->resid1 = atoi(resid);

	chain = line[29];
	this->chain2 = chain;

	strncpy(resid, &line[31], 4);
    resid[4] = '\0';
	this->resid2 = atoi(resid);
}


void savePDB(const char* filename, const SOP& sop){
	FILE* file = safe_fopen(filename, "w");
	int i;

    bool need_hex_atomid = (sop.aminos.size() + sop.additionalAminos.size() >= 100000);
	for(i = 0; i < sop.aminos.size(); i++){
        sop.aminos[i].fprint(file, i+1, need_hex_atomid);
	}
	for(i = 0; i < sop.additionalAminos.size(); i++){
        sop.additionalAminos[i].fprint(file, sop.aminos.size() + i + 1, need_hex_atomid);
	}

	int lastAmino = -1;
	for(int i = 0; i < sop.bonds.size(); i++){
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
	for(i = 0; i < this->ssbonds.size(); i++){
        this->ssbonds[i].fprint(file, i+1);
	}

	bool need_hex_atomid = (this->atoms.size() >= 100000);
	for(i = 0; i < this->atoms.size(); i++){
        this->atoms[i].fprint(file, i+1, need_hex_atomid);
    }

	if(printConnections){
		for(i = 0; i < this->atoms.size(); i++){
			fprintf(file, "CONECT%5d", i+1);
			for(j = 0; j < this->connections.connectCount[i]; j++){
				fprintf(file, "%5d", this->connections.connectMap[j*this->atoms.size() + i]+1);
			}
			fprintf(file, "\n");
		}
	}
	fprintf(file, "END");
	fclose(file);
	printf("Done saving PDB.\n");
}

/*
 * Prints PDBAtom object in a PDB format into a file.
 * Parameters:
 * 		atomData: PDBAtom object to print
 */
void PDBAtom::fprint(FILE* file, int custom_atomid, bool atomid_in_hex) const {
    char format[] = "ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
    if (atomid_in_hex) format[8] = 'x';
	fprintf(file, format,
			(custom_atomid < 0 ? id : custom_atomid),
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

void PDBSSBond::fprint(FILE* file, int bondid) const{
	fprintf(file, "SSBOND %3d CYS %c %4d    CYS %c %4d\n",
							bondid,
							this->chain1,
							this->resid1,
							this->chain2,
							this->resid2);
}


void readCoord(const char* filename, SOP& sop){
	//printf("Reading coordinates from '%s'.\n", filename);
	int is_xyz = 0;
	char buffer[BUF_SIZE+1];
	FILE* file = safe_fopen(filename, "r");
	if(file == NULL){
		DIE("ERROR: Coordinates file %s can not be found.\n", filename);
	}
	if (!strcmp(filename + strlen(filename)-4, ".xyz")) {
		is_xyz = 1;
		safe_fgets(buffer, BUF_SIZE, file);
		safe_fgets(buffer, BUF_SIZE, file);
	}
	int index = 0;
	while(fgets(buffer, BUF_SIZE, file) != NULL){
		if(is_xyz){
			char *pch = strtok(buffer, " \t\r\n");
			pch = strtok(NULL, " \t\r\n");
			sop.aminos[index].x = atof(pch);
			pch = strtok(NULL, " \t\r\n");
			sop.aminos[index].y = atof(pch);
			pch = strtok(NULL, " \t\r\n");
			sop.aminos[index].z = atof(pch);
			index++;
		}
		if(!is_xyz && strncmp(buffer, "ATOM", 4) == 0){
            PDBAtom tmp;
            tmp.parse(buffer);
			sop.aminos[index].x = tmp.x;
			sop.aminos[index].y = tmp.y;
			sop.aminos[index].z = tmp.z;
#ifdef DEBUG
			printf("%d: %f, %f, %f\n", index, sop.aminos[index].x, sop.aminos[index].y, sop.aminos[index].z);
#endif
			index++;
		}
	}
	fclose(file);
	if(index == 0){
		DIE("Can't read pdb file.\n");
	}
	if(index != sop.aminos.size()){
		DIE("Read coordinates for %d beads, yet topology has %zu beads.\n", index, sop.aminos.size());
	}
	//printf("Done reading coordinates.\n");
}

void PDB::readXYZ(const char* filename){
	printf("Reading %s.\n", filename);
	char buffer[BUF_SIZE];
	FILE* file = safe_fopen(filename, "r");
	safe_fgets(buffer, BUF_SIZE, file);
	const int n = atoi(buffer);
	atoms.resize(n);
	safe_fgets(buffer, BUF_SIZE, file); // Skip mandatory comment line
	int i;
	char* pch;
	for(i = 0; i < n; i++){
		safe_fgets(buffer, BUF_SIZE, file);
		pch = strtok(buffer, " \t\r\n");
		strncpy(atoms[i].name,pch,4);
		pch = strtok(NULL, " \t\r\n");
		atoms[i].x = atof(pch);
		pch = strtok(NULL, " \t\r\n");
		atoms[i].y = atof(pch);
		pch = strtok(NULL, " \t\r\n");
		atoms[i].z = atof(pch);
		atoms[i].id = atoms[i].resid = i+1;
		atoms[i].chain = 'X';
		strcpy(atoms[i].resName, "XXX");
		atoms[i].altLoc = 0;
		atoms[i].occupancy = atoms[i].beta = 0.0;
	}
	printf("Done reading '%s'.\n", filename);
	fclose(file);
}

