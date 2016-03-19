/*
 * main.cpp
 *
 *  Created on: Nov 28, 2012
 *      Author: kias
 */

#include "TopologyCreator/aatocg.h"
#include "TopologyCreator/pdbio.h"
#include "TopologyCreator/topio.h"
#include "TopologyCreator/psfio.h"
#include "IO/configreader.h"
#include "Util/parameters.h"

#define BUF_SIZE 1024

PARAMETER_MANDATORY(cgconfig, std::string, "path", "Path to the file with information on how to coarse-grain the system.")
PARAMETER_MANDATORY(structure, std::string, "path", "Path to the file with initial (all-atom) coordinates.")
PARAMETER(additional_bonds, std::string, "NONE", "path", "Description of additional bonds (i.e. S-S, crosslinks), if needed.")

PARAMETER_MANDATORY(coordinates, std::string, "path", "Path to save coarse-grained coordinates PDB (will be used by SOP-GPU).")
PARAMETER_MANDATORY(topology, std::string, "path", "Path to save topology in Gromacs format (will be used by SOP-GPU).")

PARAMETER(topology_psf, std::string, "NONE", "path", "Path to save topology inn PSF format (convinient for VMD).")
PARAMETER(topology_natpsf, std::string, "NONE", "path", "Path to save native contact topology in PSF format (for visual inspection of contacts).")

PARAMETER(R_limit_bond, float, 8.0f, "A", "Cutoff distance for Ca atoms to assume native contact")
PARAMETER(SC_limit_bond, float, 5.0f, "A", "Cutoff distance for side-chain atoms to assume native contact")
PARAMETER_MANDATORY(eh, std::string, "value/key", "Either value for eh or the 'O' ('B') key to take the value from occupancy(beta) column.")
PARAMETER(Ks, float, 20.0f, "kcal/mol", "Covalent spring constant")

PARAMETER(use_chains, bool, false, "true/false", "If YES, the structure will be split into chains according to the 'chain' column of the initial PDB file. By default done using segment. Chains are preferable for small systems.")

CGConfig conf;
PDB pdb;

struct SOPBead {
	char name[5], type[5], resname[5], segment[5];
	int resid;
	double x, y, z;
	double m, q;
	double occupancy, beta;
	std::vector<Atom> connectedTo;
	std::vector<PDBAtom> represents;
};

typedef struct {
	int i;
	int j;
	float Kb;
	float b0;
} SOPBond;

typedef struct {
	int i;
	int j;
	float eh;
	float r0;
} SOPNative;

typedef struct {
	int i;
	int j;
} SOPPair;

typedef struct {
	char res_name[5];
	int resid;
	std::vector<PDBAtom*> atoms;
}  PDBTreeResidue;

typedef struct {
	char segment_name[5];
	std::vector<PDBTreeResidue> residues;
} PDBTreeSegment;

std::vector<PDBTreeSegment> pdb_tree;
std::vector<SOPBead> beads;


bool bond_comparator(SOPBond b1, SOPBond b2){
	if(b1.i < b2.i){
		return true;
	} else
	if(b1.i == b2.i){
		return b1.j < b2.j;
	} else {
		return false;
	}
}

bool native_comparator(SOPNative n1, SOPNative n2){
	if(n1.i < n2.i){
		return true;
	} else
	if(n1.i == n2.i){
		return n1.j < n2.j;
	} else {
		return false;
	}
}

bool pairs_comparator(SOPPair p1, SOPPair p2){
	if(p1.i < p2.i){
		return true;
	} else
	if(p1.i == p2.i){
		return p1.j < p2.j;
	} else {
		return false;
	}
}

std::vector<SOPBond> bonds;
std::vector<SOPNative> natives;

std::vector<SOPPair> bonds13;

int** bondsList;
int* bondCount;

void createBondsList(int N){
	bondsList = (int**)calloc(N, sizeof(int*));
	bondCount = (int*)calloc(N, sizeof(int));
	int i;
	for(i = 0; i < N; i++){
		bondsList[i] = (int*)calloc(10, sizeof(int));
	}
	for(i = 0; i < bonds.size(); i++){
		SOPBond bond = bonds.at(i);
		bondsList[bond.i][bondCount[bond.i]] = bond.j;
		bondsList[bond.j][bondCount[bond.j]] = bond.i;
		bondCount[bond.i] ++;
		bondCount[bond.j] ++;
	}
}

int** bonds13List;
int* bond13Count;

void createBonds13List(int N){
	bonds13List = (int**)calloc(N, sizeof(int*));
	bond13Count = (int*)calloc(N, sizeof(int));
	int i;
	for(i = 0; i < N; i++){
		bonds13List[i] = (int*)calloc(100, sizeof(int));
	}
	for(i = 0; i < bonds13.size(); i++){
		SOPPair bond = bonds13.at(i);
		bonds13List[bond.i][bond13Count[bond.i]] = bond.j;
		bonds13List[bond.j][bond13Count[bond.j]] = bond.i;
		bond13Count[bond.i] ++;
		bond13Count[bond.j] ++;
	}
}

void addConnection(char* segment1, int resid1, char* bead1,
		char* segment2, int resid2, char* bead2);
unsigned findBead(const char* segment, int resid, const char* bead_name);
void addConnection(SOPBead bead, Atom conn, unsigned int j);

void addBond(int i, int j, float Kb, float b0);
void addNative(int i, int j, float eh, float r0);
int checkBonded(int i, int j);
int check13Bonded(int i, int j);
int checkNatived(int i, int j);
double getDistance(SOPBead bead1, SOPBead bead2);
double getDistance(PDBAtom atom1, PDBAtom atom2);


using namespace configreader;

int main(int argc, char* argv[]){

	printf("==========================\n");
	printf("SOP-GPU Topology creator 2.0\n");
	printf("==========================\n");

	if(argc < 1){
		printf("ERROR: Configuration file should be specified.\n");
		exit(-1);
	}

	parseParametersFile(argv[1]); //, argc, argv);

	readCGConfiguration(parameters::cgconfig.get().c_str(), &conf);

	readPDB(parameters::structure.get().c_str(), &pdb);


	printf("Building PDB tree...\n");

	int i, b;
	unsigned int j, k, l, m, n, o;

	if(parameters::use_chains.get()){
		printf("Using chain entries to separate polypeptide chains.\n");
		for(i = 0; i < pdb.atomCount; i++){
			sprintf(pdb.atoms[i].segment, "%c", pdb.atoms[i].chain);
		}
	}
	printf("Using segment entries to separate polypeptide chains.\n");

	for(i = 0; i < pdb.atomCount; i++){
		bool found = false;
		for(j = 0; j < pdb_tree.size(); j++){
			if(strcmp(pdb_tree.at(j).segment_name, pdb.atoms[i].segment) == 0){
				found = true;
			}
		}
		if(!found){
			PDBTreeSegment segment;
			sprintf(segment.segment_name, "%s", pdb.atoms[i].segment);
			pdb_tree.push_back(segment);
		}
	}
	printf("Found %ld segments.\n", pdb_tree.size());
	for(i = 0; i < pdb.atomCount; i++){
		for(j = 0; j < pdb_tree.size(); j++){
			if(strcmp(pdb_tree.at(j).segment_name, pdb.atoms[i].segment) == 0){
				bool found = false;
				for(k = 0; k < pdb_tree.at(j).residues.size(); k++){
					if(pdb_tree.at(j).residues.at(k).resid == pdb.atoms[i].resid){
						found = true;
					}
				}
				if(!found){
					PDBTreeResidue residue;
					sprintf(residue.res_name, "%s", pdb.atoms[i].resName);
					residue.resid = pdb.atoms[i].resid;
					pdb_tree.at(j).residues.push_back(residue);
				}
			}
		}
	}
	for(i = 0; i < pdb.atomCount; i++){
		for(j = 0; j < pdb_tree.size(); j++){
			if(strcmp(pdb_tree.at(j).segment_name, pdb.atoms[i].segment) == 0){
				for(k = 0; k < pdb_tree.at(j).residues.size(); k++){
					if(pdb_tree.at(j).residues.at(k).resid == pdb.atoms[i].resid){
						pdb_tree.at(j).residues.at(k).atoms.push_back(&pdb.atoms[i]);
					}
				}
			}
		}
	}

	for(j = 0; j < pdb_tree.size(); j++){
		printf("\n\nSegment %s:\n", pdb_tree.at(j).segment_name);
		for(k = 0; k < pdb_tree.at(j).residues.size(); k++){
			printf("\nResid %s%d:\n",
					pdb_tree.at(j).residues.at(k).res_name, pdb_tree.at(j).residues.at(k).resid);
			for(l = 0; l < pdb_tree.at(j).residues.at(k).atoms.size(); l++){
				printAtom(*pdb_tree.at(j).residues.at(k).atoms.at(l));
			}
		}
	}

	printf("PDB tree completed.\n");

	printf("Creating coarse-grained beads.\n");

	for(j = 0; j < pdb_tree.size(); j++){
		for(k = 0; k < pdb_tree.at(j).residues.size(); k++){
			Residue resid;
			PDBTreeResidue pdb_tree_resid;
			for(m = 0; m < conf.residues.size(); m++){
				if(strcmp(conf.residues.at(m).resname, pdb_tree.at(j).residues.at(k).res_name) == 0){
					resid = conf.residues.at(m);
					pdb_tree_resid = pdb_tree.at(j).residues.at(k);
				}
			}
			for(n = 0; n < resid.beads.size(); n++){
				SOPBead bead;
				bead.m = resid.beads.at(n).mass;
				bead.q = resid.beads.at(n).charge;
				sprintf(bead.name, "%s", resid.beads.at(n).name);
				sprintf(bead.type, "%s", resid.beads.at(n).type);
				sprintf(bead.resname, "%s", pdb_tree_resid.res_name);
				sprintf(bead.segment, "%s", pdb_tree.at(j).segment_name);
				bead.resid = pdb_tree_resid.resid;
				bead.x = 0.0;
				bead.y = 0.0;
				bead.z = 0.0;
				bead.connectedTo = resid.beads.at(n).connectedTo;
				int count = 0;
				for(l = 0; l < pdb_tree_resid.atoms.size(); l++){
					for(o = 0; o < resid.beads.at(n).atomsCM.size(); o++){
						if(strcmp(pdb_tree_resid.atoms.at(l)->name,
								resid.beads.at(n).atomsCM.at(o).name) == 0){
							bead.x += pdb_tree_resid.atoms.at(l)->x;
							bead.y += pdb_tree_resid.atoms.at(l)->y;
							bead.z += pdb_tree_resid.atoms.at(l)->z;
							count ++;
						}
					}
				}
				if(count != 0){
					bead.x /= (double)count;
					bead.y /= (double)count;
					bead.z /= (double)count;
				}
				count = 0;
				bead.beta = 0.0;
				bead.occupancy = 0.0;
				for(l = 0; l < pdb_tree_resid.atoms.size(); l++){
					for(o = 0; o < resid.beads.at(n).atomsRepresents.size(); o++){
						if(strcmp(pdb_tree_resid.atoms.at(l)->name,
								resid.beads.at(n).atomsRepresents.at(o).name) == 0){
							//PDBAtom atom;
							//memcpy(&atom, &pdb_tree_resid.atoms.at(l), sizeof(PDBAtom));
							bead.represents.push_back(*pdb_tree_resid.atoms.at(l));
							bead.beta += pdb_tree_resid.atoms.at(l)->beta;
							bead.occupancy += pdb_tree_resid.atoms.at(l)->occupancy;
							count ++;
						}
					}
				}
				bead.beta /= count;
				bead.occupancy /= count;
				beads.push_back(bead);
			}
		}
	}

	printf("Coarse-grained beads are added.\n");

	printf("Adding covalent bonds...\n");

	for(j = 0; j < beads.size(); j++){
		SOPBead bead = beads.at(j);
		for(k = 0; k < bead.connectedTo.size(); k++){
			Atom conn = bead.connectedTo.at(k);
			addConnection(bead, conn, j);
		}
	}

	if(strncmp(parameters::additional_bonds.get().c_str(), "NONE", 4) != 0){
		FILE* file = fopen(parameters::additional_bonds.get().c_str(), "r");
		char* segment1;
		int resid1;
		char* bead1;
		char* segment2;
		int resid2;
		char* bead2;
		if(file != NULL){
			char* pch;
			char buffer[BUF_SIZE];
			while(fgets(buffer, BUF_SIZE, file) != NULL){
				if(strncmp(buffer, "CONN", 4) == 0){
					pch = strtok(buffer, " \t\n\r");
					pch = strtok(NULL, " \t\n\r");
					segment1 = pch;
					pch = strtok(NULL, " \t\n\r");
					resid1 = atoi(pch);
					pch = strtok(NULL, " \t\n\r");
					bead1 = pch;
					pch = strtok(NULL, " \t\n\r");
					segment2 = pch;
					pch = strtok(NULL, " \t\n\r");
					resid2 = atoi(pch);
					pch = strtok(NULL, " \t\n\r");
					bead2 = pch;
					addConnection(segment1, resid1, bead1, segment2, resid2, bead2);
				}
			}
		}
		fclose(file);
	} else {
		printf("No additional covalent bonds (S-S, crosslinks, etc.) have been specified.\n");
	}

	std::sort(bonds.begin(), bonds.end(), bond_comparator);

	createBondsList(beads.size());

	int b1, b2;
	for(i = 0; i < beads.size(); i++){
		for(b1 = 0; b1 < bondCount[i]; b1++){
			j = bondsList[i][b1];
			for(b2 = 0; b2 < bondCount[j]; b2++){
				k = bondsList[j][b2];
				if(i < k){
					SOPPair bond13;
					bond13.i = i;
					bond13.j = k;
					bonds13.push_back(bond13);
					//printf("%d-%d-%d\n", i, j, k);
				}
			}
		}
	}

	std::sort(bonds13.begin(), bonds13.end(), pairs_comparator);

	createBonds13List(beads.size());

	printf("Covalent bonds are added.\n");

	printf("Adding native contacts...\n");

	double eh, cutoff, cutoffAtomistic;
	std::string ehstring = parameters::eh.get();
	cutoff = parameters::R_limit_bond.get();
	cutoffAtomistic = parameters::SC_limit_bond.get();
	for(j = 0; j < beads.size(); j++){
		if(j % 100 == 0){
			printf("Bead %d out of %ld\n", j, beads.size());
		}
		for(k = j + 1; k < beads.size(); k++){
			double r0 = getDistance(beads.at(j), beads.at(k));

			if(ehstring.compare("O") == 0){
				eh = sqrt(beads.at(j).occupancy*beads.at(k).occupancy);
			} else if(ehstring.compare("B") == 0){
				eh = sqrt(beads.at(j).beta*beads.at(k).beta);
			} else if(atof(ehstring.c_str()) != 0){
				eh = atof(ehstring.c_str());
				//printf("%f\n", eh);
			} else {
				exit(0);
			}
			bool added = false;
			if(r0 < cutoff){
				if((!checkBonded(j, k)) && (!check13Bonded(j, k))){
					addNative(j, k, eh, r0);
					added = true;
				}
			}
			if(!added){
				for(l = 0; l < beads.at(j).represents.size(); l++){
					for(m = 0; m < beads.at(k).represents.size(); m++){
						double r1 = getDistance(beads.at(j).represents.at(l), beads.at(k).represents.at(m));
						if((!added) && r1 < cutoffAtomistic){
							if((!checkBonded(j, k)) && (!check13Bonded(j, k))){
								addNative(j, k, eh, r0);
								added = true;
							}
						}
					}
				}
			}
		}
	}

	printf("Native contacts are added.\n");

	std::sort(natives.begin(), natives.end(), native_comparator);

	double pot = 0.0;
	for(j = 0; j < natives.size(); j++){
		SOPNative native = natives.at(j);
		pot += native.eh;
	}
	printf("Total energy: %f\n", pot);

	printf("Saving PDB, PSF and TOP files for coarse-grained system...\n");

	PDB cgpdb;
	cgpdb.atomCount = beads.size();
	cgpdb.ssCount = 0;
	cgpdb.atoms = (PDBAtom*)calloc(cgpdb.atomCount, sizeof(PDBAtom));
	for(j = 0; j < beads.size(); j++){
		cgpdb.atoms[j].altLoc = ' ';
		cgpdb.atoms[j].beta = beads.at(j).m;
		cgpdb.atoms[j].chain = beads.at(j).segment[0];
		cgpdb.atoms[j].id = j+1;
		sprintf(cgpdb.atoms[j].name, "%s", beads.at(j).name);
		cgpdb.atoms[j].occupancy = beads.at(j).q;
		sprintf(cgpdb.atoms[j].resName, "%s", beads.at(j).resname);
		sprintf(cgpdb.atoms[j].segment, "%s", beads.at(j).segment);
		cgpdb.atoms[j].resid = beads.at(j).resid;
		cgpdb.atoms[j].x = beads.at(j).x;
		cgpdb.atoms[j].y = beads.at(j).y;
		cgpdb.atoms[j].z = beads.at(j).z;
	}
	writePDB(parameters::coordinates.get().c_str(), &cgpdb);


	if(strncmp(parameters::topology_psf.get().c_str(), "NONE", 4) != 0){
		PSF cgpsf;
		cgpsf.natom = cgpdb.atomCount;
		cgpsf.atoms = (PSFAtom*)calloc(cgpsf.natom, sizeof(PSFAtom));
		cgpsf.nbond = 0;
		cgpsf.nbond = bonds.size();
		cgpsf.bonds = (PSFBond*)calloc(cgpsf.nbond, sizeof(PSFBond));
		cgpsf.ncmap = 0;
		cgpsf.nimphi = 0;
		cgpsf.nnb = 0;
		cgpsf.nphi = 0;
		cgpsf.ntheta = 0;
		for(i = 0; i < cgpdb.atomCount; i++){
			cgpsf.atoms[i].id = i+1;
			sprintf(cgpsf.atoms[i].name, "%s", cgpdb.atoms[i].name);
			sprintf(cgpsf.atoms[i].type, "%s", cgpdb.atoms[i].name);
			sprintf(cgpsf.atoms[i].segment, "%s", cgpdb.atoms[i].segment);

			cgpsf.atoms[i].m = cgpdb.atoms[i].beta;
			cgpsf.atoms[i].q = cgpdb.atoms[i].occupancy;


			cgpsf.atoms[i].resid = cgpdb.atoms[i].resid;
			sprintf(cgpsf.atoms[i].resName, "%s", cgpdb.atoms[i].resName);
		}
		for(b = 0; b < (int)bonds.size(); b++){
			cgpsf.bonds[b].i = bonds.at(b).i + 1;
			cgpsf.bonds[b].j = bonds.at(b).j + 1;
		}
		writePSF(parameters::topology_psf.get().c_str(), &cgpsf);
	}

	if(strncmp(parameters::topology_natpsf.get().c_str(), "NONE", 4) != 0){
		PSF cgpsfnat;
		cgpsfnat.natom = cgpdb.atomCount;
		cgpsfnat.atoms = (PSFAtom*)calloc(cgpsfnat.natom, sizeof(PSFAtom));
		cgpsfnat.nbond = 0;
		cgpsfnat.nbond = natives.size();
		cgpsfnat.bonds = (PSFBond*)calloc(cgpsfnat.nbond, sizeof(PSFBond));
		cgpsfnat.ncmap = 0;
		cgpsfnat.nimphi = 0;
		cgpsfnat.nnb = 0;
		cgpsfnat.nphi = 0;
		cgpsfnat.ntheta = 0;
		for(i = 0; i < cgpdb.atomCount; i++){
			cgpsfnat.atoms[i].id = i+1;
			sprintf(cgpsfnat.atoms[i].name, "%s", cgpdb.atoms[i].name);
			sprintf(cgpsfnat.atoms[i].type, "%s", cgpdb.atoms[i].name);
			sprintf(cgpsfnat.atoms[i].segment, "%s", cgpdb.atoms[i].segment);

			cgpsfnat.atoms[i].m = cgpdb.atoms[i].beta;
			cgpsfnat.atoms[i].q = cgpdb.atoms[i].occupancy;


			cgpsfnat.atoms[i].resid = cgpdb.atoms[i].resid;
			sprintf(cgpsfnat.atoms[i].resName, "%s", cgpdb.atoms[i].resName);
		}
		for(b = 0; b < (int)natives.size(); b++){
			cgpsfnat.bonds[b].i = natives.at(b).i + 1;
			cgpsfnat.bonds[b].j = natives.at(b).j + 1;
		}
		writePSF(parameters::topology_natpsf.get().c_str(), &cgpsfnat);
	}

	TOPData cgtop;
	cgtop.atomCount = cgpdb.atomCount;
	cgtop.bondCount = bonds.size();
	cgtop.nativesCount = natives.size();
	cgtop.pairsCount = 0;//pairs.size();
	cgtop.angleCount = 0;
	cgtop.dihedralCount = 0;
	cgtop.atoms = (TOPAtom*)calloc(cgtop.atomCount, sizeof(TOPAtom));
	cgtop.bonds = (TOPPair*)calloc(cgtop.bondCount, sizeof(TOPPair));
	cgtop.natives = (TOPPair*)calloc(cgtop.nativesCount, sizeof(TOPPair));
	cgtop.pairs = (TOPPair*)calloc(cgtop.pairsCount, sizeof(TOPPair));
	for(i = 0; i < cgpdb.atomCount; i++){
		cgtop.atoms[i].id = i;
		sprintf(cgtop.atoms[i].name, "%s", cgpdb.atoms[i].name);
		cgtop.atoms[i].resid = cgpdb.atoms[i].resid;
		sprintf(cgtop.atoms[i].resName, "%s", cgpdb.atoms[i].resName);
		cgtop.atoms[i].chain = cgpdb.atoms[i].chain;
		sprintf(cgtop.atoms[i].type, "%s", cgpdb.atoms[i].name);
		cgtop.atoms[i].charge = beads.at(i).q;
		cgtop.atoms[i].mass = beads.at(i).m;
	}
	for(b = 0; b < (int)bonds.size(); b++){
		cgtop.bonds[b].i = bonds.at(b).i;
		cgtop.bonds[b].j = bonds.at(b).j;
		cgtop.bonds[b].func = 1;
		cgtop.bonds[b].c0 = bonds.at(b).b0;
		cgtop.bonds[b].c1 = 0.0f;
		cgtop.bonds[b].c2 = 0.0f;
		cgtop.bonds[b].c3 = 0.0f;
	}

	for(b = 0; b < (int)natives.size(); b++){
		cgtop.natives[b].i = natives.at(b).i;
		cgtop.natives[b].j = natives.at(b).j;
		cgtop.natives[b].func = 1;
		cgtop.natives[b].c0 = natives.at(b).r0;
		cgtop.natives[b].c1 = natives.at(b).eh;
		cgtop.natives[b].c2 = 0.0f;
		cgtop.natives[b].c3 = 0.0f;
	}

	writeTOP(parameters::topology.get().c_str(), &cgtop);

	printf("Files are completed.\n");

	return 0;
}

void addConnection(char* segment1, int resid1, char* bead1,
		char* segment2, int resid2, char* bead2){
	unsigned i = findBead(segment1, resid1, bead1);
	unsigned j = findBead(segment2, resid2, bead2);
	printf("Connecting: %s-%d-%s to %s-%d-%s\n",
			segment1, resid1, bead1, segment2, resid2, bead2);
	if(i < beads.size() && j < beads.size()){
		addBond(i, j, parameters::Ks.get(), getDistance(beads.at(i), beads.at(j)));
		printf("Added: %d-%d\n", i, j);
	}

}

unsigned findBead(const char* segment, int resid, const char* bead_name){
	unsigned i = 0;
	for(i = 0; i < beads.size(); i++){
		SOPBead bead = beads.at(i);
		if(strncmp(bead.segment, segment, 3) == 0 && bead.resid == resid && strncmp(bead.name, bead_name, 2) == 0){
			return i;
		}
	}
	return beads.size();
}

void addConnection(SOPBead bead, Atom conn, unsigned int j){
	unsigned int l;
	printf("%s %s to %s\n", bead.resname, bead.name, conn.name);
	if(conn.name[0] == '+'){
		l = j + 1;
		while(l < beads.size() && beads.at(l).resid != bead.resid + 1){
			l ++;
		}
		while(l < beads.size() && strcmp(beads.at(l).name, &conn.name[1]) != 0 &&
				beads.at(l).resid != bead.resid + 1){
			l ++;
		}
		if(l < beads.size() && strcmp(beads.at(j).segment, beads.at(l).segment) == 0){
			addBond(j, l, parameters::Ks.get(), getDistance(bead, beads.at(l)));
			printf("Added: + %d-%d\n", j, l);
		}
	} else
	if(conn.name[0] == '-'){
		l = j - 1;
		while(l > 0 && beads.at(l).resid != bead.resid - 1){
			l --;
		}
		while(l > 0 && strcmp(beads.at(l).name, &conn.name[1]) != 0){
			l --;
		}
		if(l > 0 && strcmp(beads.at(j).segment, beads.at(l).segment) == 0){
			addBond(j, l, parameters::Ks.get(), getDistance(bead, beads.at(l)));
			printf("Added: - %d-%d\n", j, l);
		}
	} else {
		l = j + 1;
		bool added = false;
		while(l < beads.size() && bead.resid == beads.at(l).resid &&
				 strcmp(beads.at(l).name, conn.name) != 0){
			l ++;
		}
		if(l < beads.size() && strcmp(beads.at(j).segment, beads.at(l).segment) == 0 &&
				beads.at(l).resid == bead.resid){
			addBond(j, l, parameters::Ks.get(), getDistance(bead, beads.at(l)));
			added = true;
			printf("Added: %d-%d\n", j, l);
		}
		if(!added){
			if(j - 1 >= 0){
				l = j - 1;
				while(l >= 0 && bead.resid == beads.at(l).resid &&
						strcmp(beads.at(l).name, conn.name) != 0){
					l --;
				}
				if(l >= 0 && strcmp(beads.at(j).segment, beads.at(l).segment) == 0 &&
						beads.at(l).resid == bead.resid){
					addBond(j, l, parameters::Ks.get(), getDistance(bead, beads.at(l)));
					printf("Added: %d-%d\n", j, l);
				}
			}
		}
	}
}

void addBond(int i, int j, float Kb, float b0){
	int b;
	int found = 0;
	for(b = 0; b < (int)bonds.size(); b++){
		if((i == bonds.at(b).i && j == bonds.at(b).j) ||
				(i == bonds.at(b).j && j == bonds.at(b).i)){
			found = 1;
		}
	}
	if(i != j && found != 1){
		SOPBond bond;
		bond.i = i;
		bond.j = j;
		bond.Kb = Kb;
		bond.b0 = b0;
		bonds.push_back(bond);
	}
}

void addNative(int i, int j, float eh, float r0){
	int b;
	int found = 0;
	for(b = 0; b < (int)natives.size(); b++){
		if((i == natives.at(b).i && j == natives.at(b).j) ||
				(i == natives.at(b).j && j == natives.at(b).i)){
			found = 1;
		}
	}
	if(i != j && found != 1){
		SOPNative native;
		native.i = i;
		native.j = j;
		native.eh = eh;
		native.r0 = r0;
		natives.push_back(native);
	}
}

int checkBonded(int i, int j){
	/*int b;
	for(b = 0; b < (int)bonds.size(); b++){
		if((i == bonds.at(b).i && j == bonds.at(b).j) ||
				(i == bonds.at(b).j && j == bonds.at(b).i)){
			return 1;
		}
	}
	return 0;*/
	int b;
	for(b = 0; b < bondCount[i]; b++){
		if(bondsList[i][b] == j){
			return 1;
		}
	}
	return 0;
}

int check13Bonded(int i, int j){
	int b;
	for(b = 0; b < bond13Count[i]; b++){
		if(bonds13List[i][b] == j){
			return 1;
		}
	}
	return 0;
}

int checkNatived(int i, int j){
	int b;
	for(b = 0; b < (int)natives.size(); b++){
		if((i == natives.at(b).i && j == natives.at(b).j) ||
				(i == natives.at(b).j && j == natives.at(b).i)){
			return 1;
		}
	}
	return 0;
}

double getDistance(SOPBead bead1, SOPBead bead2){
	double dx = bead1.x - bead2.x;
	double dy = bead1.y - bead2.y;
	double dz = bead1.z - bead2.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

double getDistance(PDBAtom atom1, PDBAtom atom2){
	double dx = atom1.x - atom2.x;
	double dy = atom1.y - atom2.y;
	double dz = atom1.z - atom2.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

