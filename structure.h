#ifndef STRUCTURE_H_
#define STRUCTURE_H_


struct Atom{

  int id;
  char   name[5], chain, resName[4];
  int    resid;
  float x, y, z;

  float occupancy;
  float beta;

  int ssbond_aminoid; //TO BE REMOVED

};

struct SSBond{

	int resid1;
	char chain1;

	int resid2;
	char chain2;

};

struct PDB{
	int atomCount;
	Atom* atoms;
	int ssCount;
	SSBond* ssbonds;
};

struct CovalentBond{
	int i;
	int j;
	float r0;
};

struct NativeContact{
	int i;
	int j;
	float r0;
	float eh;
};

struct PossiblePair{
	int i;
	int j;
};

struct SOP{
	int aminoCount;
	int bondCount;
	int nativeCount;
	int pairCount;
	Atom* aminos;
	CovalentBond* bonds;
	NativeContact* natives;
	PossiblePair* pairs;

	int additionalAminosCount;
	Atom* additionalAminos;
};

#endif
