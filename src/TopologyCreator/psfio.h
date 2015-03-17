/*
 * psfio.h
 *
 *  Created on: Mar 19, 2009
 *      Author: rabab
 */

#ifndef PSFIO_H_
#define PSFIO_H_

#define PSF_NAME_SIZE 5
#define PSF_TYPE_SIZE 5
#define PSF_SEGMENT_SIZE 5
#define PSF_RESNAME_SIZE 4

/*
 * Structures
 */

typedef struct {
  int id;
  char name[PSF_NAME_SIZE+1];
  char type[PSF_TYPE_SIZE+1];
  char segment[PSF_SEGMENT_SIZE+1];
  char resName[PSF_RESNAME_SIZE+1];
  int resid;
  float m;
  float q;
} PSFAtom;

typedef struct {
	int i;
	int j;
} PSFBond;

typedef struct {
	int i;
	int j;
	int k;
} PSFAngle;

typedef struct {
	int i;
	int j;
	int k;
	int l;
} PSFDihedral;

typedef struct {
	int i;
	int j;
	int k;
	int l;
} PSFImproper;

typedef struct {
	int i1;
	int j1;
	int k1;
	int l1;
	int i2;
	int j2;
	int k2;
	int l2;
} PSFCMAP;

typedef struct {
	int natom;
	int nbond;
	int ntheta;
	int nphi;
	int nimphi;
	int ncmap;
	int nnb;
	PSFAtom* atoms;
	PSFBond* bonds;
	PSFAngle* angles;
	PSFDihedral* dihedrals;
	PSFImproper* impropers;
	PSFCMAP* cmaps;
	int* nbExclusions;
	int* nbExclusionsCounts;
} PSF;

/*
 * Public methods
 */

void readPSF(const char* filename, PSF* psfdata);
void writePSF(const char* filename, PSF* psfdata);

#endif
