#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "structure.h"
#include "config_reader.h"

#define version "1.07"
#define MODE_CAPSID 1

//#include "sop.h"

//#define DEBUG
//#define DEBUG1

extern int BLOCK_SIZE;

extern int device;
extern int mode;

extern int Ntr;
extern int firstrun;

extern long int numsteps;
extern int DCDfreq;
extern int restartfreq;

extern float R_sigma;
extern float eh;
extern int nav;

extern int pairs_freq;
extern int possiblepairs_freq;

extern int restart;

extern int stepnum;

extern int pullingOn;
extern int heatingOn;
extern int indentationOn;
extern int minimizationOn;

/*
 * Model parameters
 */

extern long long int initialTime;
extern long long int lastTime;
extern long long int switchTime;

extern int readInitialCoord;

extern char top_filename[100];
extern char coord_filename[100];
extern char ref_filename[100];
extern char dcd_filename[100];
extern char dat_filename[100];
extern char restartpdb_filename[100];
extern char restartconf_filename[100];
extern char final_filename[100];

extern SSBond *SSbonds;
extern PDB pdbdata;
extern SOP sop;

extern int    seed, run;
extern long int step;
extern float x_R, R, rg;

extern float forcex, forcey, forcez, xt, x0;
extern float epot_fene, epot_LJ, tempav, epot_LJ_att, epot_LJ_nei, epot_LJ_rep;

extern int* pdbRef;
