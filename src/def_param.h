#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "structure.h"
#include "IO/configreader.h"
#include "IO/topio.h"
#include "IO/pdbio.h"

#define version "trunk"

const int MODE_CAPSID = 1;

//#include "sop.h"

//#define DEBUG
//#define DEBUG1

extern int BLOCK_SIZE;

extern int mode;

extern long int numsteps;

extern int nav;

/*
 * Model parameters
 */

extern long long int initialTime;
extern long long int lastTime;

extern char top_filename[100];
extern char coord_filename[100];
extern char ref_filename[100];
extern char dcd_filename[100];
extern char dat_filename[100];
extern char restartpdb_filename[100];
extern char restartconf_filename[100];
extern char final_filename[100];

extern PDB pdbdata;
extern SOP sop;

extern long int step;

