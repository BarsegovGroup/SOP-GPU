/*
 *  dcdio.c
 *
 * Code to work with dcd files. On Oct 30, 2008 only writing is implemented.
 *
 *  Created on: Oct 29, 2008
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#ifdef UNWRAP
# define safe_fopen fopen
# define safe_fgets fgets
# define safe_fread fread
# define DIE(format, ...) do{ printf(format, ##__VA_ARGS__); exit(-1); }while(0);
#else
# include "../Util/wrapper.h"
#endif
#include "dcdio.h"

/*
 * Extend the string to be 80 bytes length (for remarks in header)
 * Taken from dcdlib.C (NAMD source)
 */
void pad(char *s, int len);

/*
 * Open DCD file to write data into.
 * The FILE* pointer is returned and must be used for further calls
 */
void DCD::open_write(const char *dcd_filename){
	this->file = safe_fopen(dcd_filename, "w");
}

void DCD::open_append(const char *dcd_filename){
	this->file = safe_fopen(dcd_filename, "a");
}

void DCD::open_read(const char *dcd_filename){
	this->file = safe_fopen(dcd_filename, "r");
}

void DCD::close(){
	fclose(this->file);
}

void DCD::allocate(){
    this->X = (float*) malloc(N * sizeof(float));
    this->Y = (float*) malloc(N * sizeof(float));
    this->Z = (float*) malloc(N * sizeof(float));
}

void DCD::deallocate(){
    if(this->X) free(this->X);
    if(this->Y) free(this->Y);
    if(this->Z) free(this->Z);
}


/*
 * Write header to dcd file.
 * Inputs are:
 *  N - Number of atoms/beads in a system
 *	NFILE - Number of frames (will not be updated in current implementation)
 *	NPRIV - Starting timestep of DCD file - NOT ZERO
 *	NSAVC - Timesteps between DCD saves
 *	NSTEP - Number of timesteps
 *	DELTA - length of a timestep
 *
 *  DCD header format is the following:
 *
 *  Position	Length			data
 *  (bytes/4)	(bytes)
 *  ------------------------------------
 *  1			4				'84'
 *  2			4				'CORD'
 *  3			4			  	NFILE
 *  4			4				NPRIV
 *  5			4				NSAVC
 *  6			4				NPRIV-NSAVC or NSTEP
 *  7-11		5*4				Zeros
 *  12			4				DELTA
 *  13			4				With/without unit cell
 *  14-21		8				Unit cell description (or zeros)
 *  22			4				'24'
 *  23			4				'84'
 *  24			4				'164'
 *  25			4				'2'
 *  26-45		80				Remarks #1 ("REMARK CREATED BY NAMD FILENAME='...'")
 *  46-65		80				Remarks #2 ("REMARK" + date + username)
 *  66			4				'164'
 *  67			4				'4'
 *  68			4				N
 *  69			4				'4'
 *
 */
void DCD::write_header() const {
	int iout;
	float fout;
	char cout[5];
	iout = 84;
	fwrite(&iout, 4, 1, this->file);
	sprintf(cout, "CORD");
	fwrite(&cout, 4, 1, this->file);
	iout = NFILE;
	fwrite(&iout, 4, 1, this->file);
	iout = NPRIV;
	fwrite(&iout, 4, 1, this->file);
	iout = NSAVC;
	fwrite(&iout, 4, 1, this->file);
	iout = NPRIV-NSAVC;
	fwrite(&iout, 4, 1, this->file);
	iout = 0;
	for(int i = 0; i < 5; i++){
		fwrite(&iout, 4, 1, this->file);
	}
	fout = DELTA;
	fwrite(&fout, 4, 1, this->file);
	iout = 0;
	for(int i = 0; i < 9; i++){
		fwrite(&iout, 4, 1, this->file);
	}
	iout = 24;
	fwrite(&iout, 4, 1, this->file);
	iout = 84;
	fwrite(&iout, 4, 1, this->file);
	iout = 164;
	fwrite(&iout, 4, 1, this->file);
	iout = 2;
	fwrite(&iout, 4, 1, this->file);
	char title[81];
	sprintf(title, "REMARKS CREATED BY dcdio.cpp");
	pad(title, 80);
	fwrite(&title, 80, 1, this->file);
	time_t rawtime;
	struct tm *timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	sprintf(title, "REMARKS DATE: %s", asctime(timeinfo));
	pad(title, 80);
	fwrite(&title, 80, 1, this->file);
	iout = 164;
	fwrite(&iout, 4, 1, this->file);
	iout = 4;
	fwrite(&iout, 4, 1, this->file);
	iout = N;
	fwrite(&iout, 4, 1, this->file);
	iout = 4;
	fwrite(&iout, 4, 1, this->file);
}

/*
 * Writes frame into dcd
 * Input
 *  X, Y, Z - pointers to arrays with coordinates
 *
 *  Writing format is following
 *
 *  Length			Data
 *  (bytes)
 *  ---------------------------(If unit cell is defined - not implemented in this code)
 *  4				'48'
 *  12*4			Unit cell data
 *  4				'48'
 *  ------------------------(If unit cell is defined - not implemented in this code)
 *  4				N*4
 *  N*4				X
 *  4				N*4
 *  4				N*4
 *  N*4				Y
 *  4				N*4
 *  4				N*4
 *  N*4				Z
 *  4				N*4
 *
 */
void DCD::write_frame(float *X, float *Y, float *Z) const {
	int iout;
	iout = N*4;
	fwrite(&iout, 4, 1, this->file);
	fwrite(X, 4*N, 1, this->file);
	fwrite(&iout, 4, 1, this->file);

	fwrite(&iout, 4, 1, this->file);
	fwrite(Y, 4*N, 1, this->file);
	fwrite(&iout, 4, 1, this->file);

	fwrite(&iout, 4, 1, this->file);
	fwrite(Z, 4*N, 1, this->file);
	fwrite(&iout, 4, 1, this->file);
}

void DCD::read_header(){
	int iin;
	char cin[5];
	safe_fread(&iin, 4, 1, this->file);
	if(iin != 84){
		DIE("Error! Wrong DCD file: DCD supposed to have '84' in first 4 bytes, but it hasn't.");
	}
	safe_fread(&cin, 4, 1, this->file);
	char cord_char[5];
	sprintf(cord_char, "CORD");
	if(strncmp(cin, cord_char, 4)){
		DIE("Error! Wrong DCD file: no 'CORD' sequence at the beginning of the file. Found: %s.", cin);
	}

	safe_fread(&NFILE, 4, 1, this->file);
	safe_fread(&NPRIV, 4, 1, this->file);
	safe_fread(&NSAVC, 4, 1, this->file);
	safe_fread(&iin, 4, 1, this->file);
	for(int i = 0; i < 5; i++){
		safe_fread(&iin, 4, 1, this->file);
	}
	safe_fread(&DELTA, 4, 1, this->file);
	for(int i = 0; i < 9; i++){
		safe_fread(&iin, 4, 1, this->file);
	}
	//24
	safe_fread(&iin, 4, 1, this->file);
	//84;
	safe_fread(&iin, 4, 1, this->file);
	//164;
	safe_fread(&iin, 4, 1, this->file);
	//2;
	safe_fread(&iin, 4, 1, this->file);
	char title[80];
	safe_fread(&title, 80, 1, this->file);
	//printf("Title1: %s\n", title);
	safe_fread(&title, 80, 1, this->file);
	//printf("Title2: %s\n", title);
	//164
	safe_fread(&iin, 4, 1, this->file);
	//4
	safe_fread(&iin, 4, 1, this->file);
	//N
	safe_fread(&N, 4, 1, this->file);
	//4
	safe_fread(&iin, 4, 1, this->file);
}

int DCD::read_frame(float *X, float *Y, float *Z){
	int iin;
	safe_fread(&iin, 4, 1, this->file);
	safe_fread(X, 4*N, 1, this->file);
	safe_fread(&iin, 4, 1, this->file);

	safe_fread(&iin, 4, 1, this->file);
	safe_fread(Y, 4*N, 1, this->file);
	safe_fread(&iin, 4, 1, this->file);

	safe_fread(&iin, 4, 1, this->file);
	safe_fread(Z, 4*N, 1, this->file);
	safe_fread(&iin, 4, 1, this->file);

	if(feof(this->file) == 0){
		return 0;
	} else {
		return -1;
	}
}

void DCD::goto_header(){
    rewind(this->file);
}

void DCD::goto_frame(int frame){
    goto_header();
    read_header();
    const int framesize = 3 * (4 + 4*N + 4);
    fseek(file, framesize*frame, SEEK_CUR);
}

void pad(char *s, int len){
	int curlen;
	int i;
	curlen = strlen(s);
	for (i = curlen; i < len; i++){
		s[i] = ' ';
	}
	s[len] = '\0';
}

