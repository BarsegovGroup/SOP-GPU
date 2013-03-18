#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
/*
 *  dcdio.c
 *
 * Code to work with dcd files. On Oct 30, 2008 only writing is implemented.
 *
 *  Created on: Oct 29, 2008
 *      Author: zhmurov
 */

FILE* dcd_open_write(FILE* dcd_file, char *dcd_filename);
FILE* dcd_open_append(FILE* dcd_file, char *dcd_filename);
FILE* dcd_open_read(FILE* dcd_file, char *dcd_filename);

void dcd_write_header(FILE* dcd_file, char *dcd_filename, int N, int NFILE, int NPRIV, int NSAVC, double DELTA);
void dcd_write_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z);

void pad(char *s, int len);
void dcd_read_header(FILE* dcd_file, char *dcd_filename, int* N, int* NFILE, int* NPRIV, int* NSAVC, float* DELTA);
int dcd_read_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z);

void dcd_close(FILE* dcd_file);

/*
 * Open DCD file to write data into.
 * The FILE* pointer is returned and must be used for further calls
 */
FILE* dcd_open_write(FILE* dcd_file, char *dcd_filename){
	dcd_file = fopen(dcd_filename, "w");
	return dcd_file;
}

FILE* dcd_open_append(FILE* dcd_file, char *dcd_filename){
	dcd_file = fopen(dcd_filename, "a");
	return dcd_file;
}

/*
 * Write header to dcd file.
 * Inputs are:
 *  dcd_file - FILE* object from dcd_open_write call
 *  dcd_filename - dcd filename to write into remarks
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
void dcd_write_header(FILE* dcd_file, char *dcd_filename, int N, int NFILE, int NPRIV, int NSAVC, double DELTA){
	int iout;
	float fout;
	char cout[5];
	iout = 84;
	fwrite(&iout, 4, 1, dcd_file);
	sprintf(cout, "CORD");
	fwrite(&cout, 4, 1, dcd_file);
	iout = NFILE;
	fwrite(&iout, 4, 1, dcd_file);
	iout = NPRIV;
	fwrite(&iout, 4, 1, dcd_file);
	iout = NSAVC;
	fwrite(&iout, 4, 1, dcd_file);
	iout = NPRIV-NSAVC;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 0;
	for(int i = 0; i < 5; i++){
		fwrite(&iout, 4, 1, dcd_file);
	}
	fout = DELTA;
	fwrite(&fout, 4, 1, dcd_file);
	iout = 0;
	for(int i = 0; i < 9; i++){
		fwrite(&iout, 4, 1, dcd_file);
	}
	iout = 24;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 84;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 164;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 2;
	fwrite(&iout, 4, 1, dcd_file);
	char title[81];
	sprintf(title, "REMARKS FILENAME = %s CREATED BY dcdio.c", dcd_filename);
	pad(title, 80);
	fwrite(&title, 80, 1, dcd_file);
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	sprintf(title, "REMARKS DATE: %s", asctime(timeinfo));
	pad(title, 80);
	fwrite(&title, 80, 1, dcd_file);
	iout = 164;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 4;
	fwrite(&iout, 4, 1, dcd_file);
	iout = N;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 4;
	fwrite(&iout, 4, 1, dcd_file);
}

/*
 * Writes frame into dcd
 * Input
 *  dcd_file - file to write in
 *  N - number of atoms
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
void dcd_write_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z){
	int iout;
	iout = N*4;
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(X, 4*N, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(Y, 4*N, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(Z, 4*N, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
}

/*
 * Close DCD after writing
 */
void dcd_close(FILE* dcd_file){
	fclose(dcd_file);
}


FILE* dcd_open_read(FILE* dcd_file, char *dcd_filename){
	dcd_file = fopen(dcd_filename, "r");
	return dcd_file;
}

void dcd_read_header(FILE* dcd_file, char *dcd_filename, int* N, int* NFILE, int* NPRIV, int* NSAVC, float* DELTA){
	int iin;
	char cin[5];
	fread(&iin, 4, 1, dcd_file);
	if(iin != 84){
		printf("Error! Wrong DCD file: DCD supposed to have '84' in first 4 bytes, but it hasn't.\n");
		exit(-1);
	}
	fread(&cin, 4, 1, dcd_file);
	char cord_char[5];
	sprintf(cord_char, "CORD");
	if(strcmp(cin, cord_char)){
		printf("Error! Wrong DCD file: no 'CORD' sequence at the beginning of the file. Found: %s.\n", cin);
		exit(-1);
	}

	fread(NFILE, 4, 1, dcd_file);
	fread(NPRIV, 4, 1, dcd_file);
	fread(NSAVC, 4, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	for(int i = 0; i < 5; i++){
		fread(&iin, 4, 1, dcd_file);
	}
	fread(DELTA, 4, 1, dcd_file);
	for(int i = 0; i < 9; i++){
		fread(&iin, 4, 1, dcd_file);
	}
	//24
	fread(&iin, 4, 1, dcd_file);
	//84;
	fread(&iin, 4, 1, dcd_file);
	//164;
	fread(&iin, 4, 1, dcd_file);
	//2;
	fread(&iin, 4, 1, dcd_file);
	char title[81];
	fread(&title, 80, 1, dcd_file);
	printf("Title1: %s\n", title);
	fread(&title, 80, 1, dcd_file);
	printf("Title2: %s\n", title);
	//164
	fread(&iin, 4, 1, dcd_file);
	//4
	fread(&iin, 4, 1, dcd_file);
	//N
	fread(N, 4, 1, dcd_file);
	//4
	fread(&iin, 4, 1, dcd_file);
}

int dcd_read_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z){
	int iin;
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(X, 4*N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(Y, 4*N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(Z, 4*N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	if(feof(dcd_file) == 0){
		return 0;
	} else {
		return -1;
	}
}

/*
 * Extend the string to be 80 bytes length (for remarks in header)
 * Taken from dcdlib.C (NAMD source)
 */
void pad(char *s, int len){
	int curlen;
	int i;
	curlen = strlen(s);
	if (curlen > len){
		s[len] = '\0';
		return;
	}
	for (i = curlen; i < len; i++){
		s[i] = ' ';
	}
	s[i] = '\0';
}

