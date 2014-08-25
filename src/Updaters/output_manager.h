/*
 * output_manager.cuh
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#pragma once

#define OUTPUT_FREQUENCY_STRING			"outputtiming"
#define OUTPUT_FILENAME_STRING			"outputname"

#define DEFAULT_OUTPUT_FREQUENCY		10000
#define DEFAULT_OUTPUT_FILENAME			"energy.<name>_<run>_<stage>.dat"

struct OutputData{
	long int step;
	float tempav;
	float temp;
	//float endToEnd;
	//float endToEnd_x;
	//float f;
	//float fx;
	//float fy;
	//float fz;
	float epot_tot;
	float epot_native;
	float epot_longrange;
	float epot_LJ;
	float epot_fene;
	int nat_num;
	//float xt;
	float rg;
	//float V;
	//float R;
	float tea_eps;
};

void createOutputManager();
void initOutputManager();
void closeDAT();

