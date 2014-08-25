/*
 * dcd_manager.cuh
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#pragma once

#define DCD_FREQUENCY_STRING		"dcdfreq"
#define DCD_FILENAME_STRING			"DCDfile"

#define DEFAULT_DCD_FREQUENCY		10000
#define DEFAULT_DCD_FILENAME		"<name>_<run>_<stage>.dcd"

void createDCDOutputManager();
void initDCD();
void saveCoord();

