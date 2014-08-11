/*
 * dcd_manager.cuh
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#include "../IO/dcdio.h"
#include "../def_param.h"
#include "../gsop.cuh"

#ifndef DCD_MANAGER_CUH_
#define DCD_MANAGER_CUH_

#define DCD_FREQUENCY_STRING		"dcdfreq"
#define DCD_FILENAME_STRING			"DCDfile"

#define DEFAULT_DCD_FREQUENCY		10000
#define DEFAULT_DCD_FILENAME		"<name>_<run>_<stage>.dcd"

void createDCDOutputManager();
void initDCD();
void saveCoord();
void closeDCD();

extern void savePDB(char* filename);

SOPUpdater dcdOutputManager;

#endif /* DCD_MANAGER_CUH_ */