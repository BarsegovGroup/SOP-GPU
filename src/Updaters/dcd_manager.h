/*
 * dcd_manager.cuh
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#pragma once

#include "../gsop.h"

#define DCD_FREQUENCY_STRING		"dcdfreq"
#define DCD_FILENAME_STRING			"DCDfile"

#define DEFAULT_DCD_FREQUENCY		10000
#define DEFAULT_DCD_FILENAME		"<name>_<run>_<stage>.dcd"

class DcdOutputManager : public SOPUpdater{
public:
    DcdOutputManager();
    virtual ~DcdOutputManager() { }
    virtual void update();
};

void createDCDOutputManager();

