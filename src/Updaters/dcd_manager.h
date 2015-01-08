/*
 * dcd_manager.cuh
 *
 *  Created on: Apr 8, 2010
 *      Author: zhmurov
 */

#pragma once

#include "../gsop.h"
#include "../Util/parameters.h"

PARAMETER(dcdfreq, int, 10000, "steps", "...")
PARAMETER(DCDfile, std::string, "<name>_<author><run>_<stage>.dcd", "path", "...")
PARAMETER(restartfreq, int, 100000, "steps", "...")
PARAMETER(restartname, std::string, "<name>_<author><run>_restart", "path", "...")

class DcdOutputManager : public SOPUpdater{
public:
    DcdOutputManager();
    virtual ~DcdOutputManager() { }
    virtual void update();
};

void createDCDOutputManager();

