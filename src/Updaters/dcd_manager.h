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
PARAMETER(DCDfile, std::string, "<name>_<run>.dcd", "path", "...")
PARAMETER(restartfreq, int, 100000, "steps", "...")
PARAMETER(restartname, std::string, "restart.<name>_<run>.pdb", "path", "...")

class DcdOutputManager : public SOPUpdater{
public:
    DcdOutputManager();
    virtual ~DcdOutputManager() { }
    virtual void update();
};

void createDCDOutputManager();

