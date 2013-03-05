/*
 * sop.h
 *
 *  Created on: Feb 27, 2009
 *      Author: zhmurov
 */

extern int checkCovalent(int i, int j);
extern int checkNative(int i, int j);
extern int checkPairs(int i, int j);
extern int checkPossiblePairs(int i, int j);

extern float getDistanceBeads(int i, int j);
extern float getOccupancy(int i);
extern float getBeta(int i);

extern float getX(int i);
extern float getY(int i);
extern float getZ(int i);

extern float getR0(int i, int j);
extern float getEh(int i, int j);
extern int getBeadMask(int i);
