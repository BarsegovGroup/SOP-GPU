/*
 * config_reader.h
 *
 *  Created on: Feb 27, 2009
 *      Author: zhmurov
 */

extern void parseFile(char* filename);

extern int getParameter(char* paramValue, const char* paramName, const char* defaulValue, int allowDefault);
extern int getIntegerParameter(const char* paramName, int defaultValue, int allowDefault);
extern long long int getLongIntegerParameter(const char* paramName, long defaultValue, int allowDefault);
extern float getFloatParameter(const char* paramName, float defaulValue, int allowDefault);
extern int getYesNoParameter(const char* paramName, int defaultValue, int allowDefault);
extern int getVectorParameter(const char* paramName, float* x, float* y, float* z, float defaultX, float defaultY, float defaultZ, int allowDefault);

extern int getMaskedParameter(char* result, const char* paramName, const char* defaultValue, int allowDefault);
extern void applyMask(char* result, const char* parameterMask);
