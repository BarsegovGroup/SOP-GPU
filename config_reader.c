/*
 * config_reader.c
 *
 *  Created on: Jan 23, 2009
 *      Author: zhmurov
 */
#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define buf_size 256
#define name_length 50
#define value_length 100

int paramCount;
char** paramNames;
char** paramValues;

void parseFile(char* filename);

int getParameter(char* paramValue, const char* paramName, const char* defaulValue, int allowDefault);
int getIntegerParameter(const char* paramName, int defaultValue, int allowDefault);
long long int getLongIntegerParameter(const char* paramName, long defaultValue, int allowDefault);
float getFloatParameter(const char* paramName, float defaulValue, int allowDefault);
int getYesNoParameter(const char* paramName, int defaultValue, int allowDefault);
int getVectorParameter(const char* paramName, float* x, float* y, float* z, float defaultX, float defaultY, float defaultZ, int allowDefault);

int getMaskedParameter(char* result, const char* paramName, const char* defaultValue, int allowDefault);
void applyMask(char* result, const char* parameterMask);
void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace);

void parseFile(char* filename){
	FILE* file = fopen(filename, "r");
	if(file != NULL){
		printf("Parsing '%s' parameters file...\n", filename);
	} else {
		printf("ERROR: Parameters file '%s' can not be found. Program will exit.\n", filename);
		exit(-1);
	}
	paramCount = 0;
	char buffer[buf_size];
	while(fgets(buffer, buf_size, file) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			//printf("%s\n", buffer);
			paramCount++;
		}
	}
	paramNames = (char**)malloc(name_length*paramCount*sizeof(char));
	paramValues = (char**)malloc(value_length*paramCount*sizeof(char));
	int i;
	for(i = 0; i < paramCount; i++){
		paramNames[i] = (char*)malloc(name_length*sizeof(char));
		paramValues[i] = (char*)malloc(value_length*sizeof(char));
	}
	rewind(file);
	i = 0;
	while(fgets(buffer, buf_size, file) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			char* pch = strtok(buffer, " \t");
			strcpy(paramNames[i], pch);
			pch = strtok(NULL, "#\n");
			strcpy(paramValues[i], pch);
#ifdef DEBUG
			printf("%s\t%s\n", paramNames[i], paramValues[i]);
#endif
			i++;
		}
	}
	fclose(file);
}

int getParameter(char* paramValue, const char* paramName, const char* defaultValue, int allowDefault){
	int i;
	for(i = 0; i < paramCount; i++){
		if(strcmp(paramName, paramNames[i]) == 0){
			strcpy(paramValue, paramValues[i]);
#ifdef DEBUG
			printf("'%s' = '%s'\n", paramName, paramValue);
#endif
			return 0;
		}
	}
	if(allowDefault){
		strcpy(paramValue, defaultValue);
		printf("Using default value for parameter %s: '%s' = '%s'\n", paramName, paramName, defaultValue);
	} else {
		printf("Parameter '%s' should be specified in a configuration file. Program will exit.\n", paramName);
		exit(-1);
	}
	return 0;
}

int getIntegerParameter(const char* paramName, int defaultValue, int allowDefault){
	char paramValue[value_length];
	char defaultString[value_length];
	sprintf(defaultString, "%d\0", defaultValue);
	int error = getMaskedParameter(paramValue, paramName, defaultString, allowDefault);
	int result = atoi(paramValue);
	if(result == 0 && strcmp(paramValue, "0") != 0){
		printf("ERROR: Wrong value of %s in a configuration file ('%s'). Should be integer. Program will exit.\n", paramName, paramValue);
		exit(-1);
	}
	if(error != 0){
		return 0;
	}
	return result;
}

long long int getLongIntegerParameter(const char* paramName, long defaultValue, int allowDefault){
	char paramValue[value_length];
	char defaultString[value_length];
	sprintf(defaultString, "%ld\0", defaultValue);
	int error = getMaskedParameter(paramValue, paramName, defaultString, allowDefault);
	long long int result = atol(paramValue);
	if(result == 0 && strcmp(paramValue, "0") != 0){
		printf("ERROR: Wrong value of %s in a configuration file ('%s'). Should be long integer. Program will exit.\n", paramName, paramValue);
		exit(-1);
	}
	if(error != 0){
		return 0;
	}
	return result;
}

float getFloatParameter(const char* paramName, float defaultValue, int allowDefault){
	char paramValue[value_length];
	char defaultString[value_length];
	sprintf(defaultString, "%f\0", defaultValue);
	int error = getMaskedParameter(paramValue, paramName, defaultString, allowDefault);
	float result = atof(paramValue);
	if(result == 0.0 && strcmp(paramValue, "0") != 0 && strcmp(paramValue, "0.0") != 0 && strcmp(paramValue, "0.0f") != 0 &&
			strcmp(paramValue, "0.000000") != 0){
		printf("ERROR: Wrong value of %s in a configuration file ('%s'). Should be float. Program will exit.\n", paramName, paramValue);
		exit(-1);
	}
	if(error != 0){
		return 0.0;
	}
	return result;
}

int getYesNoParameter(const char* paramName, int defaultValue, int allowDefault){
	char paramValue[value_length];
	char defaultString[value_length];
	if(defaultValue){
		sprintf(defaultString, "YES\0");
	} else {
		sprintf(defaultString, "NO\0");
	}
	int error = getMaskedParameter(paramValue, paramName, defaultString, allowDefault);
	if(error != 0){
		return 0;
	}
	if(strcmp(paramValue, "YES") == 0 || strcmp(paramValue, "Yes") == 0 || strcmp(paramValue, "yes") == 0
			|| strcmp(paramValue, "Y") == 0 || strcmp(paramValue, "y") == 0
			|| strcmp(paramValue, "ON") == 0 || strcmp(paramValue, "On") == 0 || strcmp(paramValue, "on") == 0
			|| strcmp(paramValue, "TRUE") == 0 || strcmp(paramValue, "True") == 0 || strcmp(paramValue, "true") == 0){
		return 1;
	} else {
		return 0;
	}
}

int getVectorParameter(const char* paramName, float* x, float* y, float* z, float defaultX, float defaultY, float defaultZ, int allowDefault){
		char paramValue[value_length];
		char defaultString[value_length];
		sprintf(defaultString, "%f %f %f\0", defaultX, defaultY, defaultZ);
		int error = getParameter(paramValue, paramName, defaultString, allowDefault);
		if(error != 0){
			return error;
		}
		char* pch = strtok(paramValue, " \t");
		x[0] = atof(pch);
		pch = strtok(NULL, " \t");
		y[0] = atof(pch);
		pch = strtok(NULL, " \t");
		z[0] = atof(pch);
		return 0;
}

int getMaskedParameter(char* result, const char* paramName, const char* defaultValue, int allowDefault){
	char parameterMask[value_length];
	int error = getParameter(parameterMask, paramName, defaultValue, allowDefault);
	char* parameterMaskTrim = strtok(parameterMask, " \t");
	if(error == 0){
		applyMask(result, parameterMaskTrim);
	}
	return error;
}

void applyMask(char* result, const char* parameterMask){
	char tempstring[100];
	strcpy(tempstring, parameterMask);
	int i;
	for(i = 0; i < paramCount; i++){
		char paramName[name_length];
		sprintf(paramName, "<%s>", paramNames[i]);
		char replacementString[value_length];
		strcpy(replacementString, paramValues[i]);
		replaceString(result, tempstring, strtok(replacementString, " \t"), paramName);
		strcpy(tempstring, result);
	}
}

void replaceString(char* resultString, const char* initialString, const char* replacementString, const char* stringToReplace){
	//printf("Looking for %s in %s.\n", stringToReplace, initialString);
	int i;
	int len1, len2;
	resultString[0] = 0;
	//printf("Result string: %s\n", resultString);
	if(strstr(initialString, stringToReplace) != NULL){
		len1 = strlen(initialString) - strlen(strstr(initialString, stringToReplace));
		//printf("len1: %d\n", len1);
		strncpy(resultString, initialString, len1);
		//printf("Result string: %s\n", resultString);
		strncpy(&resultString[len1], replacementString, strlen(replacementString));
		len2 = len1 + strlen(stringToReplace);
		len1 += strlen(replacementString);
		//printf("Result string: %s\n", resultString);

		strcpy(&resultString[len1], &initialString[len2]);
		strcat(resultString, "\0");
		//printf("Found. Result: %s\n", resultString);
	} else {
		strcpy(resultString, initialString);
	}
}
