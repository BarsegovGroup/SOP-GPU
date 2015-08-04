/*
 * config_reader.c
 *
 *  Created on: Jan 23, 2009
 *      Author: zhmurov
 */
#define IO_CONFIG_DEBUG

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "../Util/wrapper.h"
#include "configreader.h"
#include <algorithm>
#include <map>

#define BUF_SIZE 4096

namespace configreader {

std::map<std::string,std::string> parameters;

void parseParametersFile(const char* filename, int argc, char *argv[]){
	FILE* file = safe_fopen(filename, "r");
    char buffer[BUF_SIZE];
	while(fgets(buffer, BUF_SIZE, file) != NULL){
		if(buffer[0] != '#' && buffer[0] != ' ' && buffer[0] != '\n' && buffer[0] != '\t'){
			char* pch = strtok(buffer, " \t");
            const std::string paramName(pch);
			pch = strtok(NULL, "#\n");
            std::string paramValue(pch);
            paramValue = trim(paramValue);

#ifdef IO_CONFIG_DEBUG
			printf("%s\t%s\n", paramName.c_str(), paramValue.c_str());
#endif
			parameters[paramName] = paramValue;
		}
	}
	fclose(file);

    // And now patch parameters according to argc/argv (if we have ones)
    for (int i = 2; i < argc; ++i) { // argv[0] = binary name, argv[1] = config name
        for (int c = 0; argv[i][c] != '\0'; ++c) {
            if (argv[i][c] == '=') {
                argv[i][c] = '\0';
                const std::string paramName(&argv[i][0]);
                const std::string paramValue(&argv[i][c+1]);
                argv[i][c] = '=';
#ifdef IO_CONFIG_DEBUG
                printf("Enforced via argv: %s\t=%s\n", paramName.c_str(), paramValue.c_str());
#endif
                parameters[paramName] = paramValue;
                break;
            }
        }
    }
}

bool hasParameter(const char* paramName) {
    return (parameters.find(paramName) != parameters.end());
}

void setParameter(const char* paramName, const char* paramValue){
    parameters[paramName] = paramValue;
}

const std::string getParameter(const char* paramName){
    std::map<std::string,std::string>::const_iterator param = parameters.find(paramName);
    if (param == parameters.end()){ // paramName not found
		DIE("Parameter '%s' not specified", paramName);
    }else{
        return param->second;
    }
}

const std::string getMaskedParameter(const char* paramName){
    std::string paramValue = getParameter(paramName);
    for (std::map<std::string,std::string>::const_iterator param = parameters.begin(),
                                                     param_end = parameters.end();
                                                     param != param_end; ++param){
        const std::string mask = std::string("<") + param->first + ">";
        paramValue = string_replace(paramValue, mask, param->second);
    }
    return paramValue;
}

template <typename T>
std::vector<T> getArrayParameter(const char* paramName){
    std::string paramValue = to_lower(getParameter(paramName)); // Immediately convert to lowercase to simplify comparisons to magic words later
    std::vector<T> ret;

    // Magical values for empty array
    if ((paramValue == "empty") || (paramValue == "none")){
        return ret;
    }
    std::istringstream iss(paramValue);

    // Split string
    std::string ts;
    while (iss >> ts){
        if (ts == "to"){ // We're working with range
            if (ret.empty()){
                DIE("In parameters of array type, 'to' keyword must be preceeded by number");
            }
            T from = *(ret.end()-1);
            T to;
            iss >> to;
            for (T i = from + 1; i <= to; ++i){
                ret.push_back(i);
            }
        }else{
            T t = str2any<T>(ts);
            ret.push_back(t);
        }
    }
	return ret;
}
template std::vector<int> getArrayParameter<int>(const char* paramName);
template std::vector<float> getArrayParameter<float>(const char* paramName);

bool getYesNoParameter(const char* paramName){
    std::string paramValue = to_lower(getParameter(paramName));
    if (paramValue == "yes" || paramValue == "true"  || paramValue == "on"  || paramValue == "y")
        return true;
    if (paramValue == "no"  || paramValue == "false" || paramValue == "off" || paramValue == "n")
        return false;
    DIE("Wrong value '%s' for boolean parameter '%s'", paramValue.c_str(), paramName);
}

float3 getFloat3Parameter(const char* paramName) {
    const std::vector<float> paramValue = getArrayParameter<float>(paramName);
    if (paramValue.size() != 3){
        DIE("Float3 parameter '%s' has %d components; must have 3.", paramName, paramValue.size());
    }
    float3 t;
    t.x = paramValue[0];
    t.y = paramValue[1];
    t.z = paramValue[2];
    return t;
}

} // namespace configreader

