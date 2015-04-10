/*
 * configreader.h
 *
 *  Created on: Feb 27, 2009
 *      Author: zhmurov
 */

#pragma once

#include <cstring>
#include <sstream>
#include <vector>
#include <vector_types.h> // For float4
#include "../Util/mystl.h"

namespace configreader {

void parseParametersFile(const char* filename, int argc = 0, char *argv[] = NULL);

bool hasParameter(const char* paramName);
void setParameter(const char* paramName, const char* paramValue);

const std::string getParameter(const char* paramName);
const std::string getMaskedParameter(const char* paramName);

bool getYesNoParameter(const char* paramName);
float3 getFloat3Parameter(const char* paramName);

template <typename T>
std::vector<T> getArrayParameter(const char* paramName);

//    0
//  _/\\_
//  `/-=|
//  O' `O

template <typename T>
inline int setParameterAs(const char* paramName, T paramValue) {
    return setParameter(paramName, any2str(paramValue).c_str());
}

template <typename T>
inline T getRawParameterAs(const char* paramName) {
    // Get parameter without applying masks
    const std::string paramValue = getParameter(paramName);
    return str2any<T>(paramValue);
}

template <typename T>
inline T getParameterAs(const char* paramName) {
    const std::string paramValue = getMaskedParameter(paramName);
    return str2any<T>(paramValue);
}

template <>
inline bool getParameterAs<bool>(const char* paramName) {
    return getYesNoParameter(paramName);
}

template <>
inline std::vector<int> getParameterAs< std::vector<int> >(const char *paramName) {
    return getArrayParameter<int>(paramName);
}

template <>
inline float3 getParameterAs<float3>(const char *paramName) {
    return getFloat3Parameter(paramName);
}

template <typename T>
inline T getParameterAs(const char* paramName, T defVal) {
    if (hasParameter(paramName)){
        return getParameterAs<T>(paramName);
    }else{
        return defVal;
    }
}

} // namespace configreader

