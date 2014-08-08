//
// mystl.cpp
// My Little Functions: Templates Are Magic
//
// This file contains small functions which make life a bit easier
//

#include "mystl.h"
#include <algorithm>

// Function to replace inclusions of oldStr in str with newStr
// You know, like 'string::replace' in normal languages
std::string string_replace(const std::string& str, const std::string& oldStr, const std::string& newStr) {
    size_t pos = 0;
    std::string ret(str);
    while((pos = ret.find(oldStr, pos)) != std::string::npos)
    {
        ret.replace(pos, oldStr.length(), newStr);
        pos += newStr.length();
    }
    return ret;
}


std::string lower(const std::string& str) {
    std::string ret(str);
    std::transform(ret.begin(), ret.end(), ret.begin(), ::tolower);
    return ret;
}

