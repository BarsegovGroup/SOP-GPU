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

std::string to_lower(const std::string& str){
    std::string ret(str);
    std::transform(ret.begin(), ret.end(), ret.begin(), ::tolower);
    return ret;
}

std::string trim(const std::string& str){
    std::string::const_iterator it = str.begin();
    while (it != str.end() && isspace(*it)) it++;

    std::string::const_reverse_iterator rit = str.rbegin();
    while (rit.base() != it && isspace(*rit)) rit++;

    return std::string(it, rit.base());
}

bool striequals(const std::string& a, const std::string& b)
{
    // Taken from http://stackoverflow.com/a/4119881/929437
    size_t sz = a.size();
    if (b.size() != sz)
        return false;
    for (size_t i = 0; i < sz; ++i)
        if (std::tolower(a[i]) != std::tolower(b[i]))
            return false;
    return true;
}

