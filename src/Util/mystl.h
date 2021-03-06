#pragma once

//
// mystl.h
// My Little Functions: Templates Are Magic
//
// This file contains small functions which just make life a bit easier
//

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

// Get std::map element. If element does not exist, return some default value
template <typename K, typename V>
V get_map_wd(const std::map <K,V> &m, const K &key, const V &defval ) {
    typename std::map<K,V>::const_iterator it = m.find( key );
    if ( it == m.end() ) {
        return defval;
    } else {
        return it->second;
    }
}

// Convert anything to string
template <typename T>
inline std::string any2str(const T& in) {
    std::stringstream s;
    s << in;
    return s.str();
}

template <>
inline std::string any2str(const std::vector<int>& in) {
    if(in.empty()){
        return std::string("empty");
    }
    std::stringstream s;
    for(std::vector<int>::const_iterator it = in.begin(), ite = in.end(); it != ite; ++it){
        s << *it << " ";
    }
    return s.str();
}


template <typename T>
inline T str2any(const char* in) {
    std::stringstream s(in);
    T out = T();
    s >> out;
    return out;
}

template <typename T>
inline T str2any(const std::string& in) {
    std::stringstream s(in);
    T out = T();
    s >> out;
    if (s.fail()) {
        throw std::invalid_argument(std::string("Error converting string `") + in + "` to value.");
    }
    return out;
}

template <>
inline std::string str2any<std::string>(const char *in) {
    std::string ret(in);
    ret.erase( 0                                , ret.find_first_not_of( " \t" ) );
    ret.erase( ret.find_last_not_of( " \t" ) + 1, ret.max_size()                 );
    return ret;
}

// Function to replace inclusions of oldStr in str with newStr
// You know, like 'string::replace' in normal languages
std::string string_replace(const std::string& str, const std::string& oldStr, const std::string& newStr);
template <typename T>
std::string string_replace(const std::string& str, const std::string& oldStr, const T& newStr) {
    return string_replace(str, oldStr, any2str<T>(newStr));
}

// ToLower function
std::string to_lower(const std::string& str);
// Trim leading and ending spaces
std::string trim(const std::string &str);

// Case-insensitive string comparison
bool striequals(const std::string& a, const std::string& b);

// Function to calculate gretest common divisor
// Not like it belongs here or anything...
// I haven't found a better place for it
// I'm really sorry
template <typename T>
T GCD(T a, T b) {
    T t;
    while (b != 0) {
        t = b;
        b = a % b;
        a = t;
    }
    return a;
}

