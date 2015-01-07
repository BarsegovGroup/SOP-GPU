#pragma once

#include "../IO/configreader.h"
#include "../Util/mystl.h"

#define PARAMETER(name, ctype, defval, units, description) \
    namespace parameters { \
        static ParameterOptional<ctype> name(#name, defval); \
    }
    

#define PARAMETER_MANDATORY(name, ctype, units, description) \
    namespace parameters { \
        static ParameterMandatory<ctype> name(#name); \
    };

template <typename T>
class ParameterMandatory {
public:
    ParameterMandatory(const char *name) : _name(name) { }
    // Get parameter value
    virtual T get() const { return getMaskedParameterAs<T>(this->_name.c_str()); }
    // Get parameter value with replacement
    template <typename T2>
    T replace(const char *str, const T2& value) const { return string_replace(any2str<T>(this->get()), std::string(str), any2str<T2>(value)); }
    // Various introspection stuff
    const std::string& name() const { return _name; }
protected:
    std::string _name;
};

template <typename T>
class ParameterOptional : public ParameterMandatory<T> {
public:
    ParameterOptional(const char *name, T defval) : ParameterMandatory<T>(name), _defval(defval) { }
    // Get parameter value
    virtual T get() const { return getMaskedParameterAs<T>(this->_name.c_str(), this->_defval); }
protected:
    T _defval;
};

