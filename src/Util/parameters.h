#pragma once

#include "../IO/configreader.h"
#include "../Util/mystl.h"

// This file provides nice wrapper around IO/configreader

#define VERBOSE_CONFIG 1

// Parameter with default value
#define PARAMETER(name, ctype, defval, units, description) \
    namespace parameters { \
        static ParameterOptional<ctype> name(#name, defval); \
    }

// Parameter with no default value
#define PARAMETER_MANDATORY(name, ctype, units, description) \
    namespace parameters { \
        static ParameterMandatory<ctype> name(#name); \
    }

// Lazy eveluated parameter: the expression for default value is evaluated when calling ".get()"
#define PARAMETER_LAZY(name, ctype, defexpr, units, description) \
    namespace parameters { \
        class _ParameterLazy_##name : public Parameter<ctype> { \
            public: \
            virtual ctype get() const { \
                return ParameterOptional<ctype>(#name, ((defexpr))).get(); \
            } \
        }; \
        static _ParameterLazy_##name name; \
    }

namespace parameters {

template <typename T>
inline T log_parameter(const std::string& name, const T& val){
#ifdef VERBOSE_CONFIG
    printf("Using parameter '%s' = '%s'\n", name.c_str(), any2str<T>(val).c_str());
#endif
    return val;
}

template <typename T>
class Parameter {
public:
    virtual T get() const = 0;
};

template <typename T>
class ParameterMandatory : public Parameter<T> {
public:
    ParameterMandatory(const char *name) : _name(name) { }
    // Get parameter value
    virtual T get() const { return log_parameter(this->_name, configreader::getParameterAs<T>(this->_name.c_str())); }
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
    ParameterOptional(const char *name, T defval) : ParameterMandatory<T>(name), _defval(defval), _defpar(NULL) { }
    ParameterOptional(const char *name, const Parameter<T>& defpar) : ParameterMandatory<T>(name), _defpar(&defpar) { }
    // Get parameter value
    virtual T get() const { return log_parameter(this->_name, configreader::getParameterAs<T>(this->_name.c_str(), (this->_defpar ? this->_defpar->get() : this->_defval))); }
protected:
    T _defval;
    const Parameter<T>* _defpar;
};

bool _is_defined(const char *name); // TODO: probably change name

void _initialize(const char* filename, int argc = 0, char *argv[] = NULL); // TODO: probably change name

} // namespace parameters

