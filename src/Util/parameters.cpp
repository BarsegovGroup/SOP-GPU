#include "parameters.h"

// This file provides nice wrapper around IO/configreader
namespace parameters {

bool _is_defined(const char *name) { // TODO: probably change name
    return configreader::hasParameter(name);
}

void _initialize(const char* filename, int argc, char *argv[]) { // TODO: probably change name
    return configreader::parseParametersFile(filename, argc, argv);
}

} // namespace parameters

