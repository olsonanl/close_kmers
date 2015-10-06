#ifndef _DEBUG_H
#define _DEBUG_H

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "global.h"

inline bool debug_is_set(const std::string &key)
{
    return g_parameters->count(key) > 0 && (*g_parameters)[key].as<bool>();
}

#endif
