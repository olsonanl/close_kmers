#ifndef _PARALLEL_READ_H
#define _PARALLEL_READ_H

#include <string>

template <class Callback>
void parallel_read(const std::string &file, int n_threads, Callback cb);

#include "parallel_read.cc"

#endif
