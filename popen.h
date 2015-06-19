#ifndef _POPEN_H
#define _POPEN_H

using namespace std;
#include <cstdio>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>

#include <stdio.h>

#include <boost/noncopyable.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

namespace io = boost::iostreams;

class Popen: private boost::noncopyable {
public:
    explicit Popen(const char* command):
    m_stream(popen(command, "r")) {
	if (!m_stream) throw runtime_error("popen failed");
    }

    ~Popen() {
	pclose(m_stream);
    }

    FILE* stream() const {
	return m_stream;
    }

private:
    FILE* m_stream;
};

#endif
