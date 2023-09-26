/*
 * Delimited streams, per https://stackoverflow.com/questions/30073287/c-cout-auto-separator
 */

#ifndef _tabsep_h
#define _tabsep_h

#include <iostream>
#include <algorithm>

struct with_separator {
    with_separator(std::string sep)
	: sep(std::move(sep)) {}

    std::string sep;
};

struct separated_stream {
    separated_stream(std::ostream &stream, std::string sep)
	: _stream(stream), _sep(std::move(sep)), _first(true) {}

        template <class Rhs>
	separated_stream &operator << (Rhs &&rhs) {
	    if(_first)
		_first = false;
	    else
		_stream << _sep;

	    _stream << std::forward<Rhs>(rhs);
	    return *this;
	}

    separated_stream &operator << (std::ostream &(*manip)(std::ostream&)) {
	manip(_stream);
	return *this;
    }

private:
    std::ostream &_stream;
    std::string _sep;
    bool _first;
};

separated_stream operator << (std::ostream &stream, with_separator wsep) {
    return separated_stream(stream, std::move(wsep.sep));
}

#endif
