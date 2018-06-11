#ifndef _seed_utils_h
#define _seed_utils_h

#include <string>
#include <vector>
#include <boost/regex.hpp>

namespace seed_utils {

    static const boost::regex strip_func_regex("(\\s*\\#.*$)|"
					  "(^FIG\\d{6}[^:]*:\\s*)");
    static const boost::regex strip_func_comment_regex("(\\s*\\#.*$)");


    static const boost::regex split_function_regex("\\s+[/@]\\s+|\\s*;\\s+");
    
    inline std::string strip_func(const std::string &str)
    {
	return boost::regex_replace(str, strip_func_regex, "");
    }
    
    inline std::string strip_func_comment(const std::string &str)
    {
	return boost::regex_replace(str, strip_func_comment_regex, "");
    }
    
    inline std::vector<std::string> roles_of_function(const std::string function)
    {
	std::string stripped(strip_func_comment(function));
	boost::sregex_token_iterator x(stripped.begin(), stripped.end(), split_function_regex, -1);
	boost::sregex_token_iterator e;

	std::vector<std::string> ret;
	for (; x != e; x++)
	{
	    ret.emplace_back(*x);
	}
	return ret;
    }
};

#endif // _seed_utils_h
