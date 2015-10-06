#ifndef _stringutil_h
#define _stringutil_h

#include <string>

inline bool
endswith(const std::string & s, char sfx)
{
    return !s.empty() && s[s.size() - 1] == sfx;
}
   
inline bool
endswith(const std::string & s, const char * sfx, size_t len)
{
    return s.size() >= len && (std::memcmp(s.data() + s.size() - len, sfx, len) == 0);
}

inline bool
endswith(const std::string & s, const char * sfx)
{
    return endswith(s, sfx, std::strlen(sfx));
}

inline bool
endswith(const std::string & s, const std::string & sfx)
{
    return endswith(s, sfx.data(), sfx.size());
}

#endif
