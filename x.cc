#include <algorithm>
#include <functional>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <map>
#include <boost/regex.hpp>
#include "kguts.h"

void doit(KmerGuts::sig_kmer_t &x)
{
    std::cout << "got one\n";
}

void callme(std::function<void(KmerGuts::sig_kmer_t &)> hit_cb)
{
    if (hit_cb)
    {
	KmerGuts::sig_kmer_t x;
	hit_cb(x);
    }
}

template<class T>
struct less_second : std::binary_function<T,T,bool>
{
    inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
};

int main(int argc, char **argv)
{
    int n = 3;
    std::vector<double> vec(5, 1.0 / n);

    for (auto it = vec.begin(); it != vec.end(); it++)
    {
	std::cout << *it <<"\n";
	*it /= 2;
    }
    for (auto it: vec)
    {
	std::cout << it <<"\n";
	it /= 2;
    }
    
    
    const boost::regex url_regex("^(([^:/?#]+):)?(//([^/?#]*))?([^?#]*)(\\?([^#]*))?(#(.*))?");
    const boost::regex path_regex("^([^?#]*)(\\?([^#]*))?(#(.*))?");
    const boost::regex request_regex("^([A-Z]+) ([^?#]*)(\\?([^#]*))?(#(.*))? HTTP/(\\d+\\.\\d+)");

    try {
	int v = std::stoi("");
    
	std::cout << v;
    } catch (const std::invalid_argument& ia)
    {
	std::cout << " erorr: " << ia.what() << "\n";
    }

    #if 0

    std::string l;
    while (std::getline(std::cin, l))
    {
	boost::smatch what;
	if (boost::regex_match(l, what, request_regex))
	{
	    for (auto it = what.begin(); it != what.end(); it++)
	    {
		std::cout << "'" << *it << "'\n";
	    }
	}
	else
	{
	    std::cout << "No. '" << l << "'\n";
	}
	    
    }
    #endif

    std::map<int, int> mm;

    mm[3]++;
    mm[4]++;
    mm[4]++;
    mm[4]++;
    mm[3]++;
    for (auto x = mm.begin(); x != mm.end(); x++)
    {
	std::cout << x->first << " " << x->second << std::endl;
    }

    std::vector<std::pair<int, int> > mm2;
    mm2.insert(mm2.begin(), mm.begin(), mm.end());
    std::sort(mm2.begin(), mm2.end(), less_second<std::pair<int, int>>());
    for (auto x = mm2.begin(); x != mm2.end(); x++)
    {
	std::cout << "copy: " << x->first << " " << x->second << std::endl;
    }
    
    
    std::string x = "file.gz";
    std::string y = x.substr(x.length() - 3);
    std::cout << y << "\n";

    auto list = std::make_shared<std::vector<int>>();
    list->push_back(3);
    list->push_back(4);
    std::cout << (list->back()) << std::endl;

    auto list2 = list;

    list = std::make_shared<std::vector<int>>();
    list->push_back(5);
    std::cout << (list->back()) << std::endl;
    std::cout << (list2->back()) << std::endl;

    auto it = list2->begin();
    while (it != list2->end())
    {
	std::cout << *it << std::endl;
	it++;
    }

    auto khits = std::make_shared<std::vector<KmerHit>>();

    unsigned int m = 3;
    khits->push_back({ m, 1, 10, 11, 1.4, 34 });
    for (auto it2 : *khits)
    {
//	std::cout << it2.offset << " " << it2.encoded_kmer <<std::endl;
    }

    callme(doit);
}
