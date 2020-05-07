#include "tbb/parallel_for.h"
#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/compat/thread"
#include "tbb/concurrent_hash_map.h"

#include <iostream>
#include <map>
#include <algorithm>

static bool elt_compare(const std::pair<int, int> &a, const std::pair<int, int> &b)
{
    return a.second < b.second;
}

int main()
{
    tbb::concurrent_unordered_multimap<int, int> x;
    x.insert({1,2});
    x.insert({3,4});
    x.insert({3,5});
    x.insert({3,48});
    x.insert({5,6});

    tbb::parallel_for(x.range(), [](auto r) {
	    std::cout << "Start\n";
	    for (auto i = r.begin(); i != r.end(); i++)
	    {
		std::cout << i->first << " " << i->second << "\n";
	    }
	});
    return 0;
}

void xmain()
{
    std::map<int, int> x;

    x[0] = 5;

    auto elt = std::max_element(x.begin(), x.end(), [](auto a, auto b) { return a.second < b.second; });

    std::cout << elt->first << " " << elt->second << "\n";
}
