#include <iostream>
#include <map>
#include <algorithm>

static bool elt_compare(const std::pair<int, int> &a, const std::pair<int, int> &b)
{
    return a.second < b.second;
}

main()
{
    std::map<int, int> x;

    x[0] = 5;

    auto elt = std::max_element(x.begin(), x.end(), [](auto a, auto b) { return a.second < b.second; });

    std::cout << elt->first << " " << elt->second << "\n";
}
