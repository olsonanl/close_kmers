#include <algorithm>
#include <functional>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <map>
#include <boost/regex.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>

namespace acc = boost::accumulators;

int main(int argc, char **argv)
{
    acc::accumulator_set<float, acc::stats<acc::tag::mean,
						  acc::tag::median,
						  acc::tag::variance> > acc;

    acc(16450);
    acc(250);
    
    auto mean = acc::mean(acc);
    auto median = acc::median(acc);
    auto var = acc::variance(acc);

    std::cout << mean << " " << median << " " << var << " " << std::sqrt(var) << "\n";
}
