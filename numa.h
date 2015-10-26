#ifndef _numa_h
#define _numa_h

#include <hwloc.h>
#include <iostream>
#include <vector>

class Numa
{
public:
    Numa();

    hwloc_topology_t topology_;
    std::vector<hwloc_cpuset_t> chosen_cpus_;

    hwloc_cpuset_t &cpuset_for_index(int i);
    hwloc_obj_t compute_node_;

    void bind_index(int i);

    void bind_memory();
};

inline std::ostream &operator<<(std::ostream &os, const hwloc_cpuset_t &set)
{
    char x[100];
    hwloc_bitmap_taskset_snprintf(x, sizeof(x), set);
    os << x;
    return os;
}

#endif
