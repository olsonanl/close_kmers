#ifdef USE_NUMA

#include "numa.h"

#include <hwloc.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

Numa::Numa()
{
    hwloc_topology_init(&topology_);
    hwloc_topology_load(topology_);
    
    int n = hwloc_get_nbobjs_by_type(topology_, HWLOC_OBJ_NODE);
    
    compute_node_ = hwloc_get_obj_by_type(topology_, HWLOC_OBJ_NODE, n - 1);

    hwloc_cpuset_t set = compute_node_->cpuset;

    
    unsigned int id;

    /*
    hwloc_bitmap_foreach_begin(id, set)
    {
	hwloc_bitmap_t xb = hwloc_bitmap_alloc();
	hwloc_bitmap_only(xb, id);
	chosen_cpus_.push_back(xb);
    }
    hwloc_bitmap_foreach_end();
    */
    for (int id = 0; id < 12; id++)
    {
	hwloc_bitmap_t xb = hwloc_bitmap_alloc();
	hwloc_bitmap_only(xb, id);
	chosen_cpus_.push_back(xb);
    }
}

void Numa::bind_memory()
{
    int rc = hwloc_set_membind(topology_, compute_node_->cpuset, HWLOC_MEMBIND_DEFAULT, 0);
    std::cerr  << "hwloc_set_membind returns " << rc << "\n";
    if (rc < 0)
    {
	std::cerr << "hwloc_set_membind failed: " << strerror(errno) << "\n";
	exit(1);
    }
}

void Numa::bind_index(int i)
{
    hwloc_cpuset_t cpuset = cpuset_for_index(i);
    int rc = hwloc_set_cpubind(topology_, cpuset, HWLOC_CPUBIND_THREAD);
    std::cerr << "Bind " << i << " to " << cpuset << " returns " << rc << "\n";
    
}

hwloc_cpuset_t &Numa::cpuset_for_index(int i)
{
    return chosen_cpus_[i % chosen_cpus_.size()];
}

#if 0

int main()
{
    Numa numa;
    for (auto x: numa.chosen_cpus_)
    {
	std::cout << x << "\n";
    }
}

#endif

#endif
