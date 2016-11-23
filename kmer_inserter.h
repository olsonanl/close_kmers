#ifndef _KMER_INSERTER_H
#define _KMER_INSERTER_H

/*
 * Multithreaded kmer/family inserter engine.
 *
 * In order to accelerate the population of the large kmer=>family database,
 * we need to reduce synchronization overhead. We do this by forcing a synchronization
 * where during the insert phase, we guarantee only a single thread will ever
 * insert into the family lists for kmers whose encoded id  modulo N is equal
 * to the worker number, where we have N insertion workers.
 *
 * We feed the insertion workers with a TBB concurrent queue. 
 */

#include "tbb/concurrent_queue.h"
#include <boost/thread/thread.hpp>
#include <vector>
#include <utility>
#include <memory>
#include "kmer.h"

class KmerInserter
{
public:
    struct WorkElement
    {
	std::vector<std::pair<KmerPegMapping::encoded_kmer_t, KmerPegMapping::encoded_family_id_t>> work;
	bool shutdown;
    };
    typedef tbb::concurrent_bounded_queue<std::shared_ptr<WorkElement>> WorkQueue;

    KmerInserter(int n_workers, std::shared_ptr<KmerPegMapping> mapping);

    int n_workers() { return n_workers_; }

    void start();
    void stop();
    void push_work(int modulus, std::shared_ptr<WorkElement> work) {
	queues_[modulus].push(work);
    }

private:
    
    boost::mutex mtx_;
    int n_workers_;
    std::shared_ptr<KmerPegMapping> mapping_;
    std::vector<WorkQueue> queues_;

    boost::thread_group thread_pool_;
    void thread_main(int modulus);
};

#endif
