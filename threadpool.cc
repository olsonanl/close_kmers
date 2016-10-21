#include "threadpool.h"

#include <iostream>
#include <memory>
#include "kguts.h"

ThreadPool::ThreadPool(const std::string &kmer_dir) :
    next_thread_id_(0),
    kmer_dir_(kmer_dir)
{
}

ThreadPool::~ThreadPool()
{
    stop();
}

void ThreadPool::start(int n_threads)
{
    image_ = std::make_shared<KmerImage>(kmer_dir_);

    // std::cout << "image=" << image_ << "\n";
    work_ = std::make_unique<boost::asio::io_service::work>(io_service_);
    for (int thread = 0; thread < n_threads; thread++)
    {
	int thread_id = next_thread_id_++;
	thread_pool_.create_thread([this, thread_id](){
		// std::cout << "setting up thread " << thread_id << "\n";
#ifdef USE_NUMA
		numa_.bind_index(thread_id);
#endif
		
		kguts_.reset(new KmerGuts(kmer_dir_, image_));
		thread_index_.reset(new int(thread_id));
		io_service_.run();
	    });
    }	
}

void ThreadPool::add_threads(int n_threads)
{
    for (int thread = 0; thread < n_threads; thread++)
    {
	int thread_id = next_thread_id_++;
	thread_pool_.create_thread([this, thread_id](){
		// std::cout << "setting up thread " << thread_id << "\n";
#ifdef USE_NUMA
		numa_.bind_index(thread_id);
#endif
		
		kguts_.reset(new KmerGuts(kmer_dir_, image_));
		thread_index_.reset(new int(thread_id));
		io_service_.run();
	    });
    }	
}

void ThreadPool::stop()
{
    work_.reset();
    io_service_.stop();
    thread_pool_.join_all();
}
