#ifndef _nr_loader_h
#define _nr_loader_h

#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <atomic>

#include "kmer.h"
#include "threadpool.h"
#include "kmer_inserter.h"

class NRLoadState
{
public:
    int pending_count_;
    bool action_seen_;
    boost::condition_variable cond_;
    boost::mutex mut_;
    std::string name_;

    NRLoadState(const std::string &name) : pending_count_(0), action_seen_(false), name_(name) {}
    ~NRLoadState() {
	// std::cerr << name_ << " destroy with pending=" << pending_count_ << "\n";
    }
    
    void pending_inc() {
	boost::lock_guard<boost::mutex> lock(mut_);
	pending_count_++;
	// std::cerr << name_ << " incremented to pending=" << pending_count_ << "\n";
	action_seen_ = true;
    };

    void pending_dec()
    {
	{
	    boost::lock_guard<boost::mutex> lock(mut_);
	    // std::cerr << name_ << " finishing pending=" << pending_count_ << "\n";
	    pending_count_--;
	}
	cond_.notify_one();
    }

    void pending_wait()
    {
	boost::unique_lock<boost::mutex> lock(mut_);
	while (pending_count_ > 0)
	{
	    // std::cerr << name_ << " pending=" << pending_count_ << "\n";
	    cond_.wait(lock);
	}
	// std::cerr << name_ << " pending wait complete count=" << pending_count_ << "\n";
    }
};

class NRLoader
{
public:
    NRLoader(NRLoadState &load_state, const std::string &file, std::shared_ptr<KmerPegMapping> root_mapping,
	     std::shared_ptr<ThreadPool> thread_pool, size_t n_files, bool family_mode,
	     KmerInserter &inserter);

    typedef std::pair<std::string, std::string> seq_t;
    typedef std::vector<seq_t> seq_list_t;

    void start();
    void load_families();

    bool family_mode_;
    KmerInserter &inserter_;

    NRLoadState &load_state_;
    NRLoadState my_load_state_;
    std::string file_;
    std::shared_ptr<KmerPegMapping> root_mapping_;
    std::shared_ptr<ThreadPool> thread_pool_;
    size_t n_files_;
	
    std::shared_ptr<seq_list_t> cur_work_;
    size_t max_size_;
    size_t cur_size_;
    boost::mutex dbg_mut_;
    std::atomic<int> chunks_started_;
    std::atomic<int> chunks_finished_;
    int on_parsed_seq(const std::string &id, const std::string &seq);
    void thread_load(std::shared_ptr<seq_list_t> sent_work, int count);
    void on_hit(const KmerGuts::hit_in_sequence_t &hit, KmerPegMapping::encoded_id_t &enc_id, size_t seq_len);
    void on_hit_fam(const KmerGuts::hit_in_sequence_t &hit, KmerPegMapping::encoded_family_id_t &enc_id, size_t seq_len);

};

#endif // _nr_loader_h

