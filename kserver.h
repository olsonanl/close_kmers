#ifndef _KSERVER_H
#define _KSERVER_H


/*
 * Kmer-request server.
 *
 * We speak pidgin HTTP here so we can make this available over a proxy.
 * This is NOT a general purpose HTTP server.
 *
 */

#include <boost/asio.hpp>
#include <set>
#include <memory>
#include "kmer.h"
#include "krequest2.h"
#include "threadpool.h"

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

class KmerRequestServer : public std::enable_shared_from_this<KmerRequestServer>
{
public:
    KmerRequestServer(boost::asio::io_service& io_service,
		      const std::string &port,
		      const std::string &port_file,
		      std::shared_ptr<ThreadPool> thread_pool,
		      bool family_mode = false);

    void load_families_nr(const std::string &file);
    void startup();
    void deactivate(std::shared_ptr<KmerRequest2> x);

    bool family_mode() { return family_mode_; }

    std::shared_ptr<KmerPegMapping> root_mapping() { return root_mapping_; }

private:
    bool family_mode_;

    void do_accept2();
    void on_accept2(boost::system::error_code ec, std::shared_ptr<KmerRequest2>);
    
    void do_await_stop();
    std::shared_ptr<ThreadPool> thread_pool_;
    boost::asio::io_service &io_service_;
    boost::asio::ip::tcp::acceptor acceptor_;
    boost::asio::signal_set signals_;
    std::string port_;
    std::string port_file_;
    std::set<std::shared_ptr<KmerRequest2> > active_;

    std::shared_ptr<std::map<std::string, std::shared_ptr<KmerPegMapping>>> mapping_map_;
    std::shared_ptr<KmerPegMapping> root_mapping_;

    /*
     * NR loading support.
     */
public:
    typedef std::pair<std::string, std::string> seq_t;
    typedef std::vector<seq_t> seq_list_t;

private:
    NRLoadState load_state_;

};

class NRLoader
{
public:
    NRLoader(NRLoadState &load_state, const std::string &file, std::shared_ptr<KmerPegMapping> root_mapping,
	     std::shared_ptr<ThreadPool> thread_pool, size_t n_files, bool family_mode);

    void start();
    void load_families();

    bool family_mode_;

    NRLoadState &load_state_;
    NRLoadState my_load_state_;
    std::string file_;
    std::shared_ptr<KmerPegMapping> root_mapping_;
    std::shared_ptr<ThreadPool> thread_pool_;
    size_t n_files_;
	
    std::shared_ptr<KmerRequestServer::seq_list_t> cur_work_;
    size_t max_size_;
    size_t cur_size_;
    boost::mutex dbg_mut_;
    int chunks_started_;
    int chunks_finished_;
    int on_parsed_seq(const std::string &id, const std::string &seq);
    void thread_load(std::shared_ptr<KmerRequestServer::seq_list_t> sent_work, int count);
    void on_hit(KmerGuts::hit_in_sequence_t hit, KmerPegMapping::encoded_id_t &enc_id, size_t seq_len);
    void on_hit_fam(KmerGuts::hit_in_sequence_t hit, KmerPegMapping::encoded_family_id_t &enc_id, size_t seq_len);

};

#endif
