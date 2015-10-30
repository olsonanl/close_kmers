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

class KmerRequestServer : public std::enable_shared_from_this<KmerRequestServer>
{
public:
    KmerRequestServer(boost::asio::io_service& io_service,
		      const std::string &port,
		      const std::string &port_file,
		      std::shared_ptr<ThreadPool> thread_pool);

    void load_families_nr(std::shared_ptr<KmerPegMapping> &mapping,
			  const std::string &file);
    void startup();
    void deactivate(std::shared_ptr<KmerRequest2> x);

private:

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

    /*
     * NR loading support.
     */
    typedef std::pair<std::string, std::string> seq_t;
    typedef std::vector<seq_t> seq_list_t;

    std::shared_ptr<KmerPegMapping> ld_root_mapping_;
    std::shared_ptr<seq_list_t> ld_cur_work_;
    int ld_max_size_;
    int ld_cur_size_;
    boost::condition_variable ld_cond_;
    boost::mutex ld_mut_;
    int ld_pending_count_;
    int on_parsed_seq(const std::string &id, const std::string &seq);
    void thread_load(std::shared_ptr<seq_list_t> sent_work);
    void on_hit(sig_kmer_t &hit, KmerPegMapping::encoded_id_t &enc_id);

};

#endif
