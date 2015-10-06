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
#include "krequest.h"
#include "krequest2.h"
#include "threadpool.h"

class KmerRequestServer : public std::enable_shared_from_this<KmerRequestServer>
{
public:
    KmerRequestServer(boost::asio::io_service& io_service,
		      const std::string &port,
		      const std::string &port_file,
		      KmerPegMapping &mapping,
		      boost::asio::ip::tcp::endpoint &klookup_endpoint,
		      std::shared_ptr<KmerGuts> kguts,
		      std::shared_ptr<ThreadPool> thread_pool);

    void startup();
    void deactivate(std::shared_ptr<KmerRequest2> x);

private:

    void do_accept();
    void on_accept(boost::system::error_code ec, KmerRequest *);
    void do_accept2();
    void on_accept2(boost::system::error_code ec, std::shared_ptr<KmerRequest2>);
    
    void do_await_stop();
    std::shared_ptr<KmerGuts> kguts_;
    std::shared_ptr<ThreadPool> thread_pool_;
    boost::asio::io_service &io_service_;
    boost::asio::ip::tcp::acceptor acceptor_;
    boost::asio::signal_set signals_;
    std::string port_;
    std::string port_file_;
    std::set<std::shared_ptr<KmerRequest2> > active_;
    KmerPegMapping &mapping_;
    boost::asio::ip::tcp::endpoint klookup_endpoint_;

    std::shared_ptr<std::map<std::string, std::shared_ptr<KmerPegMapping>>> mapping_map_;
};

#endif
