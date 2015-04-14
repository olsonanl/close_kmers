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
#include "kmer.h"
#include "krequest.h"

class KmerRequestServer
{
public:
    KmerRequestServer(boost::asio::io_service& io_service,
		      const std::string &port,
		      KmerPegMapping &mapping,
		      boost::asio::ip::tcp::endpoint &klookup_endpoint);

private:

    void do_accept();
    
    void do_await_stop();
    boost::asio::io_service &io_service_;
    boost::asio::ip::tcp::acceptor acceptor_;
    boost::asio::signal_set signals_;
    std::string port_;
    std::set<KmerRequest *> active_;
    KmerPegMapping mapping_;
    boost::asio::ip::tcp::endpoint klookup_endpoint_;
};

#endif
