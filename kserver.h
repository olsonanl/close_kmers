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
		      const std::string &port_file,
		      KmerPegMapping &mapping,
		      boost::asio::ip::tcp::endpoint &klookup_endpoint,
		      std::shared_ptr<KmerGuts> kguts);

private:

    void do_accept();
    void on_accept(boost::system::error_code ec, KmerRequest *);
    
    void do_await_stop();
    std::shared_ptr<KmerGuts> kguts_;
    boost::asio::io_service &io_service_;
    boost::asio::ip::tcp::acceptor acceptor_;
    boost::asio::signal_set signals_;
    std::string port_;
    std::string port_file_;
    std::set<KmerRequest *> active_;
    KmerPegMapping &mapping_;
    boost::asio::ip::tcp::endpoint klookup_endpoint_;
};

#endif
