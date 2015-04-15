
#include "kserver.h"

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <iostream>
#include "global.h"

KmerRequestServer::KmerRequestServer(boost::asio::io_service& io_service,
				     const std::string &port,
				     KmerPegMapping &mapping,
				     boost::asio::ip::tcp::endpoint &klookup_endpoint) :
    io_service_(io_service),
    acceptor_(io_service_),
    port_(port),
    signals_(io_service_),
    mapping_(mapping),
    klookup_endpoint_(klookup_endpoint)
{

    /*
     * Set up for clean signal handling / termination
     */
//    signals_.add(SIGINT);
//    signals_.add(SIGTERM);
//    signals_.add(SIGQUIT);
//    do_await_stop();

    /*
     * Set up listener
     */

    boost::asio::ip::tcp::resolver resolver(io_service_);
    boost::asio::ip::tcp::endpoint endpoint = *resolver.resolve({"0.0.0.0", port_});
    acceptor_.open(endpoint.protocol());
    acceptor_.set_option(boost::asio::ip::tcp::acceptor::reuse_address(true));
    acceptor_.bind(endpoint);
    acceptor_.listen();
    std::cout << "Listening on " << endpoint << "\n";
    do_accept();
}

void KmerRequestServer::do_accept()
{
    std::cout << "start allocate request obj " << g_timer.format();
    KmerRequest *r = new KmerRequest(io_service_, mapping_, klookup_endpoint_);
    std::cout << "finish allocate request obj " << g_timer.format();
    acceptor_.async_accept(r->socket(),
			   boost::bind(&KmerRequestServer::on_accept, this,
				       boost::asio::placeholders::error, r));
}

void KmerRequestServer::on_accept(boost::system::error_code ec, KmerRequest *r)
{
    std::cout << "ACCEPT\n";
    g_timer.start();
    std::cout << "at start " << g_timer.format();
    // Check whether the server was stopped by a signal before this
    // completion handler had a chance to run.
    if (!acceptor_.is_open())
    {
	std::cout << "not open\n";
	return;
    }
    
    if (!ec)
    {
	/*
	 * Connection has come in.
	 * Begin parsing the request line and headers.
	 */
	
	std::cout << "just before do_read " << g_timer.format();
	// active_.insert(r);
	r->do_read();
    }
    
    do_accept();
}

void KmerRequestServer::do_await_stop()
{
    signals_.async_wait(
	[this](boost::system::error_code ec, int signo)
	{
	    std::cout << "Exiting with signal " << signo << "\n";
	    acceptor_.close();
	});
}
