
#include "kserver.h"

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include "global.h"

KmerRequestServer::KmerRequestServer(boost::asio::io_service& io_service,
				     const std::string &port,
				     const std::string &port_file,
				     KmerPegMapping &mapping,
				     boost::asio::ip::tcp::endpoint &klookup_endpoint,
				     std::shared_ptr<KmerGuts> kguts,
				     std::shared_ptr<ThreadPool> thread_pool
    ) :
    io_service_(io_service),
    acceptor_(io_service_),
    port_(port),
    port_file_(port_file),
    signals_(io_service_),
    mapping_(mapping),
    klookup_endpoint_(klookup_endpoint),
    kguts_(kguts),
    thread_pool_(thread_pool)
{
    mapping_map_ = std::make_shared<std::map<std::string, std::shared_ptr<KmerPegMapping> > >();

    auto x = mapping_map_->insert(std::make_pair("", std::make_shared<KmerPegMapping>()));
    std::cerr << "insert created new: " << x.second << " key='" << x.first->first << "'\n";
    x.first->second->load_genome_map((*g_parameters)["kmer-data-dir"].as<std::string>() + "/genomes");

    if (g_parameters->count("reserve-mapping"))
    {
	int count = (*g_parameters)["reserve-mapping"].as<int>();
	std::cerr << "Reserving " << count << " bytes in mapping table\n";
	x.first->second->reserve_mapping_space(count);
    }
}

/*
 * Need to split startup code from constructor due to use of enable_shared_from_this().
 */
void KmerRequestServer::startup()
{

    /*
     * Set up for clean signal handling / termination
     */
    signals_.add(SIGINT);
    signals_.add(SIGTERM);
    signals_.add(SIGQUIT);
    do_await_stop();

    /*
     * Set up listener
     */

    boost::asio::ip::tcp::resolver resolver(io_service_);
    boost::asio::ip::tcp::endpoint endpoint = *resolver.resolve({"0.0.0.0", port_});
    acceptor_.open(endpoint.protocol());
    acceptor_.set_option(boost::asio::ip::tcp::acceptor::reuse_address(true));
    acceptor_.bind(endpoint);
    acceptor_.listen();
    std::cout << "Listening on " << acceptor_.local_endpoint() << "\n";
    if (!port_file_.empty())
    {
	std::ofstream out(port_file_);
	out << acceptor_.local_endpoint().port() << "\n";
	out.close();
    }
	    
    do_accept2();
}

void KmerRequestServer::do_accept()
{
    KmerRequest *r = new KmerRequest(io_service_, mapping_, klookup_endpoint_, kguts_);
    acceptor_.async_accept(r->socket(),
			   boost::bind(&KmerRequestServer::on_accept, this,
				       boost::asio::placeholders::error, r));
}

void KmerRequestServer::on_accept(boost::system::error_code ec, KmerRequest *r)
{
    g_timer.start();
    // Check whether the server was stopped by a signal before this
    // completion handler had a chance to run.
    if (!acceptor_.is_open())
    {
	std::cout << "not open\n";
	delete r;
	return;
    }
    
    if (!ec)
    {
	/*
	 * Connection has come in.
	 * Begin parsing the request line and headers.
	 */
	
	// active_.insert(r);
	r->do_read();
    }
    
    do_accept();
}

void KmerRequestServer::do_accept2()
{
    std::shared_ptr<KmerRequest2> r = std::make_shared<KmerRequest2>(shared_from_this(), io_service_, mapping_map_, kguts_, thread_pool_);
    // std::cerr << "create " << r << "\n";
    acceptor_.async_accept(r->socket(),
			   boost::bind(&KmerRequestServer::on_accept2, this,
				       boost::asio::placeholders::error, r));
    //std::cerr << "leaving do_accept2 r use count=" << r.use_count() << "\n";
}

void KmerRequestServer::on_accept2(boost::system::error_code ec, std::shared_ptr<KmerRequest2> r)
{
    // std::cerr << "on-accept2 use=" << r.use_count() << "\n";
    g_timer.start();
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
	
	active_.insert(r);
	// std::cerr << "start read " << r << "\n";
	r->do_read();
    }
    
    do_accept2();
}

void KmerRequestServer::deactivate(std::shared_ptr<KmerRequest2> x)
{
    active_.erase(x);
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
