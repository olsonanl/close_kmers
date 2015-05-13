
#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include "klookup2.h"
#include "global.h"

using namespace boost::filesystem;
using boost::asio::ip::tcp;

KmerLookupClient2::KmerLookupClient2(boost::asio::io_service& io_service,
				     boost::asio::ip::tcp::endpoint endpoint,
				     const std::string &kmer_options,
				     std::istream &input,
				     boost::function<void ( const std::string &prot, size_t len )> on_protein,
				     boost::function<void ( unsigned long kmer )> on_hit,
				     boost::function<void ( const std::string &line )> on_call,
				     boost::function<void ( const boost::system::error_code& err )> on_completion)
    : resolver_(io_service),
      socket_(io_service),
      kmer_options_(kmer_options),
      input_(input),
      on_protein_(on_protein),
      on_hit_(on_hit),
      on_call_(on_call),
      on_completion_(on_completion)
{
    std::ostream request_stream(&request_);
    if (!kmer_options_.empty())
    {
	request_stream << kmer_options_ << "\n";
    }

    timer_.start();
    
    boost::asio::async_connect(socket_, &endpoint,
			       boost::bind(&KmerLookupClient2::handle_connect, this,
					   boost::asio::placeholders::error));

}

void KmerLookupClient2::handle_connect(const boost::system::error_code& err)
{
    if (!err)
    {

	// The connection was successful. Send the request.
	boost::asio::async_write(socket_, request_,
				 boost::bind(&KmerLookupClient2::handle_write_request, this,
					     boost::asio::placeholders::error));
	
	//
	// Also kick off a reader.
	//
	boost::asio::async_read_until(socket_, response_, "\n",
				      boost::bind(&KmerLookupClient2::handle_read, this,
						  boost::asio::placeholders::error,
						  boost::asio::placeholders::bytes_transferred));
	
	
    }
    else
    {
	std::cout << "Error: " << err.message() << "\n";
	on_completion_(err);
    }
}

void KmerLookupClient2::handle_write_request(const boost::system::error_code& err)
{
    if (!err)
    {
	//
	// check for EOF from last time
	//
	if (!input_)
	{
	    // std::cout << "all done reading\n";
	    // std::cout << "elapsed at read " << g_timer.format();
	    
	    std::ostream request_stream(&request_);
	    request_stream << ">FLUSH\n";
	    
	    boost::asio::async_write(socket_, request_,
				     boost::bind(&KmerLookupClient2::finish_write_request, this,
						 boost::asio::placeholders::error));
	    return;
	}
	
	//
	// Read next block from the file and initiate write.
	//
	input_.read(buffer_, sizeof(buffer_));
	
	if (input_.gcount() > 0)
	{
	    std::cout << "write " << input_.gcount() << "\n";
	    
	    boost::asio::async_write(socket_, boost::asio::buffer(buffer_, input_.gcount()),
				     boost::bind(&KmerLookupClient2::handle_write_request, this,
						 boost::asio::placeholders::error));
	}
	else
	{
	    
	    std::cout << "handle_write_request: read buffer return gcount " << input_.gcount() << "\n";

	    if (!input_)
	    {
		std::cout << "EOF hit\n";
		std::ostream request_stream(&request_);
		request_stream << ">FLUSH\n";
		
		boost::asio::async_write(socket_, request_,
					 boost::bind(&KmerLookupClient2::finish_write_request, this,
						     boost::asio::placeholders::error));
	    }
	    else
	    {
		std::cout << "Not eof? \n";
	    }
	       
	}
    }
    else
    {
	std::cout << "Error: " << err.message() << " at " << __FILE__ << ":" << __LINE__ << "\n";
	on_completion_(err);
    }
}

void KmerLookupClient2::finish_write_request(const boost::system::error_code& err)
{
    if (!err)
    {
	boost::system::error_code ec;
	socket_.shutdown(boost::asio::ip::tcp::socket::shutdown_send, ec);
    }
    else
    {
	std::cout << "Error: " << err.message() << "\n";
    }
}

void KmerLookupClient2::handle_read(const boost::system::error_code& err, size_t bytes_transferred)
{
    if (!err)
    {
//	std::cout << "Handle read " << bytes_transferred << ":\n";

	std::istream resp(&response_);
	std::string line;
	if (std::getline(resp, line, '\n'))
	{

//	    std::cout << line << "\n";
	    
	    std::string h;
	    unsigned long offset, kmer;
	    std::stringstream s(line);
	    s >> h;
	    
	    if (h == "HIT")
	    {
		s >> offset;
		s >> kmer;

		if (on_hit_)
		    on_hit_(kmer);
//		if (on_call_)
//		    on_call_(line);
	    }
	    else if (h == "PROTEIN-ID")
	    {
		std::string prot;
		size_t len;
		s >> prot;
		s >> len;
		if (on_protein_)
		    on_protein_(prot, len);
		if (on_call_)
		    on_call_(line);
	    }
	    else if (h == "CALL")
	    {
	        // CALL    0       226     184     17583   RNA methyltransferase, TrmH family, group 1     389.544617

		if (on_call_)
		    on_call_(line);
		
	    }
	    else
	    {
		if (on_call_)
		    on_call_(line);
	    }		
	}
	boost::asio::async_read_until(socket_, response_, "\n",
				      boost::bind(&KmerLookupClient2::handle_read, this,
						  boost::asio::placeholders::error,
						  boost::asio::placeholders::bytes_transferred));
    }
    else if (err == boost::asio::error::eof)
    {
	// std::cout << "EOF\n";
	boost::system::error_code ec;
	on_completion_(ec);
    }
    else
    {
	std::cout << "ErrorX: " << err << "\n";
	on_completion_(err);
    }
};
