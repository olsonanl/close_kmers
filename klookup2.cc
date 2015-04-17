
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
				     std::istream &input,
				     boost::function<void ( const std::string &prot )> on_protein,
				     boost::function<void ( unsigned long kmer )> on_hit,
				     boost::function<void ( const std::string &function, const std::string &count )> on_call,
				     boost::function<void ( const boost::system::error_code& err )> on_completion)
    : resolver_(io_service),
      socket_(io_service),
      input_(input),
      on_protein_(on_protein),
      on_hit_(on_hit),
      on_call_(on_call),
      on_completion_(on_completion)
{
    std::ostream request_stream(&request_);
    request_stream << "-d 1 -a\n";

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
	    
	}
	
	//
	// Read next block from the file and initiate write.
	//
	input_.read(buffer_, sizeof(buffer_));
	
	if (input_.gcount() > 0)
	{
	    // std::cout << "write " << input_.gcount() << "\n";
	    
	    boost::asio::async_write(socket_, boost::asio::buffer(buffer_, input_.gcount()),
				     boost::bind(&KmerLookupClient2::handle_write_request, this,
						 boost::asio::placeholders::error));
	}
    }
    else
    {
	std::cout << "Error: " << err.message() << "\n";
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

		on_hit_(kmer);
	    }
	    else if (h == "PROTEIN-ID")
	    {
		std::string prot;
		s >> prot;
		on_protein_(prot);
	    }
	    else if (h == "CALL")
	    {
	        // CALL    0       226     184     17583   RNA methyltransferase, TrmH family, group 1     389.544617
		
		std::string word, count, function;
		std::getline(s, word, '\t');
		std::getline(s, word, '\t');
		std::getline(s, word, '\t');
		std::getline(s, count, '\t');
		std::getline(s, word, '\t');
		std::getline(s, function, '\t');
		on_call_(function, count);
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
