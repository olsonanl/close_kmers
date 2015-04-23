
#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include "klookup3.h"
#include "global.h"

using namespace boost::filesystem;
using boost::asio::ip::tcp;

KmerLookupClient3::KmerLookupClient3(boost::asio::io_service& io_service,
				     boost::asio::ip::tcp::endpoint endpoint,
				     KmerLookupClient3::stream_queue_t &stream_queue,
				     boost::function<void ( const std::string &prot )> on_protein,
				     boost::function<void ( unsigned long kmer )> on_hit,
				     boost::function<void ( const std::string &line )> on_call,
				     boost::function<void ( const boost::system::error_code& err )> on_completion)
    : resolver_(io_service),
      socket_(io_service),
      on_protein_(on_protein),
      on_hit_(on_hit),
      on_call_(on_call),
      on_completion_(on_completion),
      stream_queue_(stream_queue),
      write_pending_(0)
{
    timer_.start();
    
    boost::asio::async_connect(socket_, &endpoint,
			       boost::bind(&KmerLookupClient3::handle_connect, this,
					   boost::asio::placeholders::error));

}

void KmerLookupClient3::handle_connect(const boost::system::error_code& err)
{
    if (!err)
    {
	// check_queue();
	//
	// Also kick off a reader.
	//
	boost::asio::async_read_until(socket_, response_, "\n",
				      boost::bind(&KmerLookupClient3::handle_read, this,
						  boost::asio::placeholders::error,
						  boost::asio::placeholders::bytes_transferred));
	
	
    }
    else
    {
	std::cout << "Error: " << err.message() << "\n";
	on_completion_(err);
    }
}

/*
 * Check for the presence of valid at at the front of our input queue.
 * if there is some, initiate a write.
 *
 * If a write is pending, do not do anything.
 */
void KmerLookupClient3::check_queue()
{
    std::cout << "check_queue " << stream_queue_.size() << "\n";
    if (write_pending_)
    {
	std::cout << "check_queue invoked with a write pending\n";
	return;
    }
    
    while (stream_queue_.size() > 0)
    {
	boost::shared_ptr<std::istream> &input = stream_queue_.front();
	if (!input)
	{
	    std::cout << "check_queue found null\n";
	    boost::system::error_code ec;
	    socket_.shutdown(boost::asio::ip::tcp::socket::shutdown_send, ec);
	    return;
	}
	else
	{
	    if (*input)
	    {
		std::cout << "check found data available\n";

		//
		// Read next block from the file and initiate write.
		//
		input->read(buffer_, sizeof(buffer_));
	
		if (input->gcount() > 0)
		{
		    std::cout << "write " << input->gcount() << "\n";

		    write_pending_ = 1;
		    boost::asio::async_write(socket_, boost::asio::buffer(buffer_, input->gcount()),
					     boost::bind(&KmerLookupClient3::handle_write_request, this,
							 boost::asio::placeholders::error));
		    return;
		}

		if (!*input)
		{
		    std::cout << "Hit EOF on input\n";
		    stream_queue_.pop_front();
		}
	    }
	    else
	    {
		std::cout << "Came in with EOF on input\n";
		stream_queue_.pop_front();
	    }

	}
    }
}

void KmerLookupClient3::handle_write_request(const boost::system::error_code& err)
{
    if (!err)
    {
	std::cout << "write completed\n";
	write_pending_ = 0;
	check_queue();
    }
    else
    {
	std::cout << "Error: " << err.message() << " at " << __FILE__ << ":" << __LINE__ << "\n";
	on_completion_(err);
    }
}


void KmerLookupClient3::handle_read(const boost::system::error_code& err, size_t bytes_transferred)
{
    if (!err)
    {
	// std::cout << "Handle read " << bytes_transferred << ":\n";

	std::istream resp(&response_);
	std::string line;
	if (std::getline(resp, line, '\n'))
	{

	    // std::cout << line << "\n";
	    
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
		if (on_call_)
		    on_call_(line);
	    }
	    else if (h == "PROTEIN-ID")
	    {
		std::string prot;
		s >> prot;
		if (on_protein_)
		    on_protein_(prot);
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
				      boost::bind(&KmerLookupClient3::handle_read, this,
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
