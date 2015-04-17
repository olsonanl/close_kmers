
#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include "klookup.h"
#include "global.h"

using namespace boost::filesystem;
using boost::asio::ip::tcp;

template<class T>
struct less_second : std::binary_function<T,T,bool>
{
    inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
};  

KmerLookupClient::KmerLookupClient(boost::asio::io_service& io_service,
				   boost::asio::ip::tcp::endpoint endpoint,
				   std::istream &input,
				   KmerPegMapping &mapping,
				   boost::function<void ( const result_t &)> on_completion)
    : resolver_(io_service),
      socket_(io_service),
      input_(input),
      mapping_(mapping),
      on_completion_(on_completion)
{
    std::ostream request_stream(&request_);
    request_stream << "-d 1 -a\n";

    timer_.start();
    
    boost::asio::async_connect(socket_, &endpoint,
			       boost::bind(&KmerLookupClient::handle_connect, this,
					   boost::asio::placeholders::error));
}

void KmerLookupClient::handle_connect(const boost::system::error_code& err)
{
    if (!err)
    {

	// The connection was successful. Send the request.
	boost::asio::async_write(socket_, request_,
				 boost::bind(&KmerLookupClient::handle_write_request, this,
					     boost::asio::placeholders::error));
	
	//
	// Also kick off a reader.
	//
	boost::asio::async_read_until(socket_, response_, "\n",
				      boost::bind(&KmerLookupClient::handle_read, this,
						  boost::asio::placeholders::error,
						  boost::asio::placeholders::bytes_transferred));
	
	
    }
    else
    {
	std::cout << "Error: " << err.message() << "\n";
    }
}

void KmerLookupClient::handle_write_request(const boost::system::error_code& err)
{
    if (!err)
    {
	//
	// check for EOF from last time
	//
	if (!input_)
	{
	    std::ostream request_stream(&request_);
	    request_stream << ">FLUSH\n";
	    
	    boost::asio::async_write(socket_, request_,
				     boost::bind(&KmerLookupClient::finish_write_request, this,
						 boost::asio::placeholders::error));
	    
	}
	
	//
	// Read next block from the file and initiate write.
	//
	input_.read(buffer_, sizeof(buffer_));
	
	if (input_.gcount() > 0)
	{
	    std::cout << "write " << input_.gcount() << "\n";
	    
	    boost::asio::async_write(socket_, boost::asio::buffer(buffer_, input_.gcount()),
				     boost::bind(&KmerLookupClient::handle_write_request, this,
						 boost::asio::placeholders::error));
	}
    }
    else
    {
	std::cout << "Error: " << err.message() << "\n";
    }
}

void KmerLookupClient::finish_write_request(const boost::system::error_code& err)
{
    if (!err)
    {
    }
    else
    {
	std::cout << "Error: " << err.message() << "\n";
    }
}

void KmerLookupClient::handle_read(const boost::system::error_code& err, size_t bytes_transferred)
{
    if (!err)
    {
//		std::cout << "Handle read " << bytes_transferred << ":\n";

	std::istream resp(&response_);
	std::string line;
	if (std::getline(resp, line, '\n'))
	{
	    std::string h;
	    unsigned long offset, kmer;
	    std::stringstream s(line);
	    s >> h;
	    
	    if (h == "HIT")
	    {
		s >> offset;
		s >> kmer;
		
		auto ki = mapping_.kmer_to_id_.find(kmer);
		if (ki != mapping_.kmer_to_id_.end())
		{
		    // std::cout << "got mapping for " << kmer << "\n";
		    for (auto it = ki->second.begin(); it != ki->second.end(); it++)
		    {
			// std::cout << "  " << *it << "\n";
			
			hit_count_[*it]++;
		    }
		}
	    }
	    else if (h == "OTU-COUNTS")
	    {
		typedef std::pair<std::string, unsigned int> data_t;

		std::vector<data_t> vec;
		for (auto it = hit_count_.begin(); it != hit_count_.end(); it++)
		{
		    std::string peg = mapping_.decode_id(it->first);
		    vec.push_back(data_t(peg, it->second));
		}

		std::sort(vec.begin(), vec.end(), less_second<data_t>()); 
		
		for (auto it = vec.begin(); it != vec.end(); it++)
		{
		    std::cout << it->first << ": " << it->second  << "\n";
		}

		std::cout << "at receipt of finished results " << g_timer.format();
		socket_.close();
		on_completion_(vec);
		return;
	    }
	    else
	    {
		std::cout << line << "\n";
	    }
	}
	boost::asio::async_read_until(socket_, response_, "\n",
				      boost::bind(&KmerLookupClient::handle_read, this,
						  boost::asio::placeholders::error,
						  boost::asio::placeholders::bytes_transferred));
    }
    else
    {
	std::cout << "Error: " << err << "\n";
    }
};
