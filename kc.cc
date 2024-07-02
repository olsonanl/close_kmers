#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind/bind.hpp>
#include <boost/filesystem.hpp>

#include "kmer.h"
#include "kserver.h"

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

class client
{
public:
    client(boost::asio::io_service& io_service,
	   const std::string& server, const std::string &port, std::istream &input,
	   KmerPegMapping &mapping)
        : resolver_(io_service),
	  socket_(io_service),
	  input_(input),
	  mapping_(mapping)
	{
	    std::ostream request_stream(&request_);
	    request_stream << "-d 1 -a\n";

	    // Start an asynchronous resolve to translate the server and service names
	    // into a list of endpoints.
	    tcp::resolver::query query(server, port);
	    resolver_.async_resolve(query,
				    boost::bind(&client::handle_resolve, this,
						boost::asio::placeholders::error,
						boost::asio::placeholders::iterator));
	}

private:
    void handle_resolve(const boost::system::error_code& err,
			tcp::resolver::iterator endpoint_iterator)
	{
	    if (!err)
	    {
		// Attempt a connection to each endpoint in the list until we
		// successfully establish a connection.
		boost::asio::async_connect(socket_, endpoint_iterator,
					   boost::bind(&client::handle_connect, this,
						       boost::asio::placeholders::error));
	    }
	    else
	    {
		std::cout << "Error: " << err.message() << "\n";
	    }
	}

    void handle_connect(const boost::system::error_code& err)
	{
	    if (!err)
	    {
		// The connection was successful. Send the request.
		boost::asio::async_write(socket_, request_,
					 boost::bind(&client::handle_write_request, this,
						     boost::asio::placeholders::error));

		//
		// Also kick off a reader.
		//
		boost::asio::async_read_until(socket_, response_, "\n",
					boost::bind(&client::handle_read, this,
						    boost::asio::placeholders::error,
						    boost::asio::placeholders::bytes_transferred));


	    }
	    else
	    {
		std::cout << "Error: " << err.message() << "\n";
	    }
	}

    void handle_write_request(const boost::system::error_code& err)
	{
	    if (!err)
	    {
		//
		// check for EOF from last time
		//
		if (!input_)
		{
		    std::cout << "all done reading\n";
		    
		    std::ostream request_stream(&request_);
		    request_stream << ">FLUSH\n";
		    
		    boost::asio::async_write(socket_, request_,
					     boost::bind(&client::finish_write_request, this,
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
					     boost::bind(&client::handle_write_request, this,
							 boost::asio::placeholders::error));
		}
	    }
	    else
	    {
		std::cout << "Error: " << err.message() << "\n";
	    }
	}

    void finish_write_request(const boost::system::error_code& err)
	{
	    if (!err)
	    {
		std::cout << "done writing\n";
	    }
	    else
	    {
		std::cout << "Error: " << err.message() << "\n";
	    }
	}

    void handle_read(const boost::system::error_code& err, size_t bytes_transferred)
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
				//std::cout << "  " << *it << "\n";
				
				hit_count_[*it]++;
			    }
			}
		    }
		    else if (h == "OTU-COUNTS")
		    {
			typedef std::pair<KmerPegMapping::encoded_id_t, unsigned int> data_t;
			std::vector< data_t > vec(hit_count_.begin(), hit_count_.end());
			std::sort(vec.begin(), vec.end(), less_second<data_t>()); 

			for (auto it = vec.begin(); it != vec.end(); it++)
			{
			    std::string peg = mapping_.decode_id(it->first);
			    std::cout << peg << ": " << it->second  << "\n";
			}

			socket_.close();
			return;
		    }
		    else
		    {
			std::cout << line << "\n";
		    }
		}

		boost::asio::async_read_until(socket_, response_, "\n",
					boost::bind(&client::handle_read, this,
						    boost::asio::placeholders::error,
						    boost::asio::placeholders::bytes_transferred));
	    }
	    else
	    {
		std::cout << "Error: " << err << "\n";
	    }
	}


    tcp::resolver resolver_;
    tcp::socket socket_;
    boost::asio::streambuf request_;
    boost::asio::streambuf response_;
    char buffer_[4096];
    std::istream &input_;
    KmerPegMapping &mapping_;
    std::map<KmerPegMapping::encoded_id_t, unsigned int> hit_count_;
};

int main(int argc, char* argv[])
{
    //KmerPegMapping mapping("/scratch/olson/core.kmers.2015-0406/Data.2");
    KmerPegMapping mapping("/Users/olson/kmer/Data.2");

    //path data_dir("/scratch/olson/core.kmers.2015-0406/peg.kmers");
    path data_dir("/Users/olson/kmer/xx");

    if (is_directory(data_dir))
    {
	directory_iterator end_iter;
	for (auto dit = directory_iterator(data_dir); dit != end_iter; dit++)
	{
	    std::cout << dit->path().string() << "\n";
	    mapping.load_compact_mapping_file(dit->path().string());
	}
    }
    else
    {
	mapping.load_compact_mapping_file(data_dir.string());
    }

      try
      {
	  if (argc != 4)
	  {
	      std::cout << "Usage: async_client <server> <path>\n";
	      std::cout << "Example:\n";
	      std::cout << "  async_client www.boost.org /LICENSE_1_0.txt\n";
	      return 1;
	  }

	  boost::asio::io_service io_service;

	  KmerRequestServer kserver(io_service, "5100");

	  std::ifstream ifile(argv[3]);
	  
	  client c(io_service, argv[1], argv[2], ifile, mapping);
	  io_service.run();
      }
      catch (std::exception& e)
      {
	  std::cout << "Exception: " << e.what() << "\n";
      }

      return 0;
}
