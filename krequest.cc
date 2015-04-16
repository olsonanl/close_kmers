#include "krequest.h"
#include <iostream>
#include <boost/bind.hpp>
#include <algorithm>

#include "klookup.h"
#include "global.h"

KmerRequest::KmerRequest(boost::asio::io_service &io_service,
			 KmerPegMapping &mapping,
			 boost::asio::ip::tcp::endpoint &klookup_endpoint) :
    io_service_(io_service),
    socket_(io_service_),
    mapping_(mapping),
    klookup_endpoint_(klookup_endpoint),
    krequest_(0),
    klookup_(0)
{
}

KmerRequest::~KmerRequest()
{
    if (krequest_)
	delete krequest_;
    if (klookup_)
	delete klookup_;
}

/*
 * Start processing our request.
 */
void KmerRequest::do_read()
{
    timer_.start();

    boost::asio::async_read_until(socket_, request_, "\n",
				  boost::bind(&KmerRequest::handle_read, this,
					      boost::asio::placeholders::error,
					      boost::asio::placeholders::bytes_transferred));
}

void KmerRequest::handle_read(boost::system::error_code err, size_t bytes)
{
    if (!err)
    {
	std::cout << "krequest handle_read " << g_timer.format();
	std::istream resp(&request_);
	std::string line;
	int done = 0;
	if (std::getline(resp, line, '\n'))
	{
	    
	    size_t l = line.length();
	    if (l && line[l - 1] == '\r')
		line.pop_back();

	    // std::cout << "|" << line << "|\n";
	    
	    std::stringstream ss(line);
	    if (line.length() == 0)
	    {
		done = 1;
	    }
	    else if (request_type_.empty())
	    {
		ss >> request_type_;
		ss >> path_;
		std::cout << path_ << " at " << request_type_ << "\n";
		std::cout << "time at hdr recv " << g_timer.format();
	    }
	    else
	    {
		size_t x = line.find(':');
		std::string k(line.substr(0,x));
		x++;
		while (line[x] == ' ')
		    x++;
		std::string v(line.substr(x));
		std::cout << k << " = " << v << "\n";
		std::transform(k.begin(), k.end(), k.begin(), ::tolower);
		headers_[k] = v;
	    }
	}

	if (done)
	{
	    handle_request();
	}
	else
	{
	    boost::asio::async_read_until(socket_, request_, "\n",
					  boost::bind(&KmerRequest::handle_read, this,
						      boost::asio::placeholders::error,
						      boost::asio::placeholders::bytes_transferred));
	}
    }
    else
    {
	std::cout << "error is " << err << "\n";
    }		      
}

void KmerRequest::handle_request()
{
    if (request_type_ == "GET")
    {
	// send_response("kmer service\n");
    }
    else if (request_type_ == "POST")
    {
	boost::asio::async_read(socket_, request_,
				    boost::asio::transfer_at_least(1),
				    boost::bind(&KmerRequest::handle_post_body, this,
						boost::asio::placeholders::error,
						boost::asio::placeholders::bytes_transferred));
    }
}

void KmerRequest::handle_post_body(boost::system::error_code err, size_t bytes)
{
    if (!err)
    {
	std::cout << "post body\n";
	std::cout << "size='" << request_.size() << "'\n";

	auto kv = headers_.find("content-length");
	if (kv != headers_.end())
	{
	    int content_size = std::stoi(kv->second);
	    std::cout << "content size=" << content_size << "\n";
	    std::cout << "size='" << request_.size() << "'\n";
	    if (request_.size() < content_size)
	    {
		boost::asio::async_read(socket_, request_,
					boost::asio::transfer_at_least(1),
					boost::bind(&KmerRequest::handle_post_body, this,
						    boost::asio::placeholders::error,
						    boost::asio::placeholders::bytes_transferred));
	    }
	    else
	    {
		process_request();
	    }
	}
	
//	std::string s((std::istreambuf_iterator<char>(&request_)),
//		      std::istreambuf_iterator<char>());
//	std::cout  << s << "\n";;
    }
    else if (err == boost::asio::error::eof)
    {
	std::cout << "got EOF\n";
    }
    else
    {
	std::cout << "got error " << err << "\n";
    }
}


void KmerRequest::process_request()
{
    /*
     * We have our request. Create a stream buffer to pass to the lookup client
     * and process the lookup.
     */

    std::cout << "initiate lookup " << g_timer.format();
    krequest_ = new std::istream(&request_);
    klookup_ = new KmerLookupClient(io_service_, klookup_endpoint_, *krequest_, mapping_,
				    boost::bind(&KmerRequest::request_complete, this, _1));

}

void KmerRequest::request_complete( const KmerLookupClient::result_t &resp)
{
    std::cout << "Got response\n";

    std::ostream resp_stream(&response_);

    resp_stream << "HTTP/1.1 200 OK\n";
    resp_stream << "Content-type: text/plain\n";
    resp_stream << "\n";
    
    // Write response here and we're done but for cleanup
    for (auto it = resp.begin(); it != resp.end(); it++)
    {
	resp_stream << it->first << "\t" << it->second << "\n";
    }
    std::cout << "Starting writing response " << g_timer.format() << "\n";

    boost::asio::async_write(socket_, response_,
			     boost::bind(&KmerRequest::write_response_complete, this,
					 boost::asio::placeholders::error));

    delete krequest_;
    delete klookup_;
    krequest_ = 0;
    klookup_ = 0;
}

void KmerRequest::write_response_complete(boost::system::error_code err)
{
    std::cout << "Finished writing response " << g_timer.format() << "\n";

    socket_.close();

    delete this;
}
