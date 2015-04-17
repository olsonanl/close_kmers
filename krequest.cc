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
    klookup_(0),
    klookup2_(0),
    response_stream_(&response_)
{
}

KmerRequest::~KmerRequest()
{
    if (krequest_)
	delete krequest_;
    if (klookup_)
	delete klookup_;
    if (klookup2_)
	delete klookup2_;
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
	    }
	    else
	    {
		size_t x = line.find(':');
		std::string k(line.substr(0,x));
		x++;
		while (line[x] == ' ')
		    x++;
		std::string v(line.substr(x));
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
	auto kv = headers_.find("content-length");
	if (kv != headers_.end())
	{
	    int content_size = std::stoi(kv->second);
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
	// std::cout << "got EOF\n";
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

    if (path_ == "/lookup")
    {
	krequest_ = new std::istream(&request_);
	klookup_ = new KmerLookupClient(io_service_, klookup_endpoint_, *krequest_, mapping_,
					boost::bind(&KmerRequest::request_complete, this, _1));
    }
    else if (path_ == "/add")
    {
	response_stream_ << "HTTP/1.1 200 OK\n";
	response_stream_ << "Content-type: text/plain\n";
	response_stream_ << "\n";

	krequest_ = new std::istream(&request_);
	klookup2_ = new KmerLookupClient2(io_service_, klookup_endpoint_, *krequest_,
					  boost::bind(&KmerRequest::on_protein, this, _1),
					  boost::bind(&KmerRequest::on_hit, this, _1),
					  boost::bind(&KmerRequest::on_call, this, _1, _2),
					  boost::bind(&KmerRequest::add_complete, this, _1));
					
    }
    else
    {
	std::cout << "not found " << path_ << "\n";

	response_stream_ << "HTTP/1.1 404 not found\n";
	response_stream_ << "Content-type: text/plain\n";
	response_stream_ << "\n";

	boost::asio::async_write(socket_, response_,
				 boost::bind(&KmerRequest::write_response_complete, this,
					     boost::asio::placeholders::error));

    }
}

void KmerRequest::on_protein(const std::string &protein)
{
    cur_protein_ = protein;
    cur_protein_id_ = mapping_.encode_id(protein);
    // std::cout << "on protein " << protein << "\n";
}

void KmerRequest::on_call(const std::string &function, const std::string &count)
{
    // std::cout << "call " << function << "\n";
    response_stream_ << cur_protein_ << "\t" << function << "\t" << count << "\n";
}

void KmerRequest::on_hit(unsigned long kmer)
{
    mapping_.add_mapping(cur_protein_id_, kmer);
}

void KmerRequest::add_complete( const boost::system::error_code& err )
{
    // std::cout << "Got add complete\n";

    boost::asio::async_write(socket_, response_,
			     boost::bind(&KmerRequest::write_response_complete, this,
					 boost::asio::placeholders::error));

    delete krequest_;
    delete klookup2_;
    krequest_ = 0;
    klookup2_ = 0;
}


void KmerRequest::request_complete( const KmerLookupClient::result_t &resp)
{
    response_stream_ << "HTTP/1.1 200 OK\n";
    response_stream_ << "Content-type: text/plain\n";
    response_stream_ << "\n";
    
    // Write response here and we're done but for cleanup
    for (auto it = resp.begin(); it != resp.end(); it++)
    {
	response_stream_ << it->first << "\t" << it->second << "\n";
    }
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
    socket_.close();

    delete this;
}
