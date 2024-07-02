#include "krequest.h"
#include <iostream>
#include <sstream>
#include <boost/bind/bind.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>

#include "klookup.h"
#include "kguts.h"
#include "global.h"

const boost::regex request_regex("^([A-Z]+) ([^?#]*)(\\?([^#]*))?(#(.*))? HTTP/(\\d+\\.\\d+)");
/* $1 = type
 * $2 = path
 * $4 = parameters
 * $6 = fragment
 */

KmerRequest::KmerRequest(boost::asio::io_service &io_service,
			 KmerPegMapping &mapping,
			 boost::asio::ip::tcp::endpoint &klookup_endpoint,
			 std::shared_ptr<KmerGuts> kguts) :
    io_service_(io_service),
    socket_(io_service_),
    mapping_(mapping),
    klookup_endpoint_(klookup_endpoint),
    krequest_(0),
    klookup_(0),
    klookup2_(0),
    response_stream_(&response_),
    bytes_left_(0),
    kguts_(kguts)
{
}

KmerRequest::~KmerRequest()
{
    // std::cout << "DESTROY KmerRequest\n";
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
	boost::asio::streambuf::const_buffers_type bufs = request_.data();
	int done = 0;

	std::string line(buffers_begin(bufs), buffers_begin(bufs) + bytes);
	request_.consume(bytes);
	
	{
	    size_t l = line.length();
	    if (l && line[l - 1] == '\n')
		line.pop_back();
	    l = line.length();
	    if (l && line[l - 1] == '\r')
		line.pop_back();

	    // std::cout << "|" << line << "| " << line.length() << "\n";
	    
	    std::stringstream ss(line);
	    if (line.length() == 0)
	    {
		// std::cout <<"done\n";
		done = 1;
	    }
	    else if (request_type_.empty())
	    {
		boost::smatch match;
		if (boost::regex_match(line, match, request_regex))
		{
		    request_type_ = match[1];
		    path_ = match[2];
		    parameters_raw_ = match[4];
		    fragment_ = match[6];

		    if (!parameters_raw_.empty())
		    {
			std::vector<std::string> parts;
			boost::split(parts, parameters_raw_, boost::is_any_of(";&"));
			for (auto it : parts)
			{
			    size_t pos = it.find('=');
			    if (pos != std::string::npos)
			    {
				parameters_[it.substr(0, pos)] = it.substr(pos+1);
			    }
			}
			/*
			for (auto it : parameters_)
			{
			    std::cout << "'" << it.first << "': '" << it.second << "'\n";
			}
			*/
		    }
		    // std::cout << "req=" << request_type_ << " path=" << path_ << "\n";
		}
		else
		{
		    std::cout << "Invalid request '" << line << "'\n";
		    socket_.close();
		    delete this;
		}
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
		// std::cout << "'" << k << "': '" << v << "'\n";
	    }
	}

	if (done)
	{
	    auto it = headers_.find("expect");
	    if (it != headers_.end() && it->second == "100-continue")
	    {
		std::ostream os(&continue_response_);
		os << "HTTP/1.1 100 Continue\n\n";
		boost::asio::async_write(socket_, continue_response_,
					 boost::bind(&KmerRequest::handle_request, this));
	    }
	    else
	    {
		handle_request();
	    }
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
	if (path_ == "/quit")
	{
	    socket_.close();
	    io_service_.stop();
	}
	if (path_ == "/dump_mapping")
	{
	    response_stream_ << "HTTP/1.1 200 OK\n";
	    response_stream_ << "Content-type: text/plain\n";
	    response_stream_ << "\n";
	    for (auto it: mapping_.kmer_to_id_)
	    {
		response_stream_ << it.first << "\t";
		for (auto elt: it.second)
		{
		    response_stream_ << " " << elt;
		}
		response_stream_ << "\n";
	    }
	    for (auto it: mapping_.family_mapping_)
	    {
		response_stream_ << it.first << "\t" << it.second.pgf << "\t" << it.second.plf << "\t" << it.second.function << "\n";
	    }
	}
	else
	{
	    response_stream_ << "HTTP/1.1 404 not found\n";
	    response_stream_ << "Content-type: text/plain\n";
	    response_stream_ << "\n";
	}
	    
	boost::asio::async_write(socket_, response_,
				 boost::bind(&KmerRequest::write_response_complete, this,
					     boost::asio::placeholders::error));
	return;
    }
    else if (request_type_ == "POST")
    {
	if (path_ == "/matrix")
	{
	    std::string opts("");
	    auto it = headers_.find("kmer-options");
	    if (it != headers_.end())
	    {
		opts = it->second;

		boost::shared_ptr<std::istream> sp(new std::istringstream(opts + "\n"));
		stream_queue_.push_back(sp);
	    }

	    boost::shared_ptr<KmerLookupClient3>
		kc(new KmerLookupClient3(io_service_,
					 klookup_endpoint_,
					 stream_queue_,
					 boost::bind(&KmerRequest::on_matrix_protein, this, _1, _2),
					 boost::bind(&KmerRequest::on_matrix_hit, this, _1),
					 0,
					 boost::bind(&KmerRequest::matrix_complete, this, _1),
					 parameters_));
	    klookup3_ = kc;
    
	    auto kv = headers_.find("content-length");
	    if (kv != headers_.end())
	    {
		size_t content_size = std::stoul(kv->second);
// 		std::cout << "reading " << content_size << "\n";
		
		bytes_left_ = content_size;
		cur_sequence_ = "";
		boost::system::error_code ec;
		handle_matrix(ec, 0);
	    }
	}
	else
	{
	    // std::cout << "handling post buf size " << request_.size() << "\n";
	    auto kv = headers_.find("content-length");
	    if (kv != headers_.end())
	    {
		size_t content_size = std::stoul(kv->second);
		// std::cout << "reading " << content_size << "\n";
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
	    else
	    {
		std::cout << "no content length provided\n";
		socket_.close();
		delete this;
	    }
	}
    }
}

/*
 * Process a matrix distance compute request. Input is a set of sequences
 * in fasta format.
 * The plan is to parse line by line, building up a list of strings one per
 * sequence.
 */
void KmerRequest::handle_matrix(boost::system::error_code err, size_t bytes)
{
    if (err)
    {
	std::cout << "handle_matrix() Err" << err << "\n";
	socket_.close();
	delete this;
	return;
    }
    std::string cur;
    std::istream istr(&request_);

    // std::cout << "req size is " << request_.size() << " "  << bytes << " " << bytes_left_ << "\n";
    while (request_.size() > 0 && istr && !istr.fail())
    {
	getline(istr, cur);
	bytes_left_ -= cur.size();
    
	if (istr.rdstate() && std::ios_base::badbit)
	{
	    // std::cout << "got partial line " << istr.rdstate() << " '" << cur<< "'\n";
	    cur_sequence_ += cur;
	    break;
	}
	bytes_left_--;
	// std::cout << "got line " << istr.rdstate() << " '" << cur<< "'\n";
	// request_.consume(cur.size());

	if (cur[0] == '>')
	{
	    if (!cur_sequence_.empty())
	    {
		handle_matrix_sequence(cur_sequence_);
	    }
	    cur_sequence_ = cur + "\n";
	}
	else
	{
	    cur_sequence_ += cur + "\n";
	}
	
	if (!istr)
	{
	    std::cout << "not istr!\n";
	    return;
	}
	else if (istr.rdstate() & std::ios_base::failbit)
	{
	    std::cout << "got failbit\n";
	}
    }
    if (bytes_left_)
    {
	// std::cout << "Issue read -- " << bytes_left_ << " bytes left\n";
	boost::asio::async_read(socket_, request_,
				boost::asio::transfer_at_least(1),
				boost::bind(&KmerRequest::handle_matrix, this,
					    boost::asio::placeholders::error,
					    boost::asio::placeholders::bytes_transferred));
    }
    else
    {
	handle_matrix_sequence(cur_sequence_);

	/*
	 * At the end; write the finis markers.
	 */

	boost::shared_ptr<std::istream> sp(new std::istringstream(">FLUSH\n"));
	boost::shared_ptr<std::istream> null;
	stream_queue_.push_back(sp);
	stream_queue_.push_back(null);

	/*
	 * Now we can start processing.
	 */
	klookup3_->check_queue();

	
//	socket_.close();
//	delete this;
    }
}

void KmerRequest::handle_matrix_sequence(const std::string &seq)
{
    // std::cout << "handle sequence " << seq << "\n";
//    std::cout << seq;

    // Find identifier
    size_t i = seq.find_first_of(" \t\n", 1);
    std::string ident = seq.substr(1, i-1);
    // std::cout << "found '" << ident << "'\n";
    
    boost::shared_ptr<std::istream> sp(new std::istringstream(seq));
    stream_queue_.push_back(sp);
    // std::cout << "queue size " << stream_queue_.size() << "\n";
    // Don't start compute until we have built our entire list of wanted pegs.
    // klookup3_->check_queue();
}

void KmerRequest::handle_post_body(boost::system::error_code err, size_t bytes)
{
    if (!err)
    {
	auto kv = headers_.find("content-length");
	if (kv != headers_.end())
	{
	    size_t content_size = std::stoul(kv->second);
	    // std::cout << "reading " << content_size << "\n";
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

    std::string opts("");
    auto it = headers_.find("kmer-options");
    if (it != headers_.end())
    {
	opts = it->second;
    }

    if (path_ == "/lookup")
    {
	response_stream_ << "HTTP/1.1 200 OK\n";
	response_stream_ << "Content-type: text/plain\n";
	response_stream_ << "\n";
	krequest_ = new std::istream(&request_);
	klookup_ = new KmerLookupClient(kguts_, response_stream_, *krequest_, mapping_,
					boost::bind(&KmerRequest::request_complete, this),
					parameters_);
    }
    else if (path_ == "/matrix2")
    {
	krequest_ = new std::istream(&request_);

	klookup2_ = new KmerLookupClient2(kguts_, *krequest_,
					  boost::bind(&KmerRequest::on_matrix_protein, this, _1, _2),
					  boost::bind(&KmerRequest::on_matrix_hit, this, _1),
					  0,
					  boost::bind(&KmerRequest::matrix_complete, this, _1),
					  parameters_);
    }
    else if (path_ == "/add")
    {
	response_stream_ << "HTTP/1.1 200 OK\n";
	response_stream_ << "Content-type: text/plain\n";
	response_stream_ << "\n";

	krequest_ = new std::istream(&request_);

	klookup2_ = new KmerLookupClient2(kguts_, *krequest_,
					  boost::bind(&KmerRequest::on_protein, this, _1, _2),
					  boost::bind(&KmerRequest::on_hit, this, _1),
					  boost::bind(&KmerRequest::on_call, this, _1),
					  boost::bind(&KmerRequest::add_complete, this, _1),
					  parameters_);
    }
    else if (path_ == "/query")
    {
	response_stream_ << "HTTP/1.1 200 OK\n";
	response_stream_ << "Content-type: text/plain\n";
	response_stream_ << "\n";

	krequest_ = new std::istream(&request_);
	klookup2_ = new KmerLookupClient2(kguts_, *krequest_,
					  0,
					  0,
					  boost::bind(&KmerRequest::on_call, this, _1),
					  boost::bind(&KmerRequest::add_complete, this, _1),
					  parameters_);
					
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

void KmerRequest::on_protein(const std::string &protein, size_t len)
{
    cur_protein_ = protein;
    cur_protein_id_ = mapping_.encode_id(protein);
    // std::cout << "on protein " << protein << "\n";
}

void KmerRequest::on_call(const std::string &line)
{
    // std::cout << "call " << function << "\n";
    response_stream_ << line << "\n";
}

void KmerRequest::on_hit(unsigned long kmer)
{
    mapping_.add_mapping(cur_protein_id_, kmer);
}

void KmerRequest::add_complete( const boost::system::error_code& err )
{
//     std::cout << "Got add complete\n";
    
    boost::asio::async_write(socket_, response_,
			     boost::bind(&KmerRequest::write_response_complete, this,
					 boost::asio::placeholders::error));

    delete krequest_;
    delete klookup2_;
    krequest_ = 0;
    klookup2_ = 0;
}

void KmerRequest::on_matrix_protein(const std::string &protein, size_t len)
{
    cur_protein_ = protein;
    cur_protein_id_ = mapping_.encode_id(protein);
    cur_protein_len_ = len;
    // std::cout << "on protein " << protein << " len=" << len << "\n";
    matrix_proteins_[cur_protein_id_] = len;
}

void KmerRequest::on_matrix_hit(unsigned long kmer)
{
//    std::cout << "got hit " << kmer << "\n";

    auto ki = mapping_.kmer_to_id_.find(kmer);
    if (ki != mapping_.kmer_to_id_.end())
  {
//	std::cout << "got mapping for " << kmer << "\n";
	for (auto it = ki->second.begin(); it != ki->second.end(); it++)
	{
//	    std::cout << "  " << *it << " " << mapping_.decode_id(*it) << "\n";

	    if (*it != cur_protein_id_ && matrix_proteins_.find(*it) != matrix_proteins_.end())
	    {
//		std::cout << "Add " << cur_protein_ << " " << mapping_.decode_id(*it) << "\n";
		distance_[std::make_pair(cur_protein_id_, *it)]++;
	    }
	    else
	    {
//		std::cout << "drop " << cur_protein_ << " " << mapping_.decode_id(*it) << "\n";
	    }
	}
    }

}



void KmerRequest::matrix_complete( const boost::system::error_code& err )
{
    response_stream_ << "HTTP/1.1 200 OK\n";
    response_stream_ << "Content-type: text/plain\n";
    response_stream_ << "\n";
    
//    std::cout << "Got matrix complete\n";

    for (auto it = distance_.begin(); it != distance_.end(); it++)
    {
	std::string p1 = mapping_.decode_id(it->first.first);
	std::string p2 = mapping_.decode_id(it->first.second);
	size_t l1 = matrix_proteins_[it->first.first];
	size_t l2 = matrix_proteins_[it->first.second];
	float score = (float) it->second / ((float) (l1 + l2));
	// std::cout << p1 << "\t" << p2 << "\t" << it->second << "\n";
	response_stream_ << p1 << "\t" << p2 << "\t" << it->second << "\t" << score << "\n";
    }
    boost::asio::async_write(socket_, response_,
			     boost::bind(&KmerRequest::write_response_complete, this,
					 boost::asio::placeholders::error));
}


void KmerRequest::request_complete()
{
    boost::asio::async_write(socket_, response_,
			     boost::bind(&KmerRequest::write_response_complete, this,
					 boost::asio::placeholders::error));

//    std::cout << "request complete, deleting klookup\n";
    delete krequest_;
    delete klookup_;
    krequest_ = 0;
    klookup_ = 0;
}

void KmerRequest::write_response_complete(boost::system::error_code err)
{
//    std::cout << "write_response_complete close()\n";
    socket_.close();

    delete this;
}
