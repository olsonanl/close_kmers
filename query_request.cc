#include "query_request.h"

#include <boost/bind.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/buffers_iterator.hpp>

inline std::string make_string(boost::asio::streambuf& streambuf)
{
    return {buffers_begin(streambuf.data()),
	    buffers_end(streambuf.data())};
}

QueryRequest::QueryRequest(std::shared_ptr<KmerRequest2> owner, int content_length, bool chunked) :
    owner_(owner),
    parser_(std::bind(&QueryRequest::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2)),
    content_length_(content_length),
    chunked_(chunked),
    header_written_(false)
{
}

void QueryRequest::run()
{
    // std::cerr << "ar run, len=" << content_length_ <<  " request_size=" << owner_->request().size() << "\n";

    if (owner_->request().size() > 0)
    {
	boost::system::error_code err;
	on_data(err, 0);
    }
    else
    {
	boost::asio::async_read(owner_->socket(), owner_->request(),
				boost::asio::transfer_at_least(content_length_),
				boost::bind(&QueryRequest::on_data, this,
					    boost::asio::placeholders::error,
					    boost::asio::placeholders::bytes_transferred));
    }
}

void QueryRequest::on_data(boost::system::error_code err, size_t bytes)
{
    if (!err || err == boost::asio::error::eof)
    {
 	// std::string r = make_string(owner_->request());
	// std::cerr << "bytes=" << bytes << " content_length_=" << content_length_ << " err=" << err <<  "\n";
	// std::cerr << "Read buffer contains: "  << r << std::endl;
	
	current_work_ = std::make_shared<work_list_t>();

	boost::asio::streambuf::const_buffers_type bufs = owner_->request().data();
	int n = 0;
	for (auto x = boost::asio::buffers_begin(bufs); x != boost::asio::buffers_end(bufs); x++)
	{
	    parser_.parse_char(*x);
	    n++;
	}

	content_length_ -= n;

	// std::cerr << "Processed " << n << " remaining " << content_length_ << "\n";
	if (err == boost::asio::error::eof || content_length_ == 0)
	{
	    parser_.parse_complete();
	}

	owner_->request().consume(n);

	/*
	 * We have parsed out this packet.
	 * We post the data to the thread pool for processing, and the
	 * results are to be posted back here so we can push them to the
	 * outgoing socket and queue the next read.
	 */

	// std::cerr << "Process set of " << current_work_->size() << "in " << pthread_self() << "\n";

	std::shared_ptr<work_list_t> cur = current_work_;
	current_work_ = 0;
	owner_->thread_pool()->post([this, cur, err ]() {
		KmerGuts *kguts = owner_->thread_pool()->kguts_.get();
		kguts->set_parameters(owner_->parameters());
		auto sbuf = std::make_shared<boost::asio::streambuf>();
		std::ostream os(sbuf.get());

		if (!header_written_)
		{
		    owner_->write_header(os, 200, "OK");
		    os << "\n";
		    header_written_ = true;
		}

		for (auto x: *cur)
		{
		    auto id = x.first;
		    auto seq = x.second;
		    
		    auto calls = std::make_shared<std::vector<KmerCall> >();
		    auto stats = std::make_shared<KmerOtuStats>();
		    kguts->process_aa_seq(id, seq, calls, 0, stats);
		    os << "PROTEIN-ID\t" << id << "\t" << seq.size() << "\n";
		    for (auto c: *calls)
		    {
			os << kguts->format_call(c);
		    }
		    os << kguts->format_otu_stats(x.first, x.second.size(), *stats);
		}
		owner_->io_service().post([this, sbuf, err](){
			/*
			 * Back in the main thread here. We can write our response.
			 */

			boost::asio::async_write(owner_->socket(), *sbuf,
						 [this, err, sbuf](const boost::system::error_code &err2, const long unsigned int &bytes2){
						     // std::cerr << "write done in " << pthread_self() << " err=" << err << " content_length=" << content_length_ <<"\n";
						     if (err == boost::asio::error::eof || content_length_ == 0)
						     {
							 // std::cerr << "packet complete\n";
							 owner_->socket().close();
							 owner_->exit_request();
						     }
						     else
						     {
							 boost::asio::async_read(owner_->socket(), owner_->request(),
										 boost::asio::transfer_at_least(content_length_),
										 boost::bind(&QueryRequest::on_data, this,
											     boost::asio::placeholders::error,
											     boost::asio::placeholders::bytes_transferred));
						     }
						 });
		    });
	    });
    }
    else
    {
	std::cerr << "ERROR in query_request: " << err << "\n";
    }
	
}

int QueryRequest::on_parsed_seq(const std::string &id, const std::string &seq)
{
    current_work_->push_back(std::make_pair(id, seq));
}
