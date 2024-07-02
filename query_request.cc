#include "query_request.h"
#include <boost/bind/bind.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/buffers_iterator.hpp>

inline std::string make_string(boost::asio::streambuf& streambuf)
{
    return {buffers_begin(streambuf.data()),
	    buffers_end(streambuf.data())};
}

QueryRequest::QueryRequest(std::shared_ptr<KmerRequest2> owner, size_t content_length, bool chunked) :
    content_length_(content_length),
    chunked_(chunked),
    owner_(owner),
    header_written_(false)
{
    parser_.set_callback(std::bind(&QueryRequest::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2));
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

		int details = 0;
		try {
		    details = std::stoi(owner_->parameters()["details"]);
		} catch (...) {}

		int find_best_call = 0;
		try {
		    find_best_call = std::stoi(owner_->parameters()["find_best_call"]);
		} catch (...) {}


		for (auto x: *cur)
		{
		    auto id = x.first;
		    auto seq = x.second;
		    
		    auto calls = std::make_shared<std::vector<KmerCall> >();
		    auto stats = std::make_shared<KmerOtuStats>();

		    std::shared_ptr<std::vector<KmerGuts::hit_in_sequence_t> > hits = 0;
		    if (details)
		    {
			hits = std::make_shared<std::vector<KmerGuts::hit_in_sequence_t> >();
			kguts->process_aa_seq_hits(id, seq, calls, hits, stats);
		    }
		    else
		    {
			kguts->process_aa_seq(id, seq, calls, 0, stats);
		    }

		    

		    if (find_best_call)
		    {
			FunctionIndex best_call_fi;
			float best_call_score, best_call_weighted_score, best_call_score_offset;
			std::string best_call_function;
			kguts->find_best_call(*calls, best_call_fi, best_call_function, best_call_score,
					      best_call_weighted_score, best_call_score_offset);
			if (!best_call_function.empty())
			{
			    os << id << "\t" << best_call_function << "\t" << best_call_score << "\t" << best_call_weighted_score << "\n";
			}
		    }
		    else
		    {
			os << "PROTEIN-ID\t" << id << "\t" << seq.size() << "\n";
			for (auto c: *calls)
			{
			    os << kguts->format_call(c);
			}
			if (details)
			{
			    for (auto h: *hits)
			    {
				os << kguts->format_hit(h);
			    }
			}
			os << kguts->format_otu_stats(x.first, x.second.size(), *stats);
		    }
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
    return 0;
}
