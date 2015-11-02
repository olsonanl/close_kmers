#include "add_request.h"

#include <boost/bind.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/buffers_iterator.hpp>
#include <array>
#include <cctype>
#include "md5.h"

inline std::string make_string(boost::asio::streambuf& streambuf)
{
    return {buffers_begin(streambuf.data()),
	    buffers_end(streambuf.data())};
}

AddRequest::AddRequest(std::shared_ptr<KmerRequest2> owner, std::shared_ptr<KmerPegMapping> mapping, int content_length, bool chunked) :
    owner_(owner),
    mapping_(mapping),
    parser_(std::bind(&AddRequest::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2)),
    content_length_(content_length),
    chunked_(chunked), silent_(0)
{
    silent_ = 0;
    try {
	silent_ = std::stoi(owner_->parameters()["silent"]);
    } catch (const std::invalid_argument& ia)
    {
    }


}

void AddRequest::run()
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
				boost::bind(&AddRequest::on_data, this,
					    boost::asio::placeholders::error,
					    boost::asio::placeholders::bytes_transferred));
    }
}

void AddRequest::on_data(boost::system::error_code err, size_t bytes)
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
	/*
	for (auto x: *current_work_)
	{
	    std::cerr << x.first << " ";
	}
	std::cerr << "\n";
	*/

	typedef std::shared_ptr<std::vector<KmerGuts::hit_in_sequence_t> > hlist_t;
	std::shared_ptr<work_list_t> cur = current_work_;
	//std::vector<std::pair<std::string, hlist_t> > seq_hits;
	auto seq_hits = std::make_shared<std::vector<std::pair<std::string, hlist_t> > >();
	current_work_ = 0;
	owner_->thread_pool()->post([this, cur, err, seq_hits ]() {
		KmerGuts *kguts = owner_->thread_pool()->kguts_.get();
		kguts->set_parameters(owner_->parameters());
		auto sbuf = std::make_shared<boost::asio::streambuf>();
		// auto md5s = std::make_shared<std::map<std::string, std::array<unsigned char, 16> > >();
		std::ostream os(sbuf.get());

		// std::cout << "compute in " << pthread_self() << "\n";

#ifdef USE_TBB
		KmerPegMapping &m = *mapping_;
#endif
		try {
		    for (auto x: *cur)
		    {
			std::string &id = x.first;
			std::string &seq = x.second;
			/*
			  md5_state_t mstate;
			  md5_init(&mstate);
			  for (char c: seq)
			  {
			  md5_byte_t u = toupper(c);
			  md5_append(&mstate, &u, 1);
			  }
			  std::array<unsigned char, 16> &c = (*md5s)[id];
			  md5_finish(&mstate, c.data());
			*/
			hlist_t hits = std::make_shared<std::vector<KmerGuts::hit_in_sequence_t> >();
			auto calls = std::make_shared<std::vector<KmerCall> >();
			auto stats = std::make_shared<KmerOtuStats>();
			kguts->process_aa_seq_hits(id, seq, calls, hits, stats);
			if (!silent_)
			{
			    os << "PROTEIN-ID\t" << id << "\t" << seq.size() << "\n";
			    for (auto c: *calls)
			    {
				os << kguts->format_call(c);
			    }
			    os << kguts->format_otu_stats(id, seq.size(), *stats);
			}
#ifdef USE_TBB
			KmerPegMapping::encoded_id_t enc_id = m.encode_id(id);
			std::vector<KmerGuts::hit_in_sequence_t> &hlist = *hits;
			for (auto hit: hlist)
			{
			    m.add_mapping(enc_id, hit.hit.which_kmer);
			}
#else
			seq_hits->push_back(std::make_pair(id, hits));
#endif
		    }
		}
		catch (std::exception &e)
		{
		    std::cerr << "ending add_request due to exception " << e.what() << "\n";
		    owner_->socket().close();
		    owner_->exit_request();
		    return;
		}
		catch (...)
		{
		    std::cerr << "ending add_request due to default exception\n";
		    owner_->socket().close();
		    owner_->exit_request();
		    return;
		}
		
		owner_->io_service().post([this, sbuf, seq_hits, err](){
			/*
			 * Back in the main thread here. We can write our response.
			 */
			// std::cout << "post response in " << pthread_self() << "\n";

			#ifndef USE_TBB
			for (auto result: *seq_hits)
			{
			    const std::string &id = result.first;
			    KmerPegMapping::encoded_id_t enc_id = mapping_->encode_id(id);
			    for (auto hit: *result.second)
			    {
				mapping_->add_mapping(enc_id, hit.which_kmer);
			    }
			}
			#endif
			

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
										 boost::bind(&AddRequest::on_data, this,
											     boost::asio::placeholders::error,
											     boost::asio::placeholders::bytes_transferred));
						     }
						 });
		    });
	    });
    }
    else
    {
	std::cerr << "ERROR in add_request: " << err << "\n";
    }
	
}

int AddRequest::on_parsed_seq(const std::string &id, const std::string &seq)
{
    current_work_->push_back(std::make_pair(id, seq));
    // std::cerr << "seq: " << id << "\n";
    // std::cerr << seq << "\n";
}
