#include "lookup_request.h"

#include <boost/bind.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/buffers_iterator.hpp>
#include "global.h"

template<class T>
struct less_second : std::binary_function<T,T,bool>
{
    inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
};  

inline std::string make_string(boost::asio::streambuf& streambuf)
{
    return {buffers_begin(streambuf.data()),
	    buffers_end(streambuf.data())};
}

LookupRequest::LookupRequest(std::shared_ptr<KmerRequest2> owner, std::shared_ptr<KmerPegMapping> mapping, int content_length, bool chunked) :
    owner_(owner),
    mapping_(mapping),
    content_length_(content_length),
    parser_(std::bind(&LookupRequest::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2)),
    chunked_(chunked),
    header_written_(false)
{
    kmer_hit_threshold_ = 3;
    try {
	kmer_hit_threshold_ = std::stoi(owner->parameters()["kmer_hit_threhsold"]);
    } catch (const std::invalid_argument& ia)
    {
    }
}

void LookupRequest::run()
{
    if (owner_->request().size() > 0)
    {
	boost::system::error_code err;
	on_data(err, 0);
    }
    else
    {
	boost::asio::async_read(owner_->socket(), owner_->request(),
				boost::asio::transfer_at_least(content_length_),
				boost::bind(&LookupRequest::on_data, this,
					    boost::asio::placeholders::error,
					    boost::asio::placeholders::bytes_transferred));
    }
}

void LookupRequest::on_data(boost::system::error_code err, size_t bytes)
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

	std::shared_ptr<work_list_t> cur = current_work_;
	current_work_ = 0;
	owner_->thread_pool()->post([this, cur, err]() {
		KmerGuts *kguts = owner_->thread_pool()->kguts_.get();
		kguts->set_parameters(owner_->parameters());

		auto sbuf = std::make_shared<boost::asio::streambuf>();
		{
		    std::ostream os(sbuf.get());

		    if (!header_written_)
		    {
			owner_->write_header(os, 200, "OK");
			os << "\n";
			header_written_ = true;
		    }

		    for (auto work_item: *cur)
		    {
			auto id = work_item.id();
			auto seq = work_item.seq();
		    
			hit_count_.clear();

			// std::cerr << "Lookup " << id << " " << seq << "\n";
			kguts->process_aa_seq(id, seq, 0,
					      std::bind(&LookupRequest::on_hit, this, std::placeholders::_1),
					      0);


			typedef std::pair<KmerPegMapping::encoded_id_t, unsigned int> data_t;

			std::vector<data_t> vec;
			for (auto it: hit_count_)
			{
			    vec.push_back(it);
			}

			std::sort(vec.begin(), vec.end(), less_second<data_t>()); 

			os << id << "\n";
			for (auto it: vec)
			{
			    auto eid = it.first;
			    auto score = it.second;
			
			    if (score < kmer_hit_threshold_)
				break;
			    std::string peg = mapping_->decode_id(eid);
			    os << peg << "\t" << score;
			    // os << eid << "\t" << peg << "\t" << score;

			    auto fhit = mapping_->family_mapping_.find(eid);
			    if (fhit != mapping_->family_mapping_.end())
			    {
				os << "\t" << fhit->second.pgf << "\t" << fhit->second.plf << "\t" << fhit->second.function;
			    }
			    os << "\n";
			}
			os << "//\n";

		    }
		    os.flush();
		}
		owner_->io_service().post([this, sbuf, err](){
			/*
			 * Back in the main thread here. We can write our response.
			 */
			// std::cout << "post response in " << pthread_self() << "\n";

			std::cerr << "write results size " << sbuf->size() << "\n";
			boost::asio::async_write(owner_->socket(), boost::asio::buffer(sbuf->data()),
						 [this, err, sbuf](const boost::system::error_code &err2, const long unsigned int &bytes2){
						     std::cerr << "write done in " << pthread_self() << " err=" << err.message() << " content_length=" << content_length_ <<"\n";
						     std::cerr << "   err2=" << err2.message() << " bytes2=" << bytes2 << "\n";
						     if (err == boost::asio::error::eof || content_length_ == 0)
						     {
							 process_results();
						     }
						     else
						     {
							 boost::asio::async_read(owner_->socket(), owner_->request(),
										 boost::asio::transfer_at_least(content_length_),
										 boost::bind(&LookupRequest::on_data, this,
											     boost::asio::placeholders::error,
											     boost::asio::placeholders::bytes_transferred));
						     }
						 });
		    });
	    });
    }
    else
    {
	std::cerr << "ERROR in lookup_request: " << err << "\n";
    }
	
}

int LookupRequest::on_parsed_seq(const std::string &id, const std::string &seq)
{
    current_work_->push_back(ProteinSequence(id, seq));
}

void LookupRequest::on_hit(KmerGuts::sig_kmer_t &kmer)
{
    auto ki = mapping_->kmer_to_id_.find(kmer.which_kmer);
    if (ki != mapping_->kmer_to_id_.end())
    {
	// std::cout << "got mapping for " << kmer.which_kmer << "\n";
	for (auto eid: ki->second)
	{
	    hit_count_[eid]++;
	}
    }
}

void LookupRequest::process_results()
{
    owner_->socket().close();
    owner_->exit_request();
}
