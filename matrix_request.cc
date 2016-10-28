#include "matrix_request.h"

#include <boost/bind.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/buffers_iterator.hpp>

inline std::string make_string(boost::asio::streambuf& streambuf)
{
    return {buffers_begin(streambuf.data()),
	    buffers_end(streambuf.data())};
}

MatrixRequest::MatrixRequest(std::shared_ptr<KmerRequest2> owner, std::shared_ptr<KmerPegMapping> mapping, int content_length, bool chunked) :
    mapping_(mapping),
    content_length_(content_length),
    chunked_(chunked),
    owner_(owner),
    parser_(std::bind(&MatrixRequest::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2))
{
}

void MatrixRequest::run()
{
    // std::cerr << "mat run, len=" << content_length_ <<  " request_size=" << owner_->request().size() << "\n";

    if (owner_->request().size() > 0)
    {
	boost::system::error_code err;
	on_data(err, 0);
    }
    else
    {
	boost::asio::async_read(owner_->socket(), owner_->request(),
				boost::asio::transfer_at_least(content_length_),
				boost::bind(&MatrixRequest::on_data, this,
					    boost::asio::placeholders::error,
					    boost::asio::placeholders::bytes_transferred));
    }
}

void MatrixRequest::on_data(boost::system::error_code err, size_t bytes)
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

		for (auto work_item: *cur)
		{
		    auto id = work_item.id();
		    auto seq = work_item.seq();
		    
		    KmerPegMapping::encoded_id_t eid = mapping_->encode_id(id);
		    
		    matrix_proteins_[eid] = seq.size();
		    
		    kguts->process_aa_seq(id, seq, 0,
					  std::bind(&MatrixRequest::on_hit, this, eid, std::placeholders::_1),
					  0);
		}
		owner_->io_service().post([this, sbuf, err](){
			/*
			 * Back in the main thread here. We can write our response.
			 */
			// std::cout << "post response in " << pthread_self() << "\n";

			if (err == boost::asio::error::eof || content_length_ == 0)
			{
			    process_results();
			}
			else
			{
			    boost::asio::async_read(owner_->socket(), owner_->request(),
						    boost::asio::transfer_at_least(content_length_),
						    boost::bind(&MatrixRequest::on_data, this,
								boost::asio::placeholders::error,
								boost::asio::placeholders::bytes_transferred));
			}
		    });
	    });
    }
    else
    {
	std::cerr << "ERROR in matrix_request: " << err << "\n";
    }
	
}

int MatrixRequest::on_parsed_seq(const std::string &id, const std::string &seq)
{
    current_work_->push_back(ProteinSequence(id, seq));
    return 0;
}

void MatrixRequest::on_hit(KmerPegMapping::encoded_id_t id, KmerGuts::hit_in_sequence_t kmer)
{
    auto ki = mapping_->kmer_to_id_.find(kmer.hit.which_kmer);
    if (ki != mapping_->kmer_to_id_.end())
    {
        // char kmerstr[10];
	// KmerGuts::decoded_kmer(kmer.which_kmer, kmerstr);
	// std::cout << "got mapping for " << kmerstr << "\n";
	//std::cout << "got mapping for " << kmer.which_kmer << "\n";
	for (auto eid: ki->second)
	{
	    // std::cout << "  " << eid << " " << mapping_->decode_id(eid) << "\n";
	    
	    if (eid != id)
	    {
		if (matrix_proteins_.find(eid) != matrix_proteins_.end())
		{
		    // std::cout << "    add " << eid << " " << mapping_->decode_id(eid) << "\n";
		    distance_[std::make_pair(id, eid)]++;
		}
		else
		{
		    // std::cout << "    drop " << eid << " " << mapping_->decode_id(eid) << "\n";
		}
	    }
	}
    }
    else
    {
	std::cerr << "no mapping for " << kmer.hit.which_kmer << "\n";
    }
}

void MatrixRequest::process_results()
{
    auto sbuf = std::make_shared<boost::asio::streambuf>();
    std::ostream os(sbuf.get());
    
    os << "HTTP/1.1 200 OK\n";
    os << "Content-type: text/plain\n";
    os << "\n";
    
    for (auto it = distance_.begin(); it != distance_.end(); it++)
    {
	KmerPegMapping::encoded_id_t e1 = it->first.first;
	KmerPegMapping::encoded_id_t e2 = it->first.second;
	
	std::string p1 = mapping_->decode_id(e1);
	std::string p2 = mapping_->decode_id(e2);
	size_t l1 = matrix_proteins_[e1];
	size_t l2 = matrix_proteins_[e2];
	float score = (float) it->second / ((float) (l1 + l2));
	// std::cout << p1 << "\t" << p2 << "\t" << it->second << "\n";
	os << p1 << "\t" << p2 << "\t" << it->second << "\t" << score << "\n";
    }
    boost::asio::async_write(owner_->socket(), *sbuf,
			     [this, sbuf](const boost::system::error_code &err, const long unsigned int &bytes){
				 owner_->socket().close();
				 owner_->exit_request();
			     });
}
