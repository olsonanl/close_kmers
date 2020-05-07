#include "fq_process_request.h"
#include "kserver.h"

#include <string>
#include <memory>

#include <boost/bind.hpp>
#include <boost/core/demangle.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/buffers_iterator.hpp>
#include "global.h"
#include "family_mapper.h"

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

FqProcessRequest::FqProcessRequest(std::shared_ptr<KmerRequest2> owner, std::shared_ptr<KmerPegMapping> mapping,
			     bool family_mode, size_t content_length, bool chunked) :
    family_mode_(family_mode),
    mapping_(mapping),
    content_length_(content_length),
    chunked_(chunked),
    owner_(owner),
    trans_table_(TranslationTable::make_table(11)),
    header_written_(false)
{
}

void FqProcessRequest::run()
{
    if (content_length_ == 0)
    {
	owner_->respond(200, "OK", "data done\n", [](){});
	return;
    }
    if (owner_->request().size() > 0)
    {
	boost::system::error_code err;
	on_first_data(err, owner_->request().size());
    }
    else
    {
	boost::asio::async_read(owner_->socket(), owner_->request(),
				boost::asio::transfer_at_least(content_length_),
				boost::bind(&FqProcessRequest
					    ::on_first_data, this,
					    boost::asio::placeholders::error,
					    boost::asio::placeholders::bytes_transferred));
    }
}

void FqProcessRequest::on_first_data(boost::system::error_code err, size_t bytes)
{
    if (!err || err == boost::asio::error::eof)
    {
	//
	// Determine if this is a gzip format request
	//
	zlib_support::gzip_header hdr;
 	std::string r = make_string(owner_->request());
	// std::cerr << "bytes=" << bytes << " content_length_=" << content_length_ << " err=" << err <<  "\n";
	//std::cerr << "Read buffer contains: "  << r << std::endl;
	
	boost::asio::streambuf::const_buffers_type bufs = owner_->request().data();
	content_length_ -= bytes;

	auto x = boost::asio::buffers_begin(bufs);
	std::copy_n(x, sizeof(hdr), (char *) &hdr);
	// std::cerr << "hdr: " << hdr.sig1 << " " << hdr.sig2 << std::endl;
	if (hdr.sig1 == 0x1f && hdr.sig2 == 0x8b)
	{
	    // std::cerr << "is gzip. compression=" << int(hdr.cm) << " lastmod=" << hdr.timestamp << " os=" << hdr.os << std::endl;

	    //
	    // Here, we will launch the streambuf over to a thread in the thread
	    // pool to manage the decompression.
	    //
	    gzip_decoder_ = std::make_shared<zlib_support::GzipDecoder>();

	    owner_->thread_pool()->post(std::bind(&FqProcessRequest::process_compressed_block, this));
	}
	else
	{
	    owner_->thread_pool()->post(std::bind(&FqProcessRequest::process_block, this, err == boost::asio::error::eof));
	}
    }
    else
    {
	std::cerr << "ERROR in fq on_first_data: " << err << "\n";
    }
	
}

/*
 * Process the block of compressed data in owner_->request().
 * Here, we are executing in a worker thread. We will decompress the
 * block, parse the resulting fastq, analyze, and push results before
 * scheduling the next read from input.
 */
void FqProcessRequest::process_compressed_block()
{
    gzip_decoder_->consume(owner_->request(),
			   std::bind(&FqProcessRequest::process_decompressed_block, this,
				     std::placeholders::_1, std::placeholders::_2));
}			

void FqProcessRequest::process_decompressed_block(std::shared_ptr<boost::asio::streambuf> buf, zlib_support::ReturnStatus status)
{
    auto more_data = [this]() {
	// std::cerr << "request transfer pdb " << content_length_ << std::endl;
	boost::asio::async_read(owner_->socket(), owner_->request(),
				boost::asio::transfer_at_least(content_length_),
				boost::bind(&FqProcessRequest::on_compressed_data, this,
					    boost::asio::placeholders::error,
					    boost::asio::placeholders::bytes_transferred));
    };

    if (buf)
    {
	// std::cerr << "got block decompressed status=" << int(status) << " buf is set\n";

	process_data(*buf, status == zlib_support::ReturnStatus::StreamFinished, [this, status, more_data]() {
		//std::cerr << "IN process_decompressed_block data callback status=" << int(status) << std::endl;
		/*
		 * If we need more input, dispatch the async transfer on our IO context.
		 */
		if (status == zlib_support::ReturnStatus::InputEmpty)
		{
		    more_data();
		}
	    });
    }
    else if (status == zlib_support::ReturnStatus::InputEmpty)
    {
	//std::cerr << "got block decompressed status=" << int(status) << " buf is zero\n";
	
	/* Request more data */
	owner_->io_service().post([this]() {
		//std::cerr << "request transfer pdb2 " << content_length_ << std::endl;
		boost::asio::async_read(owner_->socket(), owner_->request(),
				boost::asio::transfer_at_least(content_length_),
				boost::bind(&FqProcessRequest::on_compressed_data, this,
					    boost::asio::placeholders::error,
					    boost::asio::placeholders::bytes_transferred));
	    });
    }
    else
    {
	std::cerr << "Close socket " << __LINE__ << std::endl;
	std::cerr << "ODD STATUS buf=0, status=" << int(status) << std::endl;
	owner_->socket().close();
	owner_->exit_request();
    }
}

/*
 * process a block of data that we received during a non-compressed request
 * The data is in the owner's request_ streambuf.
 */
void FqProcessRequest::process_block(bool finished)
{
    // std::cerr << "process_block finished=" << finished << "\n";

    process_data(owner_->request(), finished,
		 [this]() {
		     boost::asio::async_read(owner_->socket(), owner_->request(),
					     boost::asio::transfer_at_least(content_length_),
					     boost::bind(&FqProcessRequest::on_uncompressed_data, this,
							 boost::asio::placeholders::error,
							 boost::asio::placeholders::bytes_transferred));
		 });
}

void FqProcessRequest::on_uncompressed_data(boost::system::error_code err, size_t bytes)
{
    // std::cerr << "on_uncompressed_data: err=" << err << " size=" << bytes << std::endl;
    if ( err == boost::asio::error::eof)
    {
	// std::cerr << "on_compressed_data: EOF\n";
	owner_->thread_pool()->post(std::bind(&FqProcessRequest::process_block, this, true));
    }
    else if (!err)
    {
	content_length_ -= bytes;
	owner_->thread_pool()->post(std::bind(&FqProcessRequest::process_block, this, bytes == 0));
    }
    else
    {
	std::cerr << "ERROR in fq on_compressed_data: " << err << "\n";
	owner_->socket().close();
	owner_->exit_request();
    }
	
}


void FqProcessRequest::on_compressed_data(boost::system::error_code err, size_t bytes)
{
    // std::cerr << "on_compressed_data: err=" << err << " size=" << bytes << std::endl;
    if ( err == boost::asio::error::eof)
    {
	// std::cerr << "on_compressed_data: EOF\n";
	owner_->socket().close();
	owner_->exit_request();
    }
    else if (!err)
    {
	content_length_ -= bytes;
	owner_->thread_pool()->post(std::bind(&FqProcessRequest::process_compressed_block, this));
    }
    else
    {
	std::cerr << "ERROR in fq on_compressed_data: " << err << "\n";
	owner_->socket().close();
	owner_->exit_request();
    }
	
}

void FqProcessRequest::process_data(boost::asio::streambuf &buf, bool finished,
				    std::function<void()> cb)

{
    size_t sz = std::distance(buffers_begin(buf.data()), buffers_end(buf.data()));
    // std::cerr << "sz=" << sz << std::endl;
    //auto s = std::string(buffers_begin(buf.data()), buffers_begin(buf.data()) + 10);
    //std::cerr << s << std::endl;

    KmerGuts *kguts = owner_->thread_pool()->kguts_.get();
    FamilyMapper mapper(kguts, mapping_);
    
    auto outbuf = std::make_shared<boost::asio::streambuf>();
    std::ostream os(outbuf.get());

    if (!header_written_)
    {
	owner_->write_header(os, 200, "OK");
	os << "\n";
	header_written_ = true;
    }

    if (sz > 0)
    {
	parser_.set_callback(std::bind(&FqProcessRequest::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2, std::ref(mapper), std::ref(os)));
	
	for (auto x = boost::asio::buffers_begin(buf.data()); x != boost::asio::buffers_end(buf.data()); x++)
	{
	    parser_.parse_char(*x);
	}
    }

    if (finished)
    {
	parser_.parse_complete();
    }

    buf.consume(sz);
    parser_.set_callback(0);

    /*
     * We have processed this block, with output in outbuf. Start a write on the owner I/O thread,
     * and invoke our completion callback when done (or close out the request if we were finished).
     *
     */     

    owner_->io_service().post([this, outbuf, finished, cb]() {
	    boost::asio::async_write(owner_->socket(), boost::asio::buffer(outbuf->data()),
				     [this, finished, cb, outbuf](const boost::system::error_code &err2, const long unsigned int &bytes2) {
					 /* completion of async write. If finished, close req. otherwise,
					  * invoke completion callback to queue appropriate read.
					  */
					 if (finished)
					 {
					     // std::cerr << "Close socket " << __LINE__ << std::endl;
					     owner_->socket().close();
					     owner_->exit_request();
					 }
					 else
					 {
					     cb();
					 }
				     });
	});
}

int FqProcessRequest::on_parsed_seq(const std::string &id, const std::string &seq, FamilyMapper &mapper,
				    std::ostream &os)
{
    if (id.empty())
	return 0;

    DNASequence dna(id, seq);

    auto prots = dna.get_possible_proteins(trans_table_);

    /*
     * We scan all 6 frames returned from get_possible_proteins.
     * We run the analysis on the proteins from each frame, saving the results
     * in a list.
     * If a newer result has a better aggregate score, replace the current best.
     *
     * bestp is a pointer to the frame data with the best score
     * best_score is the best score.
     * best_results is a list of the results from that frame.
     */
    auto bestp = prots.end();
    double best_score = 0.0;
    int best_frame = 0;
    std::vector<std::pair<size_t, FamilyMapper::best_match_t>> best_matches;
	    
    for (auto iter = prots.begin(); iter != prots.end(); iter++)
    {
	int frame = iter->first;
	std::list<std::string> &proteins = iter->second;
	double score = 0.0;
	std::vector<std::pair<size_t, FamilyMapper::best_match_t>> matches;
	int i = 0;
	for (auto prot: proteins)
	{
	    i++;
	    if (prot.length() > 10)
	    {
		matches.emplace_back(std::make_pair(prot.length(), mapper.find_best_family_match(id, prot)));
		auto &match = matches.back();
		score += match.second.score;
	    }
	    // std::cerr << "id =" << id << "i=" << i << " score=" << score << "\n";
	    if (score > best_score)
	    {
		best_score = score;
		bestp = iter;
		best_frame = frame;
		best_matches = matches;
	    }
	}
    }
    if (best_score > 0.0)
    {
	os << id << "\t" << best_frame << "\t" << best_score << "\t";
	auto it = best_matches.begin();
	os << it->first << "\t" << it->second;
	++it;
	while (it != best_matches.end())
	{
	    os << "\t" << it->first << "\t" << it->second;
	    ++it;
	}

	os << std::endl;
    }
    // std::cerr << "p: " << id << " " << seq << std::endl;
    return 0;
}
