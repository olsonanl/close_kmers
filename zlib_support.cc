#include "zlib_support.h"

#include <string.h>
#include <boost/asio/buffers_iterator.hpp>

#include <iostream>

using namespace zlib_support;

inline std::string make_string(boost::asio::streambuf& streambuf)
{
    return {buffers_begin(streambuf.data()),
	    buffers_end(streambuf.data())};
}

GzipDecoder::GzipDecoder()
    : output_buffer_(std::make_shared<boost::asio::streambuf>(OutputSize))
{
    memset(&zstream_, 0, sizeof(zstream_));
    // initialize zlib for gzip
    auto ret = inflateInit2(&zstream_, 16 + MAX_WBITS);
    
    if (ret != Z_OK)
	throw std::runtime_error("zlib initializatiion failed");

}

void GzipDecoder::consume(boost::asio::streambuf &sb, std::function<void(std::shared_ptr<boost::asio::streambuf> buf, ReturnStatus status)> cb)
{
    /*
     * Set up input from the incoming streambuf
     */
    
    boost::asio::streambuf::const_buffers_type bufs = sb.data();
    auto beg = boost::asio::buffers_begin(bufs);
    auto end = boost::asio::buffers_end(bufs);
    size_t input_size = zstream_.avail_in = std::distance(beg, end);
    #if BOOST_VERSION <= 105900
    zstream_.next_in = const_cast<unsigned char *>(boost::asio::buffer_cast<const unsigned char*>(bufs));
    #else
    zstream_.next_in = (unsigned char *) bufs.data();
    #endif

    // std::cerr << "dist=" << zstream_.avail_in << std::endl;

    int ret = -1;

    /*
     * And loop until we've generated all the output we can
     */
    do {
	/*
	 * Point our initial output buffer at the output streambuf
	 */
	
	boost::asio::streambuf::mutable_buffers_type obufs = output_buffer_->prepare(OutputSize);

	#if BOOST_VERSION <= 105900
	zstream_.next_out = boost::asio::buffer_cast<unsigned char*>(obufs);
	#else
	zstream_.next_out = (unsigned char *) obufs.data();
	#endif

	zstream_.avail_out = OutputSize;
	
	ret = inflate(&zstream_, Z_NO_FLUSH);
	sb.consume(input_size);
	
	// std::cerr << "ret=" << ret << " size " << (OutputSize - zstream_.avail_out) << std::endl;

	/*
	 * Check output.
	 * We will need an error callback to terminate the request in the main thread.
	 */
	if (ret != Z_OK && ret != Z_STREAM_END)
	{
	    std::cerr << "error " << ret << " from inflate\n";
	    inflateEnd(&zstream_);
	    cb(0, ReturnStatus::Error);
	    return;
	}
	    
	output_buffer_->commit(OutputSize - zstream_.avail_out);

	cb(output_buffer_, ret == Z_OK ? ReturnStatus::OK : ReturnStatus::StreamFinished);

	/*
	 * And make a new output buffer.
	 */
	output_buffer_ = std::make_shared<boost::asio::streambuf>(OutputSize);
	
    } while (zstream_.avail_out == 0);
    if (ret == Z_OK)
	cb(0, ReturnStatus::InputEmpty);
    else
	inflateEnd(&zstream_);
}
