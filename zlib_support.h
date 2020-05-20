#ifndef _zlib_support_h
#define _zlib_support_h

#include <zlib.h>
#include <cstdint>
#include <boost/asio/streambuf.hpp>
#include <functional>

namespace zlib_support {

    struct gzip_header
    {
	uint8_t sig1;
	uint8_t sig2;
	uint8_t cm;
	uint8_t flags;
	uint32_t timestamp;
	uint8_t extra_flags;
	uint8_t os;
    } __attribute__ ((packed));

    const int OutputSize = 1048576 * 8;

    enum class ReturnStatus {
	OK,
	    InputEmpty,
	    StreamFinished,
	    Error
	    };
	    

    class GzipDecoder
    {
    public:
	GzipDecoder();

	void consume(boost::asio::streambuf &sb, std::function<void(std::shared_ptr<boost::asio::streambuf> buf, ReturnStatus status)>);

    private:
	z_stream zstream_;
	std::shared_ptr<boost::asio::streambuf> output_buffer_;
    };
}


#endif
