#ifndef _STREAMUTIL_H
#define _STREAMUTIL_H


#include <functional>
#include <iostream>
#include <string>

#include <boost/asio.hpp>

/// @brief Type that derives from Boost.Asio streambuf and can frame the
///        input sequence to a portion of the actual input sequence.
template <typename Allocator = std::allocator<char> >
    class basic_framed_streambuf
    : public boost::asio::basic_streambuf<Allocator>
{
private:

typedef boost::asio::basic_streambuf<Allocator> parent_type;

public:

  explicit
basic_framed_streambuf(
    std::size_t maximum_size = (std::numeric_limits< std::size_t >::max)(),
    const Allocator& allocator = Allocator()
    )
: parent_type(maximum_size, allocator),
egptr_(nullptr)
{}

/// @brief Limit the current input sequence to n characters.
///
/// @remark An active frame is invalidated by any member function that
///        modifies the input or output sequence.
void frame(std::streamsize n)
{
    // Store actual end of input sequence.
    egptr_ = this->egptr();
    // Set the input sequence end to n characters from the current
    // input sequence pointer..
    this->setg(this->eback(), this->gptr(), this->gptr() + n);
}

/// @brief Restore the end of the input sequence.
void unframe()
{
    // Restore the end of the input sequence.
    this->setg(this->eback(), this->gptr(), this->egptr_);
    egptr_ = nullptr;
}

protected:

// When the end of the input sequence has been reached, underflow
// will be invoked.
typename parent_type::int_type underflow()
{
    // If the  streambuf is currently framed, then return eof
    // on underflow.  Otherwise, defer to the parent implementation.
    return egptr_ ? parent_type::traits_type::eof()
    : parent_type::underflow();
}

private:
char* egptr_;
};

typedef basic_framed_streambuf<> framed_streambuf;

/// @brief RAII type that helps frame a basic_framed_streambuf within a
///        given scope.
template <typename Streambuf>
class streambuf_frame
{
public:
    explicit streambuf_frame(Streambuf& streambuf, std::streamsize n)
	: streambuf_(streambuf)
    {
	streambuf_.frame(n);
    }

    ~streambuf_frame() { streambuf_.unframe(); }

    streambuf_frame(const streambuf_frame&) = delete;
    streambuf_frame& operator=(const streambuf_frame&) = delete;

private:
    Streambuf& streambuf_;
};

#endif
