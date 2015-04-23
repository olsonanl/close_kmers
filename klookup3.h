#ifndef _KLOOKUP3_H
#define _KLOOKUP3_H

#include <boost/asio.hpp>
#include <boost/function.hpp>
#include <boost/timer/timer.hpp>
#include <iostream>
#include <deque>
#include <boost/shared_ptr.hpp>
#include <string>
#include "kmer.h"

using namespace boost::asio::ip;

class KmerLookupClient3
{
public:
    typedef std::vector<std::pair<std::string, unsigned int>> result_t;
    typedef std::deque<boost::shared_ptr<std::istream>> stream_queue_t;

    KmerLookupClient3(boost::asio::io_service& io_service,
		      boost::asio::ip::tcp::endpoint,
		      stream_queue_t &stream_queue,
		      boost::function<void ( const std::string &prot )> on_protein,
		      boost::function<void ( unsigned long kmer )> on_hit,
		      boost::function<void ( const std::string &line )> on_call,
		      boost::function<void ( const boost::system::error_code& err)> on_completion);
    void check_queue();


private:
    void handle_resolve(const boost::system::error_code& err,
			tcp::resolver::iterator endpoint_iterator);
    void handle_connect(const boost::system::error_code& err);

    void handle_write_request(const boost::system::error_code& err);
    void finish_write_request(const boost::system::error_code& err);
    void handle_read(const boost::system::error_code& err, size_t bytes_transferred);

    tcp::resolver resolver_;
    tcp::socket socket_;
    boost::asio::streambuf request_;
    boost::asio::streambuf response_;
    char buffer_[4096];
    boost::function<void ( const std::string &prot )> on_protein_;
    boost::function<void ( unsigned long kmer )> on_hit_;
    boost::function<void ( const std::string &line )> on_call_;
    boost::function<void ( const boost::system::error_code& err )> on_completion_;
    boost::timer::cpu_timer timer_;
    stream_queue_t &stream_queue_;
    bool write_pending_;
};


#endif
