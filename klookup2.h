#ifndef _KLOOKUP2_H
#define _KLOOKUP2_H

#include <boost/asio.hpp>
#include <boost/function.hpp>
#include <boost/timer/timer.hpp>
#include <iostream>
#include <string>
#include "kmer.h"

using namespace boost::asio::ip;

class KmerLookupClient2
{
public:
    typedef std::vector<std::pair<std::string, unsigned int>> result_t;

    KmerLookupClient2(boost::asio::io_service& io_service,
		      boost::asio::ip::tcp::endpoint,
		      std::istream &input,
		      boost::function<void ( const std::string &prot )> on_protein,
		      boost::function<void ( unsigned long kmer )> on_hit,
		      boost::function<void ( const std::string &function, const std::string &count )> on_call,
		      boost::function<void ( const boost::system::error_code& err)> on_completion);

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
    std::istream &input_;
    boost::function<void ( const std::string &prot )> on_protein_;
    boost::function<void ( unsigned long kmer )> on_hit_;
    boost::function<void ( const std::string &function, const std::string &count )> on_call_;
    boost::function<void ( const boost::system::error_code& err )> on_completion_;
    boost::timer::cpu_timer timer_;
};


#endif
