#ifndef _KLOOKUP_H
#define _KLOOKUP_H

#include <boost/asio.hpp>
#include <boost/function.hpp>
#include <boost/timer/timer.hpp>
#include <iostream>
#include <string>
#include "kmer.h"

using namespace boost::asio::ip;

class KmerLookupClient
{
public:
    typedef std::vector<std::pair<std::string, unsigned int>> result_t;

    KmerLookupClient(boost::asio::io_service& io_service,
		     boost::asio::ip::tcp::endpoint,
		     std::istream &input,
		     KmerPegMapping &mapping,
		     boost::function<void ( const result_t & )> on_completion);

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
    KmerPegMapping &mapping_;
    std::map<KmerPegMapping::encoded_id_t, unsigned int> hit_count_;

    boost::function<void ( const result_t & )> on_completion_;
    boost::timer::cpu_timer timer_;
};


#endif
