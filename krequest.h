#ifndef _KREQUEST_H
#define _KREQUEST_H

/*
 * A kmer request wraps up the code for parsing out
 * the incoming HTTP request from a client and handling
 * the execution of the request itself.
 */

#include <vector>
#include <string>
#include <set>
#include <boost/asio.hpp>
#include <boost/timer/timer.hpp>
#include "kmer.h"
#include "klookup.h"
#include "klookup2.h"
#include "klookup3.h"

class KmerRequest
{
public:
    KmerRequest(boost::asio::io_service &io_service,
		KmerPegMapping &mapping,
		boost::asio::ip::tcp::endpoint &klookup_endpoint);
    ~KmerRequest();

    void do_read();

    boost::asio::ip::tcp::socket& socket()
    {
	return socket_;
    }

private:

    void handle_read(boost::system::error_code err, size_t bytes);
    void handle_post_body(boost::system::error_code err, size_t bytes);
    void handle_request();

    void process_request();
    void request_complete( const KmerLookupClient::result_t &);

    void write_response_complete(boost::system::error_code err);
    
    void on_protein(const std::string &protein);
    void on_hit(unsigned long kmer);
    void on_call(const std::string &line);
    void add_complete( const boost::system::error_code& err );

    void handle_matrix(boost::system::error_code err, size_t bytes);
    void on_matrix_protein(const std::string &protein);
    void on_matrix_hit(unsigned long kmer);
    void matrix_complete( const boost::system::error_code& err );

    std::string request_type_;
    std::string path_;
    std::map<std::string, std::string> headers_;

    boost::asio::io_service &io_service_;
    boost::asio::ip::tcp::socket socket_;
    boost::asio::streambuf request_;
    KmerPegMapping &mapping_;
    boost::asio::ip::tcp::endpoint klookup_endpoint_;
    KmerLookupClient *klookup_;
    KmerLookupClient2 *klookup2_;
    std::istream *krequest_;

    KmerLookupClient3::stream_queue_t stream_queue_;
    boost::shared_ptr<KmerLookupClient3> klookup3_;

    //
    // Distance matrix support.
    //
    void handle_matrix_sequence(const std::string &seq);
    std::string cur_sequence_;
    size_t bytes_left_;
    std::set<KmerPegMapping::encoded_id_t> matrix_proteins_;
    std::map<std::pair<KmerPegMapping::encoded_id_t, KmerPegMapping::encoded_id_t>, unsigned long> distance_;

    boost::asio::streambuf response_;
    //framed_streambuf response_;
    std::ostream response_stream_;
    boost::timer::cpu_timer timer_;

    std::string cur_protein_;
    KmerPegMapping::encoded_id_t cur_protein_id_;
};


#endif
