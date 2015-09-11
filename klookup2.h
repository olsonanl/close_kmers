#ifndef _KLOOKUP2_H
#define _KLOOKUP2_H

#include <boost/asio.hpp>
#include <boost/function.hpp>
#include <boost/timer/timer.hpp>
#include <iostream>
#include <string>
#include "kmer.h"
#include "kguts.h"
#include "fasta_parser.h"

using namespace boost::asio::ip;

class KmerLookupClient2
{
public:
    typedef std::vector<std::pair<std::string, unsigned int>> result_t;

    KmerLookupClient2(std::shared_ptr<KmerGuts> kguts,
		      std::istream &input,
		      boost::function<void ( const std::string &prot, size_t len )> on_protein,
		      boost::function<void ( unsigned long kmer )> on_hit,
		      boost::function<void ( const std::string &line )> on_call,
		      boost::function<void ( const boost::system::error_code& err)> on_completion);

private:
    int on_parsed_seq(const std::string &id, const std::string &seq);
    void on_hit(KmerGuts::sig_kmer_t &hit);

    FastaParser fasta_parser_;
    std::shared_ptr<KmerGuts> kguts_;
    boost::asio::streambuf request_;
    boost::asio::streambuf response_;
    char buffer_[1048576];
    std::istream &input_;
    boost::function<void ( const std::string &prot, size_t len )> on_protein_;
    boost::function<void ( unsigned long kmer )> on_hit_;
    boost::function<void ( const std::string &line )> on_call_;
    boost::function<void ( const boost::system::error_code& err )> on_completion_;
    boost::timer::cpu_timer timer_;
};


#endif
