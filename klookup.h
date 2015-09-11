#ifndef _KLOOKUP_H
#define _KLOOKUP_H

#include <boost/asio.hpp>
#include <boost/function.hpp>
#include <boost/timer/timer.hpp>
#include <iostream>
#include <string>
#include "kmer.h"
#include "kguts.h"
#include "fasta_parser.h"

using namespace boost::asio::ip;

class KmerLookupClient
{
public:
    typedef std::vector<std::pair<std::string, unsigned int>> result_t;

    KmerLookupClient(std::shared_ptr<KmerGuts> kguts,
		     std::ostream &response_stream,
		     std::istream &input,
		     KmerPegMapping &mapping,
		     boost::function<void ( )> on_completion);

private:
    std::shared_ptr<KmerGuts> kguts_;
    std::istream &input_;
    std::ostream &response_stream_;
    KmerPegMapping &mapping_;
    FastaParser fasta_parser_;
    std::map<KmerPegMapping::encoded_id_t, unsigned int> hit_count_;

    int on_parsed_seq(const std::string &id, const std::string &seq);
    void on_hit(KmerGuts::sig_kmer_t &hit);

    boost::function<void ( )> on_completion_;
    boost::timer::cpu_timer timer_;
};

#endif
