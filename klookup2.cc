
#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include "klookup2.h"
#include "global.h"

using namespace boost::filesystem;
using boost::asio::ip::tcp;

KmerLookupClient2::KmerLookupClient2(std::shared_ptr<KmerGuts> kguts,
				     std::istream &input,
				     boost::function<void ( const std::string &prot, size_t len )> on_protein,
				     boost::function<void ( unsigned long kmer )> on_hit,
				     boost::function<void ( const std::string &line )> on_call,
				     boost::function<void ( const boost::system::error_code& err )> on_completion)
    : input_(input),
      on_protein_(on_protein),
      on_hit_(on_hit),
      on_call_(on_call),
      on_completion_(on_completion),
      kguts_(kguts),
      fasta_parser_(std::bind(&KmerLookupClient2::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2))
{
    timer_.start();
    fasta_parser_.parse(input);
    boost::system::error_code ec;
    on_completion_(ec);
}

int KmerLookupClient2::on_parsed_seq(const std::string &id, const std::string &seq)
{
    std::cout << "on parsed " << id << std::endl;

    if (on_protein_)
	on_protein_(id, seq.size());
    if (on_call_)
    {
	std::ostringstream oss;
	oss << "PROTEIN_ID\t" << id << "\t" << seq.size();
	on_call_(oss.str());
    }

    auto calls = std::make_shared<std::vector<KmerCall>>();
    kguts_->process_aa_seq(id.c_str(), seq.c_str(), seq.size(), calls,
			   std::bind(&KmerLookupClient2::on_hit, this, std::placeholders::_1), 0);

    if (on_call_)
    {
	for (auto it = calls->begin(); it != calls->end(); it++)
	{
	    std::ostringstream oss;
	    oss << "CALL\t" << it->start << "\t" << it->end << "\t" << "\t" << it->count;
	    oss << "\t" << it->function_index << "\t" << kguts_->kmersH->function_array[it->function_index];
	    on_call_(oss.str());
	}
    }

    return 0;
}

void KmerLookupClient2::on_hit(KmerGuts::sig_kmer_t &hit)
{
    on_hit_(hit.which_kmer);
}

