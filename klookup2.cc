
#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <memory>
#include <boost/asio.hpp>
#include <boost/bind/bind.hpp>
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
				     boost::function<void ( const boost::system::error_code& err )> on_completion,
				     std::map<std::string, std::string> &parameters)
    : input_(input),
      parameters_(parameters),
      on_protein_(on_protein),
      on_hit_(on_hit),
      on_call_(on_call),
      on_completion_(on_completion),
      kguts_(kguts),
      fasta_parser_(std::bind(&KmerLookupClient2::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2))
{
    silent_ = 0;
    try {
	silent_ = std::stoi(parameters_["silent"]);
    } catch (const std::invalid_argument& ia)
    {
    }

    timer_.start();
    fasta_parser_.parse(input);
    boost::system::error_code ec;
    on_completion_(ec);
}

int KmerLookupClient2::on_parsed_seq(const std::string &id, const std::string &seq)
{
    // std::cout << "on parsed " << id << std::endl;

    if (on_protein_)
	on_protein_(id, seq.size());
    if (on_call_ && !silent_)
    {
	std::ostringstream oss;
	oss << "PROTEIN-ID\t" << id << "\t" << seq.size();
	on_call_(oss.str());
    }

    std::shared_ptr<std::vector<KmerCall>> calls = 0;
    std::shared_ptr<KmerOtuStats> otu_stats = 0;
    if (on_call_ && !silent_)
    {
	otu_stats = std::make_shared<KmerOtuStats>();
	calls = std::make_shared<std::vector<KmerCall>>();
    }

    if (on_hit_)
    {    
	kguts_->process_aa_seq(id, seq, calls, 
			       std::bind(&KmerLookupClient2::on_hit, this, std::placeholders::_1), otu_stats);
    }
    else
    {
	kguts_->process_aa_seq(id, seq, calls, 0, otu_stats);
    }	

    if (on_call_ && !silent_)
    {
	for (auto it = calls->begin(); it != calls->end(); it++)
	{
	    std::ostringstream oss;
	    oss << "CALL\t" << it->start << "\t" << it->end << "\t" << it->count;
	    oss << "\t" << it->function_index << "\t" << kguts_->kmersH->function_array[it->function_index];
	    oss << "\t" << it->weighted_hits;
	    on_call_(oss.str());
	}
	std::ostringstream oss;
	oss << "OTU-COUNTS\t" << id << "[" << seq.size() << "]";

	auto end = std::next(otu_stats->otus_by_count.begin(), std::min(otu_stats->otus_by_count.size(), (size_t) 5));
	for (auto x = otu_stats->otus_by_count.begin(); x != end; x++)
	{
	    oss << "\t" << x->second << "-" << x->first;
	}
	on_call_(oss.str());
    }

    return 0;
}

void KmerLookupClient2::on_hit(KmerGuts::sig_kmer_t &hit)
{
    on_hit_(hit.which_kmer);
}

