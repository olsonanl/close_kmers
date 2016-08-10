
#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include "klookup.h"
#include "global.h"

using namespace boost::filesystem;
using boost::asio::ip::tcp;

template<class T>
struct less_second : std::binary_function<T,T,bool>
{
    inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
};  


KmerLookupClient::KmerLookupClient(std::shared_ptr<KmerGuts> kguts,
				   std::ostream &response_stream,
				   std::istream &input,
				   KmerPegMapping &mapping,
				   boost::function<void ( )> on_completion,
				   std::map<std::string, std::string> &parameters)
    : kguts_(kguts),
      parameters_(parameters),
      response_stream_(response_stream),
      input_(input),
      mapping_(mapping),
      on_completion_(on_completion),
      fasta_parser_(std::bind(&KmerLookupClient::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2)),
      kmer_hit_threshold_(100)
{
    try {
	kmer_hit_threshold_ = std::stoi(parameters_["kmer_hit_threhsold"]);
    } catch (const std::invalid_argument& ia)
    {
    }


    timer_.start();
    fasta_parser_.parse(input);
    
    on_completion_();
}

int KmerLookupClient::on_parsed_seq(const std::string &id, const std::string &seq)
{
    // std::cout << "on parsed " << id << std::endl;
    kguts_->process_aa_seq(id, seq, 0, std::bind(&KmerLookupClient::on_hit, this, std::placeholders::_1), 0);
    // std::cout << "done\n";

    typedef std::pair<KmerPegMapping::encoded_id_t, unsigned int> data_t;

    std::vector<data_t> vec;
    for (auto it = hit_count_.begin(); it != hit_count_.end(); it++)
    {
	vec.push_back(*it);
    }

    std::sort(vec.begin(), vec.end(), less_second<data_t>()); 

    response_stream_ << id << "\n";
    for (auto it = vec.begin(); it != vec.end(); it++)
    {
	if (it->second < kmer_hit_threshold_)
	    break;
	std::string peg = mapping_.decode_id(it->first);
	response_stream_ << peg << "\t" << it->second;
	auto fhit = mapping_.family_mapping_.find(it->first);

	if (fhit != mapping_.family_mapping_.end())
	{
	    response_stream_ << "\t" << fhit->second.pgf << "\t" << fhit->second.plf << "\t" << fhit->second.function;
	}
	response_stream_ << "\n";
    }
    response_stream_ << "//\n";

    hit_count_.clear();

    return 0;
}

void KmerLookupClient::on_hit(KmerGuts::sig_kmer_t &hit)
{
    // std::cout << "on hit " << hit.which_kmer << " " << hit.function_index << "\n";

    auto ki = mapping_.kmer_to_id_.find(hit.which_kmer);
    if (ki != mapping_.kmer_to_id_.end())
    {
	// std::cout << "got mapping for " << hit.which_kmer << "\n";
	for (auto it = ki->second.begin(); it != ki->second.end(); it++)
	{
	    // std::cout << "  " << *it << "\n";

	    hit_count_[*it]++;
	}
    }
}

