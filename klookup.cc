
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
				   boost::function<void ( )> on_completion)
    : kguts_(kguts),
      response_stream_(response_stream),
      input_(input),
      mapping_(mapping),
      on_completion_(on_completion),
      fasta_parser_(std::bind(&KmerLookupClient::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2))
{
    timer_.start();
    fasta_parser_.parse(input);
    
    on_completion_();
}

int KmerLookupClient::on_parsed_seq(const std::string &id, const std::string &seq)
{
    // std::cout << "on parsed " << id << std::endl;
    kguts_->process_aa_seq(id.c_str(), seq.c_str(), seq.size(), 0, std::bind(&KmerLookupClient::on_hit, this, std::placeholders::_1), 0);
    // std::cout << "done\n";

    typedef std::pair<std::string, unsigned int> data_t;

    std::vector<data_t> vec;
    for (auto it = hit_count_.begin(); it != hit_count_.end(); it++)
    {
	std::string peg = mapping_.decode_id(it->first);
	vec.push_back(data_t(peg, it->second));
    }

    std::sort(vec.begin(), vec.end(), less_second<data_t>()); 

    response_stream_ << id << "\n";
    for (auto it = vec.begin(); it != vec.end(); it++)
    {
	response_stream_ << it->first << "\t" << it->second << "\n";
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

