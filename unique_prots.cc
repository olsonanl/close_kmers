#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#define DEFINE_GLOBALS 1
#include "global.h"

#include <unistd.h>

#include "kguts.h"
#include "fasta_parser.h"

using namespace boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options] kmer-data fasta\nAllowed options";
    po::options_description desc(x.str());

    std::string kmer_data;
    std::string fasta;
    
    desc.add_options()
	("help,h", "show this help message")
	("kmer-data-dir,d", po::value<std::string>(&kmer_data)->required(), "kmer data directory")
	("fasta-file", po::value<std::string>(&fasta), "input fasta file")
	;
    po::positional_options_description pd;
    pd.add("kmer-data-dir", 1)
	.add("fasta-file", 1);

    po::variables_map vm;
    g_parameters = &vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(pd).run(), vm);

    if (vm.count("help"))
    {
	std::cout << desc << "\n";
	return 1;
    }

    try {
	po::notify(vm);
    } catch (po::required_option &e)
    {
	std::cerr << "Invalid command line: " << e.what() << "\n";
	std::cerr << desc << "\n";
	exit(1);
    }

    auto image = std::make_shared<KmerImage>(kmer_data);
    KmerGuts kguts(kmer_data, image);

    /*
     * Parse input fasta and compute kmers present. Store kmer vector => list of ID mapping.
     */
    
    //typedef std::vector<unsigned long long> klist_t;
    typedef std::set<unsigned long long> klist_t;
    typedef std::vector<std::string> peglist_t;
    typedef std::map<klist_t, peglist_t> kmap_t;

    kmap_t kmap;

    FastaParser parser;

    parser.set_callback([&kguts, &kmap](const std::string &id, const std::string &seq) {
	    klist_t klist;
	    auto cb = [&kguts, &klist, id, seq](KmerGuts::hit_in_sequence_t hit) {
		klist.insert(hit.hit.which_kmer);
		// klist.push_back(hit.hit.which_kmer);
	    };

	    kguts.process_aa_seq(id, seq, 0, cb, 0);
	    kmap_t::iterator ent = kmap.find(klist);
	    if (ent == kmap.end())
	    {
		auto x = kmap.emplace(std::make_pair(klist, peglist_t()));
		x.first->second.push_back(id);
	    }
	    else
	    {
		ent->second.push_back(id);
	    }
	    return 0;
	});

    std::ifstream inp(fasta);
    parser.parse(inp);

    for (auto ent: kmap)
    {
	for (auto n: ent.second)
	{
	    std::cout << n << "\t";
	}
	std::cout << "\n";
    }
	 

    return 0;
}
