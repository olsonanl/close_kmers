#include "kguts.h"
#include <cctype>
#include <iostream>
#include "fasta_parser.h"

void process_hit(KmerGuts &kg, KmerGuts::sig_kmer_t &x)
{
    char kmer[10];
    kg.decoded_kmer(x.which_kmer, kmer);
    char *fn = kg.kmersH->function_array[x.function_index];
    std::cout << kmer << "\t" << fn << std::endl;
}

int process_seq(const std::string &id, const std::string &seq)
{
    std::cout << "process " << id << std::endl;
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
	std::cerr << "Usage: " << argv[0] << " data-dir\n";
	return 1;
    }

    KmerGuts kg(argv[1]);

    auto x = [&kg](const std::string &id, const std::string &seq)->int {
	    auto calls = std::make_shared<std::vector<KmerCall>>();
	    auto otu_stats = std::make_shared<KmerOtuStats>();
	    auto hit_cb = [&kg](KmerGuts::sig_kmer_t &x)->void {
		process_hit(kg, x);
	    };

	    kg.process_aa_seq(id, seq, calls, 0, otu_stats);

	    std::cout << id << std::endl;
	    for (auto x = calls->begin(); x != calls->end(); x++)
	    {
		std::cout << "\t" << x->start  << "\t" << x->end << "\t" << x->count << "\t" << kg.kmersH->function_array[x->function_index] << std::endl;
	    }
	    auto end = std::next(otu_stats->otus_by_count.begin(), std::min(otu_stats->otus_by_count.size(), (size_t) 10));
	    for (auto x = otu_stats->otus_by_count.begin(); x != end; x++)
	    {
		std::cout << x->first << "\t" << x->second << "\t" << kg.kmersH->otu_array[x->first] << std::endl;
	    }
    };

    FastaParser fp(x);
    fp.parse(std::cin);
}
