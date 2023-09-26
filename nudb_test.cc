#include "nudb_kmer_db.h"
#include "kmer_nudb.h"
#include "fasta_parser.h"
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/math/statistics/linear_regression.hpp>

#include <stdexcept>

namespace fs = boost::filesystem;

int main(int argc, char **argv)
{
    typedef NuDBKmerDb<8> KDB;
    KDB db("/dev/shm/out9");

    if (!db.exists())
    {
	std::cerr << "database does not exist\n";
	exit(1);
	// db.create();
    }
    db.open();

    KmerNudb kmer_nudb(db, "/disks/tmp/bob/out9/function.index");

    try {
	for (int i = 1; i < argc; i ++)
	{
	    FastaParser parser;

	    parser.set_callback([&kmer_nudb](const std::string &id, const std::string &seq) {

		double slen = static_cast<double>(seq.length());
		auto calls = std::make_shared<std::vector<KmerCall>>();
		// auto hits = std::make_shared<std::vector<hit_in_sequence_t>>();
		auto ostats = std::make_shared<KmerOtuStats>();

		auto hit_cb = [ slen](const hit_in_sequence_t &hit) {
		    auto kd = hit.kdata;
		    // std::cout << hit.kmer << " " << kd->function_index << " " << kd->mean << " " << "\n";
		};
		kmer_nudb.process_aa_seq(id, seq, calls, hit_cb, ostats);
		for (auto c: *calls)
		{
		    std::cout << c << "\n";
		}
		FunctionIndex fi;
		std::string func;
		float score;
		float wt;
		float offset;
		kmer_nudb.find_best_call(*calls, fi, func, score, wt, offset);
		std::cout << id << "\t" << func << "\t" << fi << "\t" << score << "\n";
		return 0;
	    });

	    fs::path p(argv[i]);
	    fs::ifstream ifstr(p);
	    parser.parse(ifstr);
	    parser.parse_complete();
	}

    }
    catch (std::runtime_error &x)
    {
	std::cerr<< "caught " << x.what() << "\n";
    }
}

