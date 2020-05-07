/*
 * Little test program to exercise FamilyMapper using fasta data.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include "family_mapper.h"
#include "kguts.h"
#include "kmer_inserter.h"
#include "nr_loader.h"
#include "threadpool.h"
#include "fasta_parser.h"

#define DEFINE_GLOBALS
#include "global.h"

namespace po = boost::program_options;

static std::string kmer_data;
static std::string fasta_file;
static po::variables_map vm;

static int n_load_threads;
static int n_kmer_threads;

static void process_parameters(int argc, char **argv);

int main(int argc, char **argv)
{
    process_parameters(argc, argv);

    std::cout << "testing with kmer_dir=" << kmer_data << " fasta_file=" << fasta_file << "\n";

    std::ifstream in(fasta_file);
    if (!in.good())
    {
	std::cerr << "error opening " << fasta_file << "\n";
	exit(1);
    }

    auto image = std::make_shared<KmerImage>(kmer_data);
    //auto kguts = std::make_shared<KmerGuts>(kmer_data, image);
    auto kguts = new KmerGuts(kmer_data, image);

    /*
     * Create a mapping. This is a copy of code from kserver.cc.
     */

    auto mapping = std::make_shared<KmerPegMapping>();
    mapping->load_genome_map(kmer_data + "/genomes");
    if (g_parameters->count("families-genus-mapping"))
    {
	std::string mapfile = (*g_parameters)["families-genus-mapping"].as<std::string>();
	mapping->load_genus_map(mapfile);
    }

    if (g_parameters->count("families-file"))
    {
	std::string ff = (*g_parameters)["families-file"].as<std::string>();

	std::cerr << "Loading (immediate) families from " << ff << "...\n";
	mapping->load_families(ff);
	std::cerr << "Loading families from " << ff << "... done\n";
    }

    if (g_parameters->count("reserve-mapping"))
    {
	unsigned long count = (*g_parameters)["reserve-mapping"].as<unsigned long>();
	std::cerr << "Reserving " << count << " bytes in mapping table\n";
	mapping->reserve_mapping_space(count);
    }

    int n_inserters = 1;
    if (g_parameters->count("n-inserter-threads"))
    {
	n_inserters = (*g_parameters)["n-inserter-threads"].as<int>();
    }

    KmerInserter inserter(n_inserters, mapping);

    inserter.start();

    std::vector<NRLoader *> active_loaders;
    NRLoadState load_state("master");

    std::shared_ptr<ThreadPool> tp = std::make_shared<ThreadPool>(kmer_data);

    tp->start(n_load_threads);

    if (g_parameters->count("families-nr"))
    {
	auto files = (*g_parameters)["families-nr"].as<std::vector<std::string> >();
	for (auto file: files)
	{
	    load_state.pending_inc();
	    std::cerr << "Queue load NR file " << file << "\n";
	    NRLoader *loader = new NRLoader(load_state, file, mapping, tp, files.size(), true, inserter);
	    active_loaders.push_back(loader);
	    loader->start();
	}
    }
	
    std::cerr << "wait for threads to finish\n";

    load_state.pending_wait();
    
    std::cerr << "done waiting\n";

    for (auto v: active_loaders)
    {
	// std::cerr << "remove loader for " << v->file_ << "\n";
	delete v;
    }

    /*
     * When we are done, we need to clear the inserters.
     */
    inserter.stop();

    auto mapper = FamilyMapper(kguts, mapping);

    FastaParser fp;

    fp.set_callback([&mapper](const std::string &id, const std::string &seq) {
	    mapper.find_all_matches(std::cout, id, seq);
	    return 0;
	});
    fp.parse(in);
    fp.parse_complete();
    
    // auto res = mapper.find_best_family_match(id, seq);


    std::cout << "stop thread pool...\n";
    tp->stop();
    std::cout << "stop thread pool...done\n";
    return 0;
}

static void process_parameters(int argc, char **argv)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options] kmer-data-dir fasta-file\nAllowed options";

    po::options_description desc(x.str());

    
    desc.add_options()
	("help,h", "show this help message")
	("n-family-file-threads", po::value<int>()->default_value(1), "number of family file reader threads")
	("n-inserter-threads", po::value<int>()->default_value(1), "number of kmer inserter threads")
	("n-load-threads", po::value<int>(&n_load_threads)->default_value(1), "number of NR load threads")
	("n-kmer-threads", po::value<int>(&n_kmer_threads)->default_value(1), "number of kmer processing threads")
	("families-genus-mapping", po::value<std::string>(), "genus name to taxid mapping file")
	("families-file", po::value<std::string>(), "families file")
	("families-nr", po::value<std::vector<std::string>>()->multitoken(), "families NR data")
	("kmer-data-dir,d", po::value<std::string>(&kmer_data)->required(), "kmer data directory")
	("fasta-file", po::value<std::string>(&fasta_file), "input fasta file")
	("reserve-mapping", po::value<unsigned long>(), "Reserve this much space in global mapping table")
	("no-populate-mmap", po::bool_switch(), "Don't populate mmap data at startup")
	;
    po::positional_options_description pd;
    pd
	.add("kmer-data-dir", 1)
	.add("fasta-file", 1);

    g_parameters = &vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(pd).run(), vm);

    if (vm.count("help"))
    {
	std::cout << desc << "\n";
	exit(0);
    }

    try {
	po::notify(vm);
    } catch (po::required_option &e)
    {
	std::cerr << "Invalid command line: " << e.what() << "\n";
	std::cerr << desc << "\n";
	exit(1);
    }
}
