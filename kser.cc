#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#ifdef GPROFILER
#include <gperftools/profiler.h>
#endif

#define DEFINE_GLOBALS 1
#include "global.h"

#include "kmer.h"
#include "kserver.h"
#include "klookup.h"
#include "threadpool.h"

using namespace boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options] listen-port kmer-data-dir\nAllowed options";
    po::options_description desc(x.str());

    std::string listen_port;
    std::string listen_port_file;
    std::string kmer_data;
    std::string peg_kmer_data;
    std::string families_file;
    int n_kmer_threads;

    desc.add_options()
	("help,h", "show this help message")
	("n-kmer-threads", po::value<int>(&n_kmer_threads)->default_value(1), "number of kmer processing threads")
	("listen-port-file", po::value<std::string>(&listen_port_file)->default_value("/dev/null"), "save the listen port to this file")
	("peg-kmer-data", po::value<std::string>(&peg_kmer_data), "precomputed PEG/kmer data file")
	("listen-port,l", po::value<std::string>(&listen_port)->required(), "port to listen on. 0 means to choose a random port")
	("kmer-data-dir,d", po::value<std::string>(&kmer_data)->required(), "kmer data directory")
	("families-file", po::value<std::string>(&families_file), "families file")
	("reserve-mapping", po::value<int>(), "Reserve this much space in global mapping table")
	("no-populate-mmap", po::bool_switch(), "Don't populate mmap data at startup")
	("debug-http", po::bool_switch(), "Debug HTTP protocol")
	;
    po::positional_options_description pd;
    pd.add("listen-port", 1)
	.add("kmer-data-dir", 1);

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

    /*
     * We never really used this. If we want to again we will need to
     * move this code to the global mapping initialization in kserver.cc.
     */
    /*
    if (vm.count("peg-kmer-data"))
    {
	path data_dir(vm["peg-kmer-data"].as<std::string>());

	if (is_directory(data_dir))
	{
	    directory_iterator end_iter;
	    for (auto dit = directory_iterator(data_dir); dit != end_iter; dit++)
	    {
		std::cout << dit->path().string() << "\n";
		mapping.load_compact_mapping_file(dit->path().string());
	    }
	}
	else
	{
	    mapping.load_compact_mapping_file(data_dir.string());
	}
    }
    */
    
    boost::asio::io_service io_service;

    std::shared_ptr<ThreadPool> tp = std::make_shared<ThreadPool>(kmer_data);

    std::shared_ptr<KmerRequestServer> kserver = std::make_shared<KmerRequestServer>(io_service, listen_port, listen_port_file,
										     tp);

    tp->start(n_kmer_threads);

    kserver->startup();

    #ifdef GPROFILER
    std::cout << "profiler enable\n";
    ProfilerStart("prof.out");
    #endif
    io_service.run();

    #ifdef GPROFILER
    ProfilerStop();
    std::cout << "profiler disable\n";
    #endif

    return 0;
}
