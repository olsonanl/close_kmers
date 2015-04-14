#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include "kmer.h"
#include "kserver.h"
#include "klookup.h"

using namespace boost::filesystem;

int main(int argc, char* argv[])
{
    if (argc != 5)
    {
	std::cout << "Usage: " << argv[0] << " kmer-guts-host kmer-guts-port kmer-data-dir peg-kmer-data\n";
	return 1;
    }

    std::string khost = argv[1];
    std::string kport = argv[2];
    std::string kmer_data = argv[3];
    std::string peg_kmer_data = argv[4];

    KmerPegMapping mapping(kmer_data);

    path data_dir(peg_kmer_data);

    if (is_directory(data_dir))
    {
	directory_iterator end_iter;
	for (auto dit = directory_iterator(data_dir); dit != end_iter; dit++)
	{
	    std::cout << dit->path().string() << "\n";
	    mapping.load_mapping_file(dit->path().string());
	}
    }
    else
    {
	mapping.load_mapping_file(data_dir.string());
    }

    boost::asio::io_service io_service;

    tcp::resolver resolver(io_service);

    auto iter = resolver.resolve({khost, kport});
    boost::asio::ip::tcp::endpoint endpoint = *iter;
    
    KmerRequestServer kserver(io_service, "5100", mapping, endpoint);

//	  std::ifstream ifile(argv[3]);
	  
//	  client c(io_service, argv[1], argv[2], ifile, mapping);
	  io_service.run();

      return 0;
}
