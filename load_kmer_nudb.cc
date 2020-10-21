#include <nudb/nudb.hpp>
#include <cstddef>
#include <cstdint>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <vector>

#include "kmer_encoder.h"
#include "kmer_types.h"

struct KData {
    OTUIndex otu_index;
    unsigned short  avg_from_end;
    FunctionIndex function_index;
    float function_wt;
};

namespace fs = boost::filesystem;

int main()
{
    using key_type = unsigned long long;
    nudb::error_code ec;
    auto const dat_path = "db.dat";
    auto const key_path = "db.key";
    auto const log_path = "db.log";
    nudb::create<nudb::xxhasher>(
	dat_path, key_path, log_path,
	1,
	nudb::make_salt(),
	sizeof(key_type),
	nudb::block_size("."),
	0.5f,
	ec);
    nudb::store db;
    db.open(dat_path, key_path, log_path, ec);
    if (ec)
    {
	std::cerr << "open failed: " << ec << "\n";
	exit(1);
    }

    fs::path file = "/vol/core-seed/kmers/core.2020-0417/Data.2/final.kmers";
    fs::ifstream f(file);
    if (f.fail())
    {
	throw std::runtime_error("kmer file open " + file.string() + " failed");
    }
    
    unsigned long line_number = 0;
    KmerEncoder encoder;
    
    std::string line;
    std::vector<std::string> cols;
    KData kdata;
    try {
	while (std::getline(f, line))
	{
	    line_number++;
	    boost::split(cols, line, boost::is_any_of("\t"));

	    unsigned long long enc = encoder.encoded_aa_kmer(cols[0]);
	    kdata.avg_from_end = std::stoi(cols[1]);
	    kdata.function_index = std::stoi(cols[2]);
	    kdata.function_wt = std::stof(cols[3]);
	    kdata.otu_index = 0;
	    if (!cols[4].empty())
		kdata.otu_index = std::stoi(cols[4]);

	    db.insert(&enc, &kdata, sizeof(kdata), ec);
	    if (ec)
	    {
		std::cerr << "insert failed: " << ec << "\n";
		exit(1);
	    }
	    //if (line_number > 100)
	    //break;
	}
    }
    catch (std::exception &e)
    {
	std::cerr << "Error loading " << file.string() << ": " << e.what() << "\n";
    }


    db.close(ec);
    if (ec)
    {
	std::cerr << "close failed: " << ec << "\n";
	exit(1);
    }
}
