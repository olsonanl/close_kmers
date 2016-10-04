#include <stdlib.h>
#include <unordered_map>
#include <array>
#include <string>
#include <map>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

#include "popen.h"

#include "kmer.h"

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
	elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

KmerPegMapping::KmerPegMapping() :
    next_genome_id_(0),
    next_family_id_(0),
    kcount_(0)
{
    std::cout << "Constructed KmerPegMapping\n";
}

void KmerPegMapping::reserve_mapping_space(size_t n)
{
    #ifdef USE_TBB
    //kmer_to_id_.rehash(n);
    kmer_to_family_id_.rehash(n);
    #else
    kmer_to_id_.reserve(n);
    #endif
}
					 

void KmerPegMapping::load_genome_map(const std::string &genome_file)
{
    std::ifstream gfile(genome_file);
    if (gfile.fail())
    {
	std::cout << "Failed to open " << genome_file << "\n";
	exit(1);
    }

    std::string line;
    while (std::getline(gfile, line))
    {
	std::vector<std::string> x = split(line, '\t');
	std::string genome = x[1];
	unsigned long id = next_genome_id_++;
	genome_to_id_[genome] = id;
	id_to_genome_[id] = genome;
    }

    gfile.close();

}

KmerPegMapping::~KmerPegMapping()
{
    std::cout << "Destroyed KmerPegMapping\n";
}


void KmerPegMapping::load_mapping_file(const std::string &mapping_file)
{
    /*
     * Read and encode the peg/kmer data.
     */
    unsigned int i = 0;

    abort();
    #if 0
    std::ifstream kfile(mapping_file);

    std::string line;
    while (std::getline(kfile, line))
    {
	size_t i1 = line.find('\t');
	unsigned long kmer = std::stoul(line.substr(0, i1));
	i1 = line.find('|', i1 + 1);
	size_t i2 = line.find('p', i1 + 1);
	i1 = line.find('.', i2);
	i2 = line.find('\t', i1);

	unsigned long enc = encode_id(line.substr(i1 + 1, i2 - i1 - 2),
				      line.substr(i1 + 1, i2 - i1 - 1));
	
	//printf("k=%lu g=%lu id=%lu enc=%lu\n", kmer, gid, id, enc);
	kmer_to_id_[kmer].push_back(enc);
	i++;
	if (i % 1000000 == 0)
	    printf("%d\n", i);

    }
    kfile.close();
    #endif
}

void KmerPegMapping::load_compact_mapping_file(const std::string &mapping_file)
{
    /*
     * Read and encode the peg/kmer data.
     */
    unsigned int i = 0;

    if (mapping_file.substr(mapping_file.length() - 3) == ".gz")
    {
	Popen pipe_wrapper(("gunzip < " + mapping_file).c_str());
	io::file_descriptor_source pipe_device(fileno(pipe_wrapper.stream()),
					       io::never_close_handle);
	io::stream<io::file_descriptor_source> pipe_stream(pipe_device, 0x1000, 0x1000);
	load_compact_mapping_file(pipe_stream);
    }
    else
    {
	std::ifstream kfile(mapping_file);
	load_compact_mapping_file(kfile);
	kfile.close();
    }
}

void KmerPegMapping::load_compact_mapping_file(std::istream &kfile)
{
    unsigned long kmer;
    id_set vals;

    while (kfile >> kmer)
    {
	unsigned long size;
	kfile >> size;

	vals.reserve(size);
	
	while (kfile.peek() != '\n')
	{
	    std::string fid;
	    kfile >> fid;

	    encoded_id_t enc = encode_id(fid);
	    vals.push_back(enc);
	}
	// std::cout << kmer << " " << size << " got n values " << vals.size() << "\n";
	kmer_to_id_[kmer] = vals;
	vals.clear();
    }
}

void KmerPegMapping::add_mapping(KmerPegMapping::encoded_id_t enc, unsigned long kmer)
{
    #ifdef USE_TBB

    // kmer_to_id_[kmer].push_back(enc);

    auto it = kmer_to_id_.find(kmer);

    if (it == kmer_to_id_.end())
    {
	auto n = kmer_to_id_.emplace(std::make_pair(kmer, id_set()));
	//std::cout << "Alloc new for " << kmer << "\n";
	n.first->second.reserve(8);
	n.first->second.push_back(enc);
    }
    else
    {
	//std::cout << "reuse " << kmer << "\n";
	it->second.push_back(enc);
    }
    #else
    auto it = kmer_to_id_.find(kmer);

    if (it == kmer_to_id_.end())
    {
	auto n = kmer_to_id_.emplace(std::make_pair(kmer, id_set()));
	//std::cout << "Alloc new for " << kmer << "\n";
	n.first->second.reserve(2);
	n.first->second.push_back(enc);
    }
    else
    {
	//std::cout << "reuse " << kmer << "\n";
	it->second.push_back(enc);
    }
    #endif
    kcount_++;
    if (kcount_ % 1000000 == 0)
	std::cerr << kmer_to_id_.size() << " entries with " << kcount_ << " values load-factor " << kmer_to_id_.load_factor() << "\n";
    
}

void KmerPegMapping::add_fam_mapping(KmerPegMapping::encoded_id_t enc, encoded_kmer_t kmer)
{

    /*
     * use id_to_family_ to determine the family mapped to this id, and store
     * the kmer relationship.
     */

    auto fam_id_iter = family_mapping_.find(enc);
    if (fam_id_iter == family_mapping_.end())
    {
	std::string dec = decode_id(enc);
	std::cout << "NO FAM FOR " << enc << " " << dec << "\n";
	return;
    }
    encoded_family_id_t fam_id = fam_id_iter->second.family_id;
    
    auto it = kmer_to_family_id_.find(kmer);

    if (it == kmer_to_family_id_.end())
    {
	auto n = kmer_to_family_id_.emplace(std::make_pair(kmer, family_counts_t()));
	// std::cout << "Alloc new for " << kmer << "\n";
	n.first->second[fam_id] = 1;
    }
    else
    {
	// std::cout << "reuse " << kmer << "\n";
	it->second[fam_id]++;
    }
    kcount_++;
    if (kcount_ % 1000000 == 0)
	std::cerr << kmer_to_family_id_.size() << " entries with " << kcount_ << " values load-factor " << kmer_to_family_id_.load_factor() << "\n";
    
}

#if 1

KmerPegMapping::encoded_id_t KmerPegMapping::encode_id(const std::string &peg)
{
    auto x = peg_to_id_.find(peg);
    if (x != peg_to_id_.end())
    {
	encoded_id_t id = x->second;
	// std::cerr << "found " << peg << " " << id << "\n";
	return id;
    }
    else
    {
	boost::lock_guard<boost::mutex> guard(mtx_);
	
	encoded_id_t id = next_peg_id_++;
	peg_to_id_[peg] = id;
	id_to_peg_[id] = peg;
	if (id % 1000000 == 0)
	    std::cerr << "assigned " << peg << " " << id << "\n";
	
	return id;
    }
}

std::string KmerPegMapping::decode_id(encoded_id_t id)
{
    auto x = id_to_peg_.find(id);
    if (x != id_to_peg_.end())
	return x->second;
    else
	return "";
}

#else

KmerPegMapping::encoded_id_t KmerPegMapping::encode_id(const std::string &peg)
{
    size_t i1 = peg.find('|');
    size_t i2 = peg.find('p', i1 + 1);
    std::string genome = peg.substr(i1 + 1, i2 - i1 - 2);
    i1 = peg.find('.', i2);
    std::string fid = peg.substr(i1 + 1);
    // std::cout << "peg='" << peg << "' genome='" << genome << "' id='" << fid << "'\n";
    return encode_id(genome, fid);
}

KmerPegMapping::encoded_id_t KmerPegMapping::encode_id(const std::string &genome, const std::string &peg)
{
    unsigned long gid;
    auto it = genome_to_id_.find(genome);
    if (it == genome_to_id_.end())
    {
	gid = next_genome_id_++;
	genome_to_id_[genome] = gid;
	id_to_genome_[gid] = genome;
    }
    else
    {
	gid = it->second;
    }
    // unsigned long gid = genome_to_id_[genome];

    return (gid << 17) | (std::stoul(peg));
}

std::string KmerPegMapping::decode_id(encoded_id_t id)
{
    std::string g = id_to_genome_[id >> 17];
    unsigned int peg = id & 0x1ffff;
    return "fig|" + g + ".peg." + std::to_string((long long) peg);
}

#endif

void KmerPegMapping::load_genus_map(const std::string &genus_file)
{
    std::ifstream gf(genus_file);

    if (gf.fail())
    {
	std::cerr << "Error opening gnus file " << genus_file << "\n";
	exit(1);
    }

    std::string line;
    while (std::getline(gf, line))
    {
	std::vector<std::string> cols;
	boost::split(cols, line, boost::is_any_of("\t"));
	genus_map_[cols[0]] = cols[1];
    }
}


/*
 * Load a PATRIC families global-fams file and set up the global and local family attributes.
 *
 * Columns:
 * 0 	global family
 * 1	fams merged
 * 2 	genera merged
 * 3	peg id
 * 4	protein length
 * 5 	family function
 * 6	local family number
 * 7	genus
 * 8 	local family number again
 *
 * GF00000000	28	21	fig|1049789.4.peg.4658	250	(2-pyrone-4,6-)dicarboxylic acid hydrolase	11358	Leptospira 11358
 *
 */
void KmerPegMapping::load_families(const std::string &families_file)
{
    std::ifstream f(families_file);
    std::string line;
    std::string zeros("00000000");

    if (f.fail())
    {
	std::cerr << "Failure opening families file " << families_file << "\n";
	exit(1);
    }


    /*
     * Read a sample of lines to estimate number of lines in file.
     */
    int lines = 0;
    int n = 100;
    size_t sample_size = 0;
    while (std::getline(f, line) && lines++ < n)
    {
	sample_size += line.size();
    }

    f.seekg(0, f.end);
    size_t fsize = f.tellg();
    f.seekg(0, f.beg);

    // fsize / (sample_size / lines )
       
    size_t size_estimate = fsize * lines / sample_size;
    std::cerr << "fsize=" << fsize << " line estimate=" << size_estimate << "\n";
    
    #ifdef USE_TBB
    family_mapping_.rehash(size_estimate * 10);
    #else
    family_mapping_.reserve(size_estimate * 10);
    #endif

    unsigned int lineno = 0;
    std::map<std::string, bool> warned;
    while (std::getline(f, line))
    {
	lineno++;
	std::vector<std::string> cols;
	boost::split(cols, line, boost::is_any_of("\t"));

	std::string pgf("PGF_");
	pgf += cols[0].substr(2);
	std::string plf("PLF_");
	auto mapped_name = genus_map_.find(cols[7]);
	if (mapped_name == genus_map_.end())
	{
	    if (!warned[cols[7]])
	    {
		std::cerr << "Cannot map genus '" << cols[7] << "' at line " << lineno << " of " << families_file << "\n";
		warned[cols[7]] = true;
	    }
	    plf += cols[7];
	}
	else
	{
	    plf += mapped_name->second;
	}
	plf += "_";
	plf += zeros.substr(0, 8 - cols[8].size());
	plf += cols[8];
	encoded_id_t id = encode_id(cols[3]);

	// assign unique key for for family
	family_key_t fkey(pgf, plf);
	encoded_family_id_t fam_id;
	
	auto fiter = family_to_id_.find(fkey);
	if (fiter == family_to_id_.end())
	{
	    boost::lock_guard<boost::mutex> guard(mtx_);

	    fam_id = next_family_id_++;
	    family_to_id_[fkey] = fam_id;
	    id_to_family_[fam_id] = fkey;
//	    std::cout << "assigned " << pgf << " " << plf << ": " << fam_id << "\n";
//	    std::cout << fkey.first << " " << fkey.second << "\n";
	}
	else
	{
	    fam_id = fiter->second;
	}

	// std::cerr << "fmap " << id << " " << fam_id << " " << plf << "\n";
	family_mapping_.emplace(std::make_pair(id, family_data(pgf, plf, cols[5], fam_id)));
    }
}

void KmerPegMapping::dump_sizes(std::ostream &os)
{
    os << "kmer_to_id_: size=" << kmer_to_id_.size() << "\n";
    size_t bytes = 0;
    for (auto it: kmer_to_id_)
    {
	bytes += it.second.size();
    }
    os << "kmer_to_id_: content size=" << bytes << "\n";
    os << "peg_to_id_: size=" << peg_to_id_.size() << "\n";
    os << "id_to_peg_: size=" << id_to_peg_.size() << "\n";
    os << "genome_to_id_: size=" << genome_to_id_.size() << "\n";
    os << "id_to_genome_: size=" << id_to_genome_.size() << "\n";
    
}
