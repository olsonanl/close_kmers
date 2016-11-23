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
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "parallel_read.h"

#include "popen.h"

#include "kmer.h"
#include "kguts.h"
#include "global.h"

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
    next_family_id_(0),
    kcount_(0),
    next_genome_id_(0)
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
	unsigned int id = next_genome_id_++;
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

    abort();
    #if 0
    unsigned int i = 0;
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

inline static void fam_map_insert(std::vector<KmerPegMapping::encoded_family_id_t> &data,
				  KmerPegMapping::encoded_family_id_t fam_id,
				  KmerPegMapping::encoded_kmer_t kmer)
{
    auto ent = std::find(data.begin(), data.end(), fam_id);
    if (ent == data.end())
    {
//	std::cerr << kmer << " add new " << fam_id << "\n";
	data.push_back(fam_id);
    }
//    else
//    {
//	std::cerr << kmer << " have " << fam_id << "\n";
//    }
}

inline static void fam_map_insert(std::unordered_map<KmerPegMapping::encoded_family_id_t, unsigned int> &data,
				  KmerPegMapping::encoded_family_id_t fam_id)
{
    data[fam_id]++;
}

inline static void fam_map_insert(std::unordered_set<KmerPegMapping::encoded_family_id_t> &data,
				  KmerPegMapping::encoded_family_id_t fam_id)
{
    data.insert(fam_id);
}

void KmerPegMapping::add_fam_mapping(KmerPegMapping::encoded_family_id_t fam_id, encoded_kmer_t kmer)
{
    /*
     * kmer_to_family_id_ is a tbb concurrent type so we don't have to lock for the
     * emplace. However, to save memory the contained set is not a concurrent type,
     * so we guard the family insert into the contained set with a mutex.
     */
    
    std::pair<family_map_type_t::iterator, bool> n = kmer_to_family_id_.emplace(std::make_pair(kmer, family_counts_t()));

    // For concurrent map don't need this
    // boost::lock_guard<boost::mutex> guard(mtx_);
    // tbb::spin_mutex::scoped_lock guard(tmtx_);

    family_map_type_t::iterator &fmt = n.first;

    fam_map_insert(fmt->second, fam_id, kmer);

    if (n.second)
    {
	unsigned long val = kcount_++;
	if (val % 1000000 == 0)
	    std::cerr << kmer_to_family_id_.size() << " entries with " << val << " values load-factor " << kmer_to_family_id_.load_factor() << "\n";
    }
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
	encoded_id_t id = assign_new_peg_id(peg);
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
    peg_to_family_.rehash(size_estimate * 10);
    #else
    peg_to_family_.reserve(size_estimate * 10);
    #endif

    f.close();

    int n_threads = (*g_parameters)["n-family-file-threads"].as<int>();
    
    std::map<std::string, bool> warned;
    parallel_read(families_file, n_threads, [this, &warned, zeros, families_file](const std::string &line) {

	    std::vector<std::string> cols;
	    boost::split(cols, line, boost::is_any_of("\t"));

	    std::string pgf("PGF_");
	    pgf += cols[0].substr(2);
	    std::string plf("PLF_");
	    auto mapped_name = genus_map_.find(cols[7]);
	    unsigned long genus_id = 0;
	    if (mapped_name == genus_map_.end())
	    {
		if (!warned[cols[7]])
		{
		    std::cerr << "Cannot map genus '" << cols[7] << "' in " << families_file << "\n";
		    warned[cols[7]] = true;
		}
		plf += cols[7];
	    }
	    else
	    {
		plf += mapped_name->second;
		genus_id = std::stoul(mapped_name->second);
	    }
	    plf += "_";
	    plf += zeros.substr(0, 8 - cols[8].size());
	    plf += cols[8];
	    // encoded_id_t id = encode_id(cols[3]);
	    encoded_id_t id = assign_new_peg_id(cols[3]);
	    
	    // assign unique key for for family
	    family_key_t fkey(pgf, plf);
	    encoded_family_id_t fam_id;
	    
	    unsigned long seqlen = std::stoul(cols[4]);
	    auto fiter = family_key_to_id_.find(fkey);
	    if (fiter == family_key_to_id_.end())
	    {
		/*
		 * First occurrence of this family.
		 * Assign an internal ID and initialize the full data entry.
		 */
		boost::lock_guard<boost::mutex> guard(mtx_);

		/*
		 * Recheck the hash in case someone else came through and updated.
		 */
		fiter = family_key_to_id_.find(fkey);
		if (fiter == family_key_to_id_.end())
		{
		    fam_id = next_family_id_++;
		    family_key_to_id_[fkey] = fam_id;
		    
		    family_data_.emplace(std::make_pair(fam_id, family_data_t { pgf, plf, genus_id, cols[5], fam_id, seqlen, 1 }));
		}
		else
		{
		    std::cerr << "Family map found on second check for " << pgf << " " << plf << "\n";
		    fam_id = fiter->second;
		    auto famit = family_data_.find(fam_id);
		    if (famit == family_data_.end())
		    {
			std::cerr << "Fail: no family data found for " << fam_id << "\n";
			abort();
		    }
		    famit->second.total_size += seqlen;
		    famit->second.count++;
		}
		
//	    std::cout << "assigned " << pgf << " " << plf << ": " << fam_id << "\n";
//	    std::cout << fkey.first << " " << fkey.second << "\n";
	    }
	    else
	    {
		fam_id = fiter->second;
		auto famit = family_data_.find(fam_id);
		if (famit == family_data_.end())
		{
		    std::cerr << "Fail: no family data found for " << fam_id << "\n";
		    abort();
		}
		famit->second.total_size += seqlen;
		famit->second.count++;
	    }
	    
	    // std::cerr << "fmap " << id << " " << fam_id << " " << plf << "\n";
	    
	    peg_to_family_.insert(std::make_pair(id, fam_id));
	    
	});
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

void KmerPegMapping::write_kmer_distribution(std::ostream&os)
{
    char kmer[10];
    for (std::pair<encoded_kmer_t, family_counts_t> entry: kmer_to_family_id_)
    {
	const family_counts_t &counts = entry.second;
	KmerGuts::decoded_kmer(entry.first, kmer);
	os << kmer << "\t" << entry.first << "\t" << counts.size();
	if (counts.size() == 1)
	{
	    //const family_data_t &fam = family_data_[counts.begin()->first];
	    const family_data_t &fam = family_data_[counts[0]];
	    os << "\t" << fam.pgf << "\t" << fam.plf << "\t" << fam.function << "\n";
	}
	else
	{
	    os << "\n";
	}
    }
}
