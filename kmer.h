#ifndef _KMER_H
#define _KMER_H

/*
 * Kmer to peg mapping database
 */

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>

#ifdef USE_TBB
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#endif

#include <boost/thread/mutex.hpp>

class KmerPegMapping
{
public:
    KmerPegMapping();
    ~KmerPegMapping();

    void reserve_mapping_space(size_t n);

    void load_genome_map(const std::string &genome_mapping_file);

    void load_mapping_file(const std::string &mapping_file);
    void load_compact_mapping_file(const std::string &mapping_file);
    void load_compact_mapping_file(std::istream &kfile);

    typedef unsigned int encoded_family_id_t;
    typedef unsigned int encoded_id_t;
    typedef unsigned long encoded_kmer_t;


    /*
     * We maintain the following data structures for holding the family data.
     *
     * encoded_family_id_t is the numeric index assigned to each (local) family.
     * family_key_t is a pair of pgf,plf that identifies a unique family.
     * family_to_id_map_t is a mapping from family_key_t to encoded_family_id_t
     * 
     * The full family data is indexed on encoded_family_id_t; it is held in the family_data
     * struct. 
     * We use id_to_family_map_t to store the mapping from encoded family ID to full family data.
     * We use peg_to_family_map_t to store the mapping from a constituent peg ID to the family ID.
     *
     */

    struct family_data_t {
	std::string pgf;
	std::string plf;
	unsigned long genus_id;
	std::string function;
	encoded_family_id_t family_id;
	unsigned long total_size; /* in aa */
	unsigned short count;
//    family_data(const std::string a, const std::string b, const std::string c, encoded_family_id_t d)
//    : pgf(a), plf(b), function(c), family_id(d), total_size(0) {};
    };
    typedef std::pair<std::string, std::string> family_key_t;
    
    /*
     * Maps string genus name from family file
     * to the taxon id used in PLF identifiers.
     */
    std::map<std::string, std::string> genus_map_;

#ifdef USE_TBB
    typedef std::unordered_map<encoded_family_id_t, unsigned int> family_counts_t;
    //typedef std::unordered_set<encoded_family_id_t> family_counts_t;
    //typedef tbb::concurrent_vector<encoded_family_id_t> family_counts_t;
    //typedef tbb::concurrent_unordered_map<encoded_family_id_t, unsigned long> family_counts_t;
    typedef tbb::concurrent_vector<encoded_id_t> id_set;
    typedef tbb::concurrent_unordered_map<encoded_kmer_t, id_set> map_type_t;
    typedef tbb::concurrent_unordered_map<encoded_kmer_t, family_counts_t> family_map_type_t;
    typedef tbb::concurrent_unordered_map<std::string, encoded_id_t> genome_to_id_map_t;
    typedef tbb::concurrent_unordered_map<encoded_id_t, std::string> id_to_genome_map_t;
    typedef tbb::concurrent_unordered_map<std::string, encoded_id_t> peg_to_id_map_t;
    typedef tbb::concurrent_unordered_map<encoded_id_t, std::string> id_to_peg_map_t;

    typedef tbb::concurrent_unordered_map<encoded_family_id_t, family_data_t> family_id_to_family_map_t;
    typedef tbb::concurrent_unordered_map<family_key_t, encoded_family_id_t> family_to_family_id_map_t;
    typedef tbb::concurrent_unordered_map<encoded_id_t, encoded_family_id_t> peg_to_family_map_t;
#else
    typedef std::vector<encoded_id_t> id_set;
    typedef std::unordered_map<encoded_id_t, id_set> map_type_t;
    typedef std::map<std::string, encoded_id_t> genome_to_id_map_t;
    typedef std::map<encoded_id_t, std::string> id_to_genome_map_t;
    typedef std::unordered_map<encoded_id_t, family_data> family_map_t;
#endif

    boost::mutex mtx_;

    // Peg ID mapping
    encoded_id_t next_peg_id_;
    peg_to_id_map_t peg_to_id_;
    id_to_peg_map_t id_to_peg_;

    // Family ID mapping
    encoded_family_id_t next_family_id_;
    family_id_to_family_map_t family_data_;
    family_to_family_id_map_t family_key_to_id_;
    peg_to_family_map_t peg_to_family_;
    
    map_type_t kmer_to_id_;
    family_map_type_t kmer_to_family_id_;

    genome_to_id_map_t genome_to_id_;
    id_to_genome_map_t id_to_genome_;

    void write_kmer_distribution(std::ostream &os);

    void load_genus_map(const std::string &genus_file);
    void load_families(const std::string &families_file);
    const std::string& lookup_genus(const std::string &genus) { return genus_map_[genus]; }

    /*
     * We also maintain a mapping from an md5 value to the pegs with that
     * md5.
     */

    std::map<std::array<char, 16>, std::vector<encoded_id_t> > md5_to_peg_;

    encoded_id_t encode_id(const std::string &peg);
    #if 0
    encoded_id_t encode_id(const std::string &genome, const std::string &peg);
    #endif

    void add_mapping(encoded_id_t enc, encoded_kmer_t kmer);
    void add_fam_mapping(encoded_family_id_t enc, encoded_kmer_t kmer);

    std::string decode_id(encoded_id_t id);

    unsigned long kcount_;
    unsigned int next_genome_id_;

    void dump_sizes(std::ostream &os);
};


#endif
