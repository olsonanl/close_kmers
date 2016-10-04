#ifndef _KMER_H
#define _KMER_H

/*
 * Kmer to peg mapping database
 */

#include <unordered_map>
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

    /* peg to peg-attributes mapping. */
    struct family_data {
	std::string pgf;
	std::string plf;
	std::string function;
	encoded_family_id_t family_id;
    family_data(const std::string a, const std::string b, const std::string c, encoded_family_id_t d)
    : pgf(a), plf(b), function(c), family_id(d) {};
    };
    typedef std::pair<std::string, std::string> family_key_t;
    
    /*
     * Maps string genus name from family file
     * to the taxon id used in PLF identifiers.
     */
    std::map<std::string, std::string> genus_map_;

#ifdef USE_TBB
    typedef tbb::concurrent_vector<encoded_family_id_t> family_counts_t;
    //typedef tbb::concurrent_unordered_map<encoded_family_id_t, unsigned long> family_counts_t;
    typedef tbb::concurrent_vector<encoded_id_t> id_set;
    typedef tbb::concurrent_unordered_map<encoded_kmer_t, id_set> map_type_t;
    typedef tbb::concurrent_unordered_map<encoded_kmer_t, family_counts_t> family_map_type_t;
    typedef tbb::concurrent_unordered_map<std::string, encoded_id_t> genome_to_id_map_t;
    typedef tbb::concurrent_unordered_map<encoded_id_t, std::string> id_to_genome_map_t;
    typedef tbb::concurrent_unordered_map<encoded_id_t, family_data> family_map_t;
    typedef tbb::concurrent_unordered_map<std::string, encoded_id_t> peg_to_id_map_t;
    typedef tbb::concurrent_unordered_map<encoded_id_t, std::string> id_to_peg_map_t;

    typedef tbb::concurrent_unordered_map<encoded_family_id_t, family_key_t> id_to_family_map_t;
    typedef tbb::concurrent_unordered_map<family_key_t, encoded_family_id_t> family_to_id_map_t;
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
    family_to_id_map_t family_to_id_;
    id_to_family_map_t id_to_family_;
    
    map_type_t kmer_to_id_;
    family_map_type_t kmer_to_family_id_;

    genome_to_id_map_t genome_to_id_;
    id_to_genome_map_t id_to_genome_;

    // family_mapping goes from encoded peg id to family data
    family_map_t family_mapping_;


    void load_genus_map(const std::string &genus_file);
    void load_families(const std::string &families_file);

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
    void add_fam_mapping(encoded_id_t enc, encoded_kmer_t kmer);

    std::string decode_id(encoded_id_t id);

    unsigned long kcount_;
    unsigned long next_genome_id_;

    void dump_sizes(std::ostream &os);
};


#endif
