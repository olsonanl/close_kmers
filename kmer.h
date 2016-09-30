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

    typedef unsigned long encoded_id_t;

    /* The kmer to peg mapping table. This is huge. */

    /* peg to peg-attributes mapping. */
    struct family_data {
	std::string pgf;
	std::string plf;
	std::string function;
	family_data(const std::string a, const std::string b, const std::string c)
	: pgf(a), plf(b), function(c) {};
    };

    /*
     * Maps string genus name from family file
     * to the taxon id used in PLF identifiers.
     */
    std::map<std::string, std::string> genus_map_;

#ifdef USE_TBB
    typedef tbb::concurrent_vector<encoded_id_t> id_set;
    typedef tbb::concurrent_unordered_map<encoded_id_t, id_set> map_type_t;
    typedef tbb::concurrent_unordered_map<std::string, encoded_id_t> genome_to_id_map_t;
    typedef tbb::concurrent_unordered_map<encoded_id_t, std::string> id_to_genome_map_t;
    typedef tbb::concurrent_unordered_map<encoded_id_t, family_data> family_map_t;
#else
    typedef std::vector<encoded_id_t> id_set;
    typedef std::unordered_map<encoded_id_t, id_set> map_type_t;
    typedef std::map<std::string, encoded_id_t> genome_to_id_map_t;
    typedef std::map<encoded_id_t, std::string> id_to_genome_map_t;
    typedef std::unordered_map<encoded_id_t, family_data> family_map_t;
#endif
    map_type_t kmer_to_id_;

    genome_to_id_map_t genome_to_id_;
    id_to_genome_map_t id_to_genome_;

    family_map_t family_mapping_;
    void load_genus_map(const std::string &genus_file);
    void load_families(const std::string &families_file);

    /*
     * We also maintain a mapping from an md5 value to the pegs with that
     * md5.
     */

    std::map<std::array<char, 16>, std::vector<encoded_id_t> > md5_to_peg_;

    encoded_id_t encode_id(const std::string &peg);
    encoded_id_t encode_id(const std::string &genome, const std::string &peg);

    void add_mapping(KmerPegMapping::encoded_id_t enc, unsigned long kmer);

    std::string decode_id(encoded_id_t id);

    unsigned long kcount_;
    unsigned long next_genome_id_;
};


#endif
