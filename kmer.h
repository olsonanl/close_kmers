#ifndef _KMER_H
#define _KMER_H

/*
 * Kmer to peg mapping database
 */

#include <unordered_map>
#include <vector>


class KmerPegMapping
{
public:
    KmerPegMapping(const std::string &data_dir);

    void load_mapping_file(const std::string &mapping_file);
    void load_compact_mapping_file(const std::string &mapping_file);
    void load_compact_mapping_file(std::istream &kfile);

    typedef unsigned long encoded_id_t;
    typedef std::vector<encoded_id_t> id_set;

    std::string data_dir_;

    std::map<std::string, encoded_id_t> genome_to_id_;
    std::map<encoded_id_t, std::string> id_to_genome_;

    /* The kmer to peg mapping table. This is huge. */
       
    std::unordered_map<encoded_id_t, id_set> kmer_to_id_;

    /* peg to peg-attributes mapping. */
    std::unordered_map<encoded_id_t, std::pair<std::string, std::string>> family_mapping_;
    void load_families(const std::string &families_file);

    encoded_id_t encode_id(const std::string peg);
    encoded_id_t encode_id(const std::string &genome, const std::string &peg);

    void add_mapping(KmerPegMapping::encoded_id_t enc, unsigned long kmer);

    std::string decode_id(encoded_id_t id);

    unsigned long kcount_;
    unsigned long next_genome_id_;
};


#endif
