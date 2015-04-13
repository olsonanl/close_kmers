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

    typedef unsigned long encoded_id_t;
    typedef std::vector<encoded_id_t> id_set;

    std::string data_dir_;

    std::map<std::string, encoded_id_t> genome_to_id_;
    std::map<encoded_id_t, std::string> id_to_genome_;
	    
    std::unordered_map<encoded_id_t, id_set> kmer_to_id_;

    std::string decode_id(encoded_id_t id);
};

