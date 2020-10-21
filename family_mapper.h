#ifndef _family_mapper_h
#define _family_mapper_h

/*
 * Process a protein to compute its family membership.
 */

#include <string>
#include "kguts.h"
#include "kmer.h"

#include <unordered_map>

class FamilyMapper
{
public:
    FamilyMapper(KmerGuts *kguts,
		 std::shared_ptr<KmerPegMapping> mapping);

    struct best_match_t {
	std::string gfam_id;
	float gfam_score;
	std::string lfam_id;
	float lfam_score;
	std::string function;
	float score;
	friend bool operator<(const best_match_t &l, const best_match_t &r) { return l.score < r.score; };
    };

    best_match_t find_best_family_match(const std::string &id, const std::string &seq);
    void find_all_matches(std::ostream &os, const std::string &id, const std::string &seq);

    struct sequence_accumulated_score_t {
	unsigned int hit_count;
	unsigned int hit_total;
	float weighted_total;

	inline void increment(KmerPegMapping::encoded_family_id_t val, float weight) {
	    hit_count++;
	    hit_total++;
	    weighted_total += weight;
	}
	    
	inline void increment(std::pair<KmerPegMapping::encoded_family_id_t, unsigned int> val, float weight) {
	    hit_count += val.second;
	    hit_total++;
	    weighted_total += weight;
	}
    } ;
private:

    void ingest_protein(const std::string &id, const std::string &seq);

    std::unordered_map<KmerPegMapping::encoded_id_t, sequence_accumulated_score_t> seq_score_;

    void on_hit(const KmerGuts::hit_in_sequence_t &kmer);

    bool find_best_match_;
    bool family_mode_;
    unsigned int kmer_hit_threshold_;
    bool allow_ambiguous_functions_;
    std::shared_ptr<KmerPegMapping> mapping_;
    unsigned long target_genus_id_;
    bool find_reps_;
    KmerGuts *kguts_;
    typedef std::vector<KmerCall> call_vector_t;
    std::shared_ptr<call_vector_t> calls_;
};

inline std::ostream &operator<<(std::ostream &os, const FamilyMapper::best_match_t &m)
{
    os << m.gfam_id << "\t" << m.gfam_score << "\t" << m.lfam_id << "\t" << m.lfam_score
       << "\t" << m.function << "\t" << m.score;
    return os;
}

#endif // _family_mapper_h
