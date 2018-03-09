#ifndef _lookup_request_h
#define _lookup_request_h

#include "compute_request.h"
#include "krequest2.h"
#include "kmer.h"
#include "fasta_parser.h"
#include <memory>
#include "prot_seq.h"

class LookupRequest : public ComputeRequest
{
public:
    LookupRequest(std::shared_ptr<KmerRequest2> owner, std::shared_ptr<KmerPegMapping> mapping, bool family_mode, int content_length, bool chunked = false);
    ~LookupRequest() { std::cerr << "Destroy lookuprequest " << this << "\n"; }
    void run();

private:
    bool family_mode_;
    std::shared_ptr<KmerPegMapping> mapping_;
    int content_length_;
    bool chunked_;
    std::shared_ptr<KmerRequest2> owner_;
    FastaParser parser_;

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
    std::unordered_map<KmerPegMapping::encoded_id_t, sequence_accumulated_score_t> seq_score_;

    unsigned int kmer_hit_threshold_;
    bool header_written_;
    bool find_best_match_;
    bool allow_ambiguous_functions_;
    enum { best_by_accumulation, best_by_count, best_by_weighted_count } best_match_method_;
    std::string target_genus_;
    unsigned long target_genus_id_;
    bool find_reps_;
    
    typedef std::vector<ProteinSequence> work_list_t;
    std::shared_ptr<work_list_t> current_work_;

    int on_parsed_seq(const std::string &id, const std::string &seq);

    void on_hit(KmerGuts::hit_in_sequence_t kmer);
    void on_data(boost::system::error_code err, size_t bytes);

    void process_results();
};

#endif
