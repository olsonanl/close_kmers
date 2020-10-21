#ifndef _kmer_nudb_h
#define _kmer_nudb_h

#define _FILE_OFFSET_BITS 64
#include <stddef.h>
#include <netinet/in.h>
#include <string>
#include <memory>
#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include <ostream>

#include "kmer_params.h"
#include "kmer_encoder.h"
#include "kmer_types.h"

class KmerNudb
{
public:

    KmerEncoder encoder_;
    KmerParameters parameters_;

    KmerNudb(const std::string &file_base);
    ~KmerNudb(); 

    void gather_hits(const char *pseq,
		     unsigned char *pIseq,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
	
    void process_seq(const char *id,const char *data,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
	
    void process_aa_seq(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			std::function<void(const hit_in_sequence_t &)> hit_cb,
			std::shared_ptr<KmerOtuStats> otu_stats);

    void process_aa_seq_hits(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			std::shared_ptr<std::vector<hit_in_sequence_t>> hits,
			std::shared_ptr<KmerOtuStats> otu_stats);

    std::string format_call(const KmerCall &c);
    std::string format_hit(const hit_in_sequence_t &h);
    std::string format_otu_stats(const std::string &id, size_t size, KmerOtuStats &otu_stats);

    void find_best_call(std::vector<KmerCall> &calls, FunctionIndex &function_index,
			std::string &function, float &score, float &weighted_score, float &score_offset);
};

inline std::ostream &operator<<(std::ostream &os, const KmerCall &c)
{
    os << "KmerCall(" << c.start << "-" << c.end << ": " << c.count << ", " << c.function_index << ", " << c.weighted_hits << ")";
    return os;
}

#endif /* _kmer_guts_h */
