#include "kmer_nudb.h"

#include <algorithm>
#include <memory>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#include <boost/program_options.hpp>
#include "global.h"

KmerNudb::KmerNudb(const std::string &file_base)
{
}

KmerNudb::~KmerNudb() {
}


/*
void KmerNudb::process_set_of_hits(std::shared_ptr<std::vector<KmerCall>> calls,
				   std::shared_ptr<KmerOtuStats> otu_stats)
{
}
*/

void KmerNudb::process_aa_seq_hits(const std::string &id, const std::string &seq,
				   std::shared_ptr<std::vector<KmerCall>> calls,
				   std::shared_ptr<std::vector<hit_in_sequence_t>> hits,
				   std::shared_ptr<KmerOtuStats> otu_stats)
{
    auto cb = [this, hits](hit_in_sequence_t k) { hits->push_back(k); };
    process_aa_seq(id, seq, calls, cb, otu_stats);
}

void KmerNudb::process_aa_seq(const std::string &idstr, const std::string &seqstr,
			      std::shared_ptr<std::vector<KmerCall>> calls,
			      std::function<void(const hit_in_sequence_t &)> hit_cb,
			      std::shared_ptr<KmerOtuStats> otu_stats)
{
//    gather_hits2(idstr, seqstr, calls, hit_cb, otu_stats);
    if (otu_stats)
	otu_stats->finalize();
}

void KmerNudb::process_seq(const char *id,const char *data,
			   std::shared_ptr<std::vector<KmerCall>> calls,
			   std::function<void(const hit_in_sequence_t &)> hit_cb,
			   std::shared_ptr<KmerOtuStats> otu_stats)

{
}

std::string KmerNudb::format_call(const KmerCall &c)
{
    std::ostringstream oss;
    oss << "CALL\t" << c.start << "\t" << c.end << "\t" << c.count;
//    oss << "\t" << c.function_index << "\t" << function_at_index(c.function_index);
    oss << "\t" << c.weighted_hits << "\n";

    return oss.str();
}

std::string KmerNudb::format_hit(const hit_in_sequence_t &h)
{
    std::ostringstream oss;

    char dc[KMER_SIZE + 1];
    encoder_.decoded_kmer(h.hit.which_kmer, dc);

//    oss << "HIT\t" << h.offset << "\t" << dc << "\t" << h.hit.avg_from_end << "\t" << function_at_index(h.hit.function_index) << "\t" << h.hit.function_wt << "\t" << h.hit.otu_index << "\n";
    
    return oss.str();
}

std::string KmerNudb::format_otu_stats(const std::string &id, size_t size, KmerOtuStats &otu_stats)
{
    std::ostringstream oss;
    oss << "OTU-COUNTS\t" << id << "[" << size << "]";

    int max_to_print = 5;
    for (auto x = otu_stats.otus_by_count.begin();
	 x != otu_stats.otus_by_count.end() && max_to_print > 0; x++, max_to_print--)
    {
	oss << "\t" << x->second << "-" << x->first;
    }
    oss << "\n";
    return oss.str();
}

template<class T>
struct less_second : std::binary_function<T,T,bool>
{
    inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
};  

struct FuncScore
{
    int count;
    float weighted;
    FuncScore() : count(0), weighted(0.0) {}
    FuncScore(int c, float w) : count(c), weighted(w) {}
    FuncScore(const FuncScore &f) : count(f.count), weighted(f.weighted) {}
    FuncScore& operator=(const FuncScore& other)
        {
	    // check for self-assignment
	    if(&other == this)
		return *this;
	    count = other.count;
	    weighted = other.weighted;
	    return *this;
	}
};

/*
 * Find the best call from this set of calls.
 *
 * This code replicates the amino acid version of the SEED pipeline
 * km_process_hits_to_regions | km_pick_best_hit_in_peg
 */
void KmerNudb::find_best_call(std::vector<KmerCall> &calls, FunctionIndex &function_index, std::string &function, float &score, float &weighted_score, float &score_offset)
{
    function_index = UndefinedFunction;
    function = "";
    score = 0.0;
    weighted_score = 0.0;

    if (calls.size() == 0)
    {
	return;
    }
    
    /*
     * First merge adjacent hits that have the same function.
     */
    std::vector<KmerCall> collapsed;

    auto comp = calls.begin();

    while (comp != calls.end())
    {
	collapsed.push_back(*comp);
	comp++;
	KmerCall &cur = collapsed.back();

	while (comp != calls.end() && cur.function_index == comp->function_index)
	{
	    cur.end = comp->end;
	    cur.count += comp->count;
	    cur.weighted_hits += comp->weighted_hits;
	    comp++;
	}
    }
#if 0
    std::cout << "after collapse:\n";
    for (auto iter = collapsed.begin(); iter != collapsed.end(); iter++)
    {
	std::cout << format_call(*iter);
    }
#endif

    /*
     *
     * Merge hits when we have a case with
     * +------+--+-------+--+-----------+
     * |  F1  |  |   F2  |  |   F1      |
     * +------+--+-------+--+-----------+
     *
     * where the score for F2 is below 5 and the combined scores for
     * the two F1 hits is 10 or more.
     *
     * If that is the case we discard the F2 hit and combine
     * the F1 hits.
     */

    std::vector<KmerCall> merged;

    int merge_interior_thresh = 5;
    int merge_exterior_thresh = 10;

    comp = collapsed.begin();
    while (comp != collapsed.end())
    {
	merged.push_back(*comp);
	comp++;
	auto comp2 = comp + 1;
	KmerCall &cur = merged.back();
	while (comp != collapsed.end() && comp2 != collapsed.end() &&
	       cur.function_index == comp2->function_index &&
	       comp->count < merge_interior_thresh &&
	       (cur.count + comp2->count) >= merge_exterior_thresh)
	{
	    cur.end = comp2->end;
	    cur.count += comp2->count;
	    cur.weighted_hits += comp2->weighted_hits;
	    comp += 2;
	    comp2 = comp + 1;
	}
    }
   
#if 0
    std::cerr << "after merge:\n";
    for (auto iter = merged.begin(); iter != merged.end(); iter++)
    {
	std::cerr << format_call(*iter);
    }
#endif

    /*
     * Determine best call.
     *
     * In the current perl kmer_search (km_pick_best_hit_in_peg) code we just take the best
     * function in terms of weighted score. However, that allows tied scores to be called
     * arbitrarily, and close-to-tied scores to be settled based on insignificant differences.
     *
     * In the original kmerv1 code, we required a score threshold difference (typically 5 hits)
     * between the best function and next best function. We resurrect that here.
     *
     */

    typedef std::map<int, FuncScore> map_t;

    map_t by_func;

    for (auto c: merged)
    {
	auto it = by_func.find(c.function_index);
	if (it == by_func.end())
	{
	    by_func.insert(std::make_pair(c.function_index, FuncScore(c.count, c.weighted_hits)));
	}
	else
	{
	    it->second.count += c.count;
	    it->second.weighted += c.weighted_hits;
	}
    }

    //typedef map_t::value_type ent_t;
    typedef std::pair<FunctionIndex, FuncScore> ent_t;

    std::vector<ent_t> vec;
    for (auto it = by_func.begin(); it != by_func.end(); it++)
	vec.push_back(*it);

//    std::cerr << "vec len " << vec.size() << "\n";
    if (vec.size() > 1)
    {
	std::partial_sort(vec.begin(), vec.begin() +  2, vec.end(), 
			  [](const ent_t& s1, const ent_t& s2) {
			      return s1.second.weighted > s2.second.weighted; });
    }
    
#if 0
    for (auto x: vec)
    {
	std::cerr << x.first << " " << x.second.count << " " << x.second.weighted << " ";
	std::cerr << function_at_index(x.first) << "\n";
    }
#endif
    
    if (vec.size() == 1)
	score_offset = (float) vec[0].second.count;
    else
	score_offset = (float) (vec[0].second.count - vec[1].second.count);
    
    // std::cerr << "Offset=" << score_offset << "\n";

    if (score_offset >= 5.0)
    {
	auto best = vec[0];
	function_index = best.first;
//	function = function_at_index(function_index);
	score = (float) best.second.count;
	weighted_score = best.second.weighted;
    }
    else
    {
	function_index = UndefinedFunction;
	function = "";
	score = 0.0;

	/*
	 * Try to compute a fallback function naming the two best hits if there are two hits within the
	 * threshold but greater than the next hit.
	 */
	if (vec.size() >= 2)
	{
	    std::string f1 = function_at_index(vec[0].first);
	    std::string f2 = function_at_index(vec[1].first);
	    if (f2 > f1)
		std::swap(f1, f2);

	    if (vec.size() == 2)
	    {
		function = f1 + " ?? " + f2;
		score = (float) vec[0].second.count;
	    }
	    else if (vec.size() > 2)
	    {
		float pair_offset = (float) (vec[1].second.count - vec[2].second.count);
		if (pair_offset > 5.0)
		{
		    function = f1 + " ?? " + f2;
		    score = (float) vec[0].second.count;
		    score_offset = pair_offset;
		    weighted_score = vec[0].second.weighted;
		}
	    }
	}
    }
}
