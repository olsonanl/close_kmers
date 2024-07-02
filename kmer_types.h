#ifndef _kmer_types_h
#define _kmer_types_h

#include <climits>
#include <map>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <iostream>

/**
 * Function indexes are use to reference the entries in function.index
 * which represent function strings assigned to proteins.
 */
typedef unsigned short FunctionIndex;

/**
 * Value representing a missing or undefined function.
 */
const FunctionIndex UndefinedFunction = USHRT_MAX;

/**
 * OTU indexes are use to reference the entries in otu.index
 * which represent OTUs associated with kmers.
 */
typedef unsigned short OTUIndex;

/**
 * Value representing a missing or undefined function.
 */
const FunctionIndex UndefinedOTU = USHRT_MAX;

class KmerResult
{
public:
    std::string id;
    size_t len;
    
};


/*
 * For compatibility with current kmer builds,
 * we use the original kguts types here rather than
 * the shorter OTUIndex and FunctionIndex.
 */
typedef struct sig_kmer {
    unsigned long long  which_kmer;
    int  otu_index;
    unsigned short  avg_from_end;
    int  function_index;
    float function_wt;
} sig_kmer_t;

// NuDB data struct

struct KData {
    OTUIndex otu_index;
    unsigned short  avg_from_end;
    FunctionIndex function_index;
    float function_wt;
};

struct hit_in_sequence_t {
    sig_kmer_t hit;
    unsigned int offset;
    
    hit_in_sequence_t(const sig_kmer_t &h, unsigned int o) : hit(h), offset(o) {}
hit_in_sequence_t(const hit_in_sequence_t &h) : hit(h.hit), offset(h.offset) {}

hit_in_sequence_t(unsigned long long enc, const KData &k, unsigned int o) :
    hit {enc, k.otu_index, k.avg_from_end, k.function_index, k.function_wt}, offset{o} {}

};

class KmerHit
{
public:
    OTUIndex   oI;
    unsigned int from0_in_prot;      /* offset from start of protein sequence */
    unsigned short avg_off_from_end;  /* average offset from the end */
    FunctionIndex fI;
    float function_wt;
    unsigned long long encoded_kmer;
};
typedef KmerHit hit_t;

class KmerCall
{
public:
    unsigned int start;
    unsigned int end;
    int count;
    FunctionIndex function_index;
    float weighted_hits;

    KmerCall() : start(0), end(0), count(0), function_index(UndefinedFunction), weighted_hits(0.0) { }
    KmerCall(unsigned int s, unsigned int e, int c, FunctionIndex f, float w)  :
        start(s), end(e), count(c), function_index(f), weighted_hits(w) { }
KmerCall(KmerCall &&k) :
        start(k.start), end(k.end), count(k.count), function_index(k.function_index), weighted_hits(k.weighted_hits) { }
KmerCall(const KmerCall &k) :
        start(k.start), end(k.end), count(k.count), function_index(k.function_index), weighted_hits(k.weighted_hits) { }

};

class KmerOtuStats
{
public:
    std::string contig_id;
    int contig_len;

    std::map<int, int> otu_map;
    std::vector<std::pair<int, int>> otus_by_count;
		    

    template<class T>
	struct less_second : std::binary_function<T,T,bool>
    {
	inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
    };

    void write(FILE *fh)
    {
	fprintf(fh, "OTU-COUNTS\t%s[%d]", contig_id.c_str(), contig_len);
	for (auto it = otus_by_count.begin(); it != otus_by_count.end(); it++)
	{
	    fprintf(fh, "\t%d-%d",it->second, it->first);
	}
	fprintf(fh, "\n");
    }

    void finalize()
    {
	otus_by_count.insert(otus_by_count.begin(), otu_map.begin(), otu_map.end());
	std::sort(otus_by_count.begin(), otus_by_count.end(), less_second<std::pair< int, int> > ());
    }
};

class KmerParameters
{
public:
KmerParameters() :
    param_map_(
	{
	    { "order_constraint", order_constraint },
	    { "min_hits", min_hits },
	    { "min_weighted_hits", min_weighted_hits },
	    { "max_gap", max_gap }
	})

    {
	set_default_parameters();
    }
    
    void set_default_parameters()
    {
	order_constraint = 0;
	min_hits = 5;
	min_weighted_hits = 0;
	max_gap  = 200;
    }

    std::map<std::string, int &> param_map_;
    void set_parameters(const std::map<std::string, std::string> &params)
    {
	set_default_parameters();
	for (auto p: params)
	{
	    auto myp = param_map_.find(p.first);
	    if (myp != param_map_.end())
	    {
		try {
		    int val = std::stoi(p.second);
		    myp->second = val;
		} catch (const std::invalid_argument& ia)
		  {
		      std::cerr << "Warning: invalid integer '" << p.second << "' passed for parameter " << p.first << "\n";
		  }
	    }
	}
    }

    int   order_constraint;
    int   min_hits;
    int   min_weighted_hits;
    int   max_gap;
};

#endif /* _kmer_types_h */
