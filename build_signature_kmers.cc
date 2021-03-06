/*
 * Thoughts.
 *
 * We need function definitions for each sequence.
 *
 * Loading from a SEED, it makes sense to have separate files of fasta
 * protein data and annotation data.
 *
 * Loading from independent fasta, it makes sense to put the annotation
 * in the def line.
 *
 * Initial plan is to allow loading of annotation definitions id <tab> func,
 * then have annotations in a def line take precedence.
 */
 

#include <cmath>
#include <array>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <iostream>
#include <fstream>
#include <functional>
#include <memory>
#include <experimental/optional>
#include <unistd.h>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread/tss.hpp>
#include <boost/thread/thread.hpp>

#include "seed_utils.h"
#include "operators.h"
#include "fasta_parser.h"
#include "kguts.h"

#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/compat/thread"
#include "tbb/concurrent_hash_map.h"

using namespace seed_utils;

#define DEFINE_GLOBALS
#include "global.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

const int K = KMER_SIZE;
const int MaxSequencesPerFile = 100000;

typedef std::array<char, K> Kmer;
namespace std {
template<class T, size_t N>
struct hash<std::array<T, N>> {
    auto operator() (const array<T, N>& key) const {
	hash<T> hasher;
	size_t result = 0;
	for(size_t i = 0; i < N; ++i) {
	    result = result * 31 + hasher(key[i]); // ??
	}
	return result;
    }
};
}

struct KmerAttributes
{
    int func_index;
    int otu_index;
    int offset;
    unsigned seq_id;
};

namespace tbb {
    // Test whether tbb_hash_compare can be partially specialized as stated in Reference manual.
    template<> struct tbb_hash_compare<Kmer> {
	size_t hash( Kmer ) const {return 0;}
	bool equal( Kmer /*x*/, Kmer /*y*/ ) {return true;}
    };
}

struct tbb_hash {
    tbb_hash() {}
    size_t operator()(const Kmer& k) const {
        size_t h = 0;
	for (auto s: k)
	    h = (h*17)^s;
	return h;
    }
};

template<class T, size_t N>
struct MyHashCompare {
    static size_t hash( const std::array<T,N> &x ) {
        size_t h = 0;
	for (auto s: x)
	    h = (h*17)^s;
	return h;
    }
    //! True if strings are equal
    static bool equal( const std::array<T,N>&  x, const std::array<T,N> & y ) {
	return x==y;
    }
};

//typedef std::unordered_multimap<Kmer, KmerAttributes> KmerAttributeMap;
typedef tbb::concurrent_unordered_multimap<Kmer, KmerAttributes, tbb_hash> KmerAttributeMap;

inline std::ostream &operator<<(std::ostream &os, const Kmer &k)
{
    std::copy(k.begin(), k.end(),
	      std::ostream_iterator<char>(os, ""));
    return os;
}

struct KmerStatistics
{
    tbb::atomic<int> distinct_signatures = 0;
    tbb::concurrent_unordered_map<int, int> distinct_functions;
    tbb::concurrent_unordered_map<int, int> seqs_with_func;
    tbb::concurrent_unordered_set<int> seqs_with_a_signature;
};

//
// Whack some globals in here. We're maintaining enough state this
// needs to be wrapped up in an object.
//

KmerStatistics kmer_stats;

std::ofstream rejected_stream;
std::ofstream kept_stream;
std::ofstream kept_function_stream;

/*
 * Object to keep state for kmers we are keeping. We hang onto
 * some statistics in order to compute weights later on.
 */
struct KeptKmer
{
    Kmer kmer;
    int median_offset;		// Median offset from the end of the protein
    int function_index;
    int otu_index;
    int seqs_containing_sig;	// Count of sequences containing this kmer
    int seqs_containing_function; // Count of sequences with the kmer that have the function
    float weight;
};

tbb::concurrent_vector<KeptKmer> kept_kmers;

struct KmerSet
{
    KmerSet() : count(0) {}
    void reset() {
	count = 0;
	func_count.clear();
	set.clear();
    }
    // ~KmerSet() { std::cerr << "destroy " << this << "\n"; }
    Kmer kmer;
    std::map<int, int> func_count;
    int count;
    std::vector<KmerAttributes> set;
};

inline std::ostream &operator<<(std::ostream &os, const KmerSet &ks)
{
    os << "KmerSet: kmer=" << ks.kmer << " count=" << ks.count << "\n";
    return os;
}

void process_set(KmerSet &set);
typedef std::vector<KmerSet> KmerSetList;
typedef std::shared_ptr<KmerSetList> KmerSetListPtr;

struct KmerProcessor
{
    KmerProcessor(int n) : n_threads(n) {}
    boost::thread_group thread_pool;
    tbb::concurrent_bounded_queue<KmerSetListPtr> queue;
    int n_threads;
    void start() {
	for (int i = 0; i < n_threads; i++)
	{
	    std::cerr << "starting " << i << "\n";
	    thread_pool.create_thread([this, i]() {
		    thread_main(i);
		});
	}
    }
    void stop() {
	for (int i = 0; i < n_threads; i++)
	{
	    std::cerr << "stopping " << i << "\n";
	    KmerSetListPtr work = std::make_shared<KmerSetList>();
	    enqueue_work(work);
	}
	std::cerr << "Awaiting threads\n";
	thread_pool.join_all();
    }
	
    void enqueue_work(KmerSetListPtr work) {
	queue.push(work);
    }
    void thread_main(int i) {
	std::cerr << "running " << i << "\n";
	int sp = 0;
	while (1)
	{
	    KmerSetListPtr work;
	    queue.pop(work);
	    if (work->size() == 0)
	    {
		std::cerr << "shutting down " << i << "\n";
		break;
	    }

	    for (auto entry: *work)
	    {
		// std::cerr << i << " process " << entry << "\n";
		sp++;
		process_set(entry);
	    }
	}
	std::cerr << "thread " << i << " processed " << sp << " sets\n";
    }
};
KmerProcessor *g_kmer_processor;

/*
 * Manage the ID to function mapping.
 *
 * FunctionMap also maintains the database of function => genome mappings
 * that is used to determine whether a given function has enough evidence
 * available to be used to create signature kmers.
 *
 * Functions occurring in genomes are determined by the presence of a
 * protein with the given function occurring in the fasta file for
 * a genome.
 *
 * The genome for a fasta file is determined using either the identifier
 * from the first sequence in the file (if it is a fig| identifier), or
 * the genome name in square brackets at the end of the definition line
 * in the case of fasta files loaded in standard genbank format from
 * outside the SEED environment.
 *
 * Members include:
 *
 *  function_genome_map_: maps from function string to the list of genomes in which it occurs.
 			  Initialized from scanning the protein fasta files.
 *
 *  id_function_map_: maps from protein ID to assigned function.
 *		      Initialized from the function def files, overriden by fasta files
 * 		
 */
class FunctionMap
{
public:
    FunctionMap() {
    }

    void load_id_assignments(const fs::path &file) {
	fs::ifstream ifstr(file);
	std::string line;
	int lineno = 0;
	while (std::getline(ifstr, line))
	{
	    lineno++;
	    size_t s = line.find('\t');
	    if (s == std::string::npos)
	    {
		std::cerr << "bad line " << lineno << " in file " << file << "\n";
		continue;
	    }
	    size_t s2 = line.find('\t', s+1);
	    std::string id = line.substr(0, s);
	    std::string func;
	    if (s2 == std::string::npos)
		func = line.substr(s+1);
	    else
		func = line.substr(s+1, s2 - s - 1);

	    std::string stripped = strip_func_comment(func);
	    id_function_map_[id] = stripped;
	    // std::cerr << "Load '" << id << "' as '" << func << "' stripped='" << stripped << "'\n";
	}
    }

    /*
     * Load assignments and sequence visibility data from a fasta file.
     *
     * If the def line has an assignment on it, update the assignment map.
     *
     * We also update the function_genome map to record which functions appear
     * in which genomes.
     *
     * We make the simplifying assumption that there is a 1:1 mapping
     * between the proteins in a genome and the contents of a given fasta file.
     */
    void load_fasta_file(const fs::path &file, bool keep_function_flag) {

	const boost::regex genome_regex("\\s+(.*)\\s+\\[([^]]+)\\]$");
	const boost::regex figid_regex("fig\\|(\\d+\\.\\d+)");
	const boost::regex genome_id_regex("\\d+\\.\\d+");
	
	fs::ifstream ifstr(file);

	FastaParser parser;

	std::string genome;

	parser.set_def_callback([this, &genome, &genome_regex, &figid_regex, &genome_id_regex, &file, keep_function_flag]
				(const std::string &id, const std::string &def, const std::string &seq) {
		if (id.empty())
		    return 0;

		boost::smatch match;

		//
		// Need to always parse for [genome]
		//

		std::string func;
		if (!def.empty())
		{
		    size_t x = def.find_first_not_of(" \t");
		    func = def.substr(x);
		}
		std::string genome_loc;
		if (boost::regex_match(def, match, genome_regex))
		{
		    func = strip_func_comment(match[1]);
		    genome_loc = match[2];
		}
		
		// Determine genome from first sequence.
		if (genome.empty())
		{
		    if (def.empty())
		    {
			if (boost::regex_search(id, match, figid_regex))
			{
			    genome = match[1];
			}
		    }
		    else
		    {
			if (!genome_loc.empty())
			{
			    genome = genome_loc;
			}
		    }
		}
		if (genome.empty())
		{
		    // default it to the file, just to have a value
		    genome = file.leaf().string();
		    
		    if (!boost::regex_match(genome, genome_id_regex))
		    {
			std::cerr << "cannot determine genome from file " << file << "\n";
		    }
		}

		/*
		 *
		 * If the current-function map has a value, use that. We assume explicit setting
		 * of functions overrides what is in the fasta.
		 *
		 * If we're assigning a function, update the id to function
		 * map.
		 *
		 * Then look up the function for this id and add to the function_genome map.
		 */
		std::string cur_func = id_function_map_[id];
		if (cur_func.empty())
		{
		    if (!func.empty())
		    {
			id_function_map_[id] = func;
		    }
		}
		else
		{
		    func = cur_func;
		}
		
		if (func.empty())
		{
		    // std::cerr << "No function found for " << id << "\n";
		}
		else
		{
		    function_genome_map_[func].insert(genome);
		    if (keep_function_flag)
		    {
			// std::cerr << "Keeping function " << func << "\n";
			good_functions_.insert(func);
		    }
		}
			
		return 0;
	    });
	parser.parse(ifstr);
	parser.parse_complete();
    }

    /*
     * Process the contents of function_genome_map to determine which
     * function for which we should build signatures.
     *
     * Criteria - one of the following must be met:
     *
     *    Genome count for a function >= min_reps_required
     *
     *    One of the roles in the function is in the good_roles set.
     *
     *    The function is in the good_functions set.
     *
     * For each function that we keep we will assign an identifier in the
     * function_index_map which is then used as the means to describe
     * which functions are to have signatures created.
     */
    void process_kept_functions(int min_reps_required) {
	std::set<std::string> kept;
	for (auto entry: function_genome_map_)
	{
	    auto function = entry.first;
	    int n_genomes = entry.second.size();
	    kept_function_stream << function << ": " << n_genomes << " genomes\n";
	    bool ok = false;

	    if (n_genomes >= min_reps_required)
	    {
		kept_function_stream << "Keeping " << function << ": enough genomes\n";
		ok = true;
	    }
	    else if (good_functions_.find(function) != good_functions_.end())
	    {
		kept_function_stream << "Keeping " << function << ": in good functions list\n";
		ok = true;
	    }
	    else
	    {
		std::vector<std::string> roles = roles_of_function(function);
		kept_function_stream << "Role check " << function << ":\n";
		for (auto role: roles)
		{
		    if (good_roles_.find(role) != good_roles_.end())
		    {
			kept_function_stream << "  Keeping " << function << ": " << role << " in good roles list\n";
			ok = true;
			break;
		    }
		    else
		    {
			kept_function_stream << "  " << function << ": " << role << " not in list\n";
		    }
		}

		if (!ok)
		{
		    kept_function_stream << "Reject " << function << "\n";
		}
	    }
	    if (ok)
		kept.insert(function);
	}

	/*
	 * Assign sequential function IDs.
	 */
	int next = 0;
	for (auto f: kept)
	{
	    int id = next++;
	    function_index_map_[f] = id;
	}
	std::cout << "kept " << next << " functions\n";
    }

    void dump() {
	for (auto x: function_genome_map_)
	{
	    std::cout << x.first << ":";
	    for (auto y: x.second)
		std::cout << " " << y;
	    std::cout << "\n";
	}
	for (auto x: id_function_map_)
	{
	    std::cout << x.first << " '" << x.second << "'\n";
	}
    }

    std::string lookup_function(const std::string &id) {
	auto it = id_function_map_.find(id);
	if (it == id_function_map_.end())
	    return "";
	else
	    return it->second;
    }
    int lookup_index(const std::string &func) {
	auto it = function_index_map_.find(func);
	if (it == function_index_map_.end())
	    return -1;
	else
	    return it->second;
    }

    void write_function_index(const fs::path &dir) {
	fs::ofstream of(dir / "function.index");

	std::map<int, std::string> by_index;
	for (auto ent: function_index_map_)
	    by_index.insert(std::make_pair(ent.second, ent.first));
	for (auto ent: by_index)
	{
	    of << ent.first << "\t" << ent.second << "\n";
	}
    }

    // Create a function-array for compatibility with kguts
    char **create_kg_function_array(int *size) {
	int n = function_index_map_.size();
	char **ret = (char **) malloc(sizeof(char *) * n);
	
	for (auto ent: function_index_map_)
	{
	    if (ent.second >= n)
	    {
		std::cerr << "Bad ent=" << ent.second << " for n=" << n << "\n";
		exit(1);
	    }
	    ret[ent.second] = const_cast<char *>(ent.first.c_str());
	}
	//for (int i = 0; i < n;  i++)
	//    std::cerr << i << ": '" << ret[i] << "'\n";
	*size = n;
	return ret;
    }
	

    void add_good_roles(const std::vector<std::string> &r) {
	std::copy(r.begin(), r.end(), std::inserter(good_roles_, good_roles_.end()));
    }
    
    void add_good_functions(const std::vector<std::string> &r) {
	std::copy(r.begin(), r.end(), std::inserter(good_functions_, good_functions_.end()));
    }
    
private:
    std::map<std::string, std::set<std::string> > function_genome_map_;
    std::map<std::string, std::string> id_function_map_;
    std::map<std::string, int> function_index_map_;

    std::set<std::string> good_roles_;
    std::set<std::string> good_functions_;
};

std::set<char> ok_prot = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
			   'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y'};

void load_sequence(FunctionMap &fm, KmerAttributeMap &m, unsigned &next_sequence_id,
		   const std::string &id, const std::string &def, const std::string &seq)
{
    if (id.empty())
	return;

    /*
     * CHECK THIS: we should be using using functions from our database.
    std::string func = def;
    if (func.empty())
    {
	func = fm.lookup_function(id);
    }
     */

    std::string func = fm.lookup_function(id);
    
    /*
     * Empty means empty (and perhaps deleted feature).
     */
    
    if (func.empty())
    {
	return;
    }
    
    unsigned seq_id = next_sequence_id++;

    int function_index = fm.lookup_index(func);
    // std::cout << "Load seq " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
    if (function_index < 0)
	return;

    kmer_stats.seqs_with_func[function_index]++;

    for (auto it = seq.begin(); it < seq.end() - K + 1; it++)
    {
        int n = std::distance(it, seq.end());
	Kmer kmer;
	Kmer::iterator kiter = kmer.begin();
	bool ok = true;
	std::copy_n(it, K, kmer.begin());
	for (auto x: kmer)
	{
	    if (ok_prot.find(x) == ok_prot.end())
	    {
		// std::cerr << "not ok prot "  << x << "\n";
		ok = false;
		break;
	    }
	}
	if (ok)
	{
	    m.insert({kmer, { function_index, -1, n, seq_id }});
	    // std::cout << kmer << " " << n << "\n";
	}
	else
	{
	    // std::cout << "reject " << kmer << " " << n << "\n";
	}
    }
}


/*
 * Load sequence data.
 *
 * If a function index is not defined, this is not a function we wish
 * to process so skip this sequence.
 */
void load_fasta(FunctionMap &fm, KmerAttributeMap &m, unsigned file_number, const fs::path &file)
{
    fs::ifstream ifstr(file);

    FastaParser parser;
    
    unsigned next_sequence_id = file_number * MaxSequencesPerFile;

    parser.set_def_callback([&fm, &m, &next_sequence_id](const std::string &id, const std::string &def, const std::string &seq) {
	    load_sequence(fm, m, next_sequence_id, id, def, seq);
	    return 0;
	});
    parser.parse(ifstr);
    parser.parse_complete();
}

void enqueue_set(KmerSetListPtr set_list)
{
    g_kmer_processor->enqueue_work(set_list);
}

void process_set(KmerSet &set)
//void process_set(Kmer &kmer, std::map<int, int> &func_count, int count, std::vector<KmerAttributes> &set)
{
    auto elt = std::max_element(set.func_count.begin(), set.func_count.end(),
				[](auto a, auto b) { return a.second < b.second; });

    int best_func = elt->first;
    int best_count = elt->second;

    float thresh = float(set.count) * 0.8;

    /*
    kept_stream << "Process set for " << kmer << " best=" << best_func << " best_count=" << best_count << " count=" << count << " thresh=" << thresh << "\n";
    for (auto x: set)
    {
	kept_stream << "  " << x.func_index << " " << x.seq_id << "\n";
    }
    */

    if ((float) best_count < thresh)
    {
	// std::cout << "discard for " << best_count << " < " << thresh << "\n";
	// rejected_stream << kmer <<  " best_count=" << best_count << " thresh=" << thresh << "\n";
	return;
    }

    int seqs_containing_func = 0;
    std::vector<int> offsets;

    for (auto item: set.set)
    {
	if (item.func_index == best_func)
	{
	    seqs_containing_func++;
	}
	offsets.push_back(item.offset);
	kmer_stats.seqs_with_a_signature.insert(item.seq_id);
    }
    std::sort(offsets.begin(), offsets.end());
    int median_offset = offsets[offsets.size() / 2];
    // std::cout << seqs_containing_func << " " << median_offset<< "\n";

    kmer_stats.distinct_signatures++;
    kmer_stats.distinct_functions[best_func]++;

    kept_kmers.push_back(KeptKmer { set.kmer, median_offset, best_func, -1,
				       (int) set.set.size(), seqs_containing_func });
}


void process_kmers(KmerAttributeMap &m)
{
    //std::map<int, int> func_count;
    //std::vector<KmerAttributes> set;
    //int count = 0;
    Kmer cur { 0 } ;
    
    KmerSetListPtr cur_set_list = std::make_shared<KmerSetList>();
    cur_set_list->emplace_back(KmerSet());
//    KmerSet cur_set;

    int sets_pushed = 0;
	
    for (auto ent: m)
    {
	const Kmer &kmer = ent.first;
	KmerAttributes &attr = ent.second;

	if (kmer != cur)
	{
	    if (cur_set_list->back().count > 0)
	    {
		sets_pushed++;
		if (cur_set_list->size() > 200000)
		{
		    /*
		    std::cerr << "push setlist size " << cur_set_list->size() << "\n";
		    for (auto x: *cur_set_list)
		    {
			std::cerr << "  " << x;
		    }
		    */
		    enqueue_set(cur_set_list);
		    cur_set_list = std::make_shared<KmerSetList>();
		}
		cur_set_list->emplace_back(KmerSet());
	    }
	    //cur_set.reset();
	    cur_set_list->back().kmer = kmer;
	    cur = kmer;
	}
	KmerSet &s = cur_set_list->back();
	s.func_count[attr.func_index]++;
	s.count++;
	s.set.emplace_back(attr);
    }
    sets_pushed++;
//    cur_set_list->emplace_back(cur_set);
    enqueue_set(cur_set_list);
    std::cout << "done, pushed " << sets_pushed << " sets\n";
}

void par_process_kmers(KmerAttributeMap &m)
{
    tbb::parallel_for(m.range(), [](auto r) {
	    KmerSet cur_set;
	    Kmer cur { 0 };
	    for (auto ent = r.begin(); ent != r.end(); ent++)
	    {
		const Kmer &kmer = ent->first;
		KmerAttributes &attr = ent->second;

		if (kmer != cur)
		{
		    if (cur_set.count > 0)
			process_set(cur_set);
		    
		    cur_set.reset();
		    cur_set.kmer = kmer;
		    cur = kmer;
		}
		cur_set.func_count[attr.func_index]++;
		cur_set.count++;
		cur_set.set.emplace_back(attr);
	    }
	    process_set(cur_set);
	});
}

void process_kmer_block(KmerAttributeMap &m)
{
    Kmer cur { 0 } ;
    
    KmerSetListPtr cur_set_list = std::make_shared<KmerSetList>();
    cur_set_list->emplace_back(KmerSet());
//    KmerSet cur_set;

    int sets_pushed = 0;
	
    for (auto ent: m)
    {
	const Kmer &kmer = ent.first;
	KmerAttributes &attr = ent.second;

	if (kmer != cur)
	{
	    if (cur_set_list->back().count > 0)
	    {
		sets_pushed++;
		if (cur_set_list->size() > 200000)
		{
		    /*
		    std::cerr << "push setlist size " << cur_set_list->size() << "\n";
		    for (auto x: *cur_set_list)
		    {
			std::cerr << "  " << x;
		    }
		    */
		    enqueue_set(cur_set_list);
		    cur_set_list = std::make_shared<KmerSetList>();
		}
		cur_set_list->emplace_back(KmerSet());
	    }
	    //cur_set.reset();
	    cur_set_list->back().kmer = kmer;
	    cur = kmer;
	}
	KmerSet &s = cur_set_list->back();
	s.func_count[attr.func_index]++;
	s.count++;
	s.set.emplace_back(attr);
    }
    sets_pushed++;
//    cur_set_list->emplace_back(cur_set);
    enqueue_set(cur_set_list);
    std::cout << "done, pushed " << sets_pushed << " sets\n";
}

void compute_weight_of_signature(KeptKmer &kk)
{
    float NSF = kmer_stats.seqs_with_a_signature.size();
    float KS  = kmer_stats.distinct_signatures;
    // float KF  = kmer_stats.seqs_with_func.size();
    float NSi = kk.seqs_containing_sig;
    float NFj = kmer_stats.seqs_with_func[kk.function_index];
    float NSiFj = kk.seqs_containing_function;

    kk.weight = log((NSiFj + 1.0) / (NSi - NSiFj + 1.0)) +
	log((NSF - NFj + KS) / (NFj + KS));

}

void write_function_index(const fs::path &dir, FunctionMap &fm)
{
    fm.write_function_index(dir);
}

KmerGuts *write_hashtable(const fs::path &dir, tbb::concurrent_vector<KeptKmer> &kmers)
{
    std::vector<unsigned long> primes {3769,6337,12791,24571,51043,101533,206933,400187,
	    821999,2000003,4000037,8000009,16000057,32000011,
	    64000031,128000003,248000009,508000037,1073741824,
	    1400303159,2147483648,1190492993,3559786523,6461346257};

    //
    // find a hash table size appropriate for the number of kmers we have.
    //
    unsigned long hashtable_size = 0;
    for (auto p: primes)
    {
	if (p > 3 * kmers.size())
	{
	    hashtable_size = p;
	    break;
	}
    }
    if (hashtable_size == 0)
    {
	std::cerr << "Cannot find prime for " << kmers.size() << "\n";
	exit(1);
    }
    std::cerr << "Using hashtable size " << hashtable_size << " for " << kmers.size() << "\n";
	
    KmerGuts *kguts = new KmerGuts(dir.string(), hashtable_size);

    for (auto k: kmers)
    {
	kguts->insert_kmer(kguts->encoder_.encoded_aa_kmer(k.kmer),
			  k.function_index, k.otu_index, k.median_offset, k.weight);
    }
    if (!dir.empty())
    {
	kguts->save_kmer_hash_table((dir / "kmer.table.mem_map").string());
    }
    return kguts;
}

struct RecallData
{
    std::string new_function;
    std::string old_function;
    float score;
    float weighted_score;
    float score_offset;
};

std::experimental::optional<RecallData> recall_sequence(FunctionMap &fm, KmerGuts *kguts,
//		    fs::ofstream &calls_stream, fs::ofstream &new_stream,
		    const std::string &id, const std::string &seq)
{

    if (id.empty())
	return {};
    typedef std::vector<KmerCall> call_vector_t;
    std::shared_ptr<call_vector_t> calls_ = std::make_shared<call_vector_t>();
    call_vector_t &calls = *calls_;

    // std::cout << id << ":\n";
    // std::cout << seq << "\n";

    kguts->process_aa_seq(id, seq, calls_, 0, 0);

    /*
    for (auto x: *calls_)
    {
	std::cout << "  " << x << "\n";
    }
    */
	    
    int best_fi;
    RecallData ret;
    kguts->find_best_call(calls, best_fi, ret.new_function, ret.score, ret.weighted_score, ret.score_offset);

    ret.old_function = fm.lookup_function(id);

    return ret;
}

/*
 * Recall a fasta file using our just-computed signatures.
 *
 * calls_dir is where we write the new calls; new_dir is where we write
 * the changed annotationsl
 */
void recall_fasta(FunctionMap &fm, const fs::path &file, KmerGuts *kguts,
		  fs::path calls_dir, fs::path new_dir)
{
    fs::ifstream ifstr(file);

    fs::path calls_path = calls_dir / file.leaf();
    fs::path new_path = new_dir / file.leaf();

    fs::ofstream calls_stream(calls_path);
    fs::ofstream new_stream(new_path);

    FastaParser parser;
    
    auto cb = [&fm, kguts, &calls_stream, &new_stream](const std::string &id, const std::string &seq) {
	auto res = recall_sequence(fm, kguts, id, seq);
	if (res)
	{
	    if (res->new_function != res->old_function)
	    {
		new_stream << id << "\t" << res->old_function << "\t" << res->new_function << "\n";
	    }
	
	    calls_stream << id << "\t" << res->new_function << "\t" << res->score << "\t" << res->weighted_score << "\n";
	}
	return 0;
    };
	    
    parser.set_callback(cb);

    parser.parse(ifstr);
    parser.parse_complete();
}

/*
 * Validate a fasta file using our just-computed signatures.
 *
 */
void validate_fasta(FunctionMap &fm, const fs::path &file, KmerGuts *kguts,
		    FunctionMap &correct_calls, bool verbose)
{
    fs::ifstream ifstr(file);

    FastaParser parser;
    int n_correct = 0, n_incorrect = 0, n_missing = 0;
    int count = 0;
    
    auto cb = [&fm, kguts, &correct_calls, &n_correct, &n_incorrect, &n_missing, &count, verbose](const std::string &id, const std::string &seq) {
	auto res = recall_sequence(fm, kguts, id, seq);
	std::string correct_function = correct_calls.lookup_function(id);

	count++;
	
	if (res)
	{
	    if (res->new_function == correct_function)
	    {
		n_correct++;
	    }
	    else
	    {
		if (verbose)
		    std::cout << "incorrect\t" << id << "\t" << correct_function << "\t" << res->new_function << std::endl;
		n_incorrect++;
	    }
	}
	else
	{
	    if (!correct_function.empty())
		n_missing++;
	}
	return 0;
    };
	    
    parser.set_callback(cb);

    parser.parse(ifstr);
    parser.parse_complete();

    std::cout << file << ": count=" << count << " correct=" << n_correct << " incorrect=" << n_incorrect << " missing=" << n_missing << std::endl;
}


void show_ps()
{
    std::string cmd = "ps uwww" + std::to_string(getpid());
    system(cmd.c_str());
}

void populate_path_list(const std::vector<std::string> &dirs, std::vector<fs::path> &paths)
{
    for (auto dir: dirs)
    {
	fs::path p(dir);
	for (auto dit: fs::directory_iterator(dir))
	{
	    if (fs::is_regular_file(dit.path()))
	    {
		paths.emplace_back(dit);
	    }
	}
    }
}    

void load_strings(const std::vector<std::string> &files, std::vector<std::string> &strings)
{
    for (auto f: files)
    {
	std::ifstream ifstr(f);
	if (ifstr.good())
	{
	    // std::cout << "load " << f << "\n";
	    std::string line;
	    while (std::getline(ifstr, line, '\n'))
	    {
		strings.emplace_back(line);
	    }
	}
	else
	{
	    std::cerr << "could not open " << f << "\n";
	}
    }
}

bool process_command_line_options(int argc, char *argv[],
				  std::vector<fs::path> &function_definitions,
				  std::vector<fs::path> &fasta_data,
				  std::vector<fs::path> &fasta_data_kept_functions,
				  std::vector<std::string> &good_functions,
				  std::vector<std::string> &good_roles,
				  int &min_reps_required,
				  fs::path &recall_output_path,
				  fs::path &validation_folder_path,
				  bool &validation_verbose,
				  int &recall_min_hits,
				  int &recall_max_gap,
				  fs::path &kmer_data_dir,
				  fs::path &final_kmers,
				  int &n_threads)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options]\nAllowed options";
    po::options_description desc(x.str());

    std::vector<std::string> definition_dirs;
    std::vector<std::string> fasta_dirs;
    std::vector<std::string> fasta_keep_dirs;
    std::vector<std::string> good_function_files;
    std::vector<std::string> good_role_files;

    std::string recall_output;
    std::string validation_folder;

    recall_min_hits = 5;
    recall_max_gap = 200;

    n_threads = 1;

    desc.add_options()
	("definition-dir,D", po::value<std::vector<std::string>>(&definition_dirs)->multitoken(), "Directory of function definition files")
	("fasta-dir,F", po::value<std::vector<std::string>>(&fasta_dirs)->multitoken(), "Directory of fasta files of protein data")
	("fasta-keep-functions-dir,K", po::value<std::vector<std::string>>(&fasta_keep_dirs), "Directory of fasta files of protein data (keep functions defined here)")
	("good-functions", po::value<std::vector<std::string>>(&good_function_files), "File containing list of functions to be kept")
	("good-roles", po::value<std::vector<std::string>>(&good_role_files), "File containing list of roles to be kept")
	("kmer-data-dir", po::value<fs::path>(&kmer_data_dir), "Write kmer data files to this directory")
	("min-reps-required", po::value<int>(&min_reps_required), "Minimum number of genomes a function must be seen in to be considered for kmers")
	("final-kmers", po::value<fs::path>(&final_kmers), "Write final.kmers file to be consistent with km_build_Data")
	("recall-output", po::value<std::string>(&recall_output), "Recall proteins and write output to this path")
	("validation-folder", po::value<std::string>(&validation_folder), "Call proteins from this folder and write statistics")
	("validation-verbose", po::bool_switch(&validation_verbose)->default_value(false), "Emit verbose output during validation")
	("recall-min-hits", po::value<int>(&recall_min_hits), "min-hits parameter to use in recall")
	("recall-max-gap", po::value<int>(&recall_max_gap), "max_gaps parameter to use in recall")
	("n-threads", po::value<int>(&n_threads), "Number of threads to use")
	("help,h", "show this help message");

    po::variables_map vm;

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
	std::cout << desc << "\n";
	return false;
    }

    po::notify(vm);

    /*
     * Read definition and fasta dirs to populate the path lists.
     */
    populate_path_list(definition_dirs, function_definitions);
    populate_path_list(fasta_dirs, fasta_data);
    populate_path_list(fasta_keep_dirs, fasta_data_kept_functions);

    std::cout << "definitions: ";
    for (auto x: definition_dirs)
	std::cout << x << " ";
    std::cout << std::endl;
    std::cout << "fasta: ";
    for (auto x: fasta_dirs)
	std::cout << x << " ";
    std::cout << std::endl;
    std::cout << "keep: ";
    for (auto x: fasta_keep_dirs)
	std::cout << x << " ";
    std::cout << std::endl;

    load_strings(good_function_files, good_functions);
    load_strings(good_role_files, good_roles);

    if (!recall_output.empty())
    {
	recall_output_path = recall_output;
    }
    
    if (!validation_folder.empty())
    {
	validation_folder_path = validation_folder;
    }
    
    return true;
}

int main(int argc, char *argv[])
{
    KmerAttributeMap m;
    FunctionMap fm;
    unsigned next_sequence_id = 0;

    rejected_stream.open("rejected.by.build_signature_kmers");
    kept_stream.open("kept.by.build_signature_kmers");
    kept_function_stream.open("function.kept.log");

    std::vector<fs::path> function_definitions;
    std::vector<fs::path> fasta_data;
    std::vector<fs::path> fasta_data_kept_functions;

    std::vector<std::string> good_functions;
    std::vector<std::string> good_roles;

    fs::path recall_output_path;
    fs::path final_kmers;
    fs::path kmer_data_dir;

    fs::path validation_folder;
    bool validation_verbose = false;

    int min_reps_required = 5;

    int recall_min_hits;
    int recall_max_gap;
    
    int n_threads;

    if (!process_command_line_options(argc, argv,
				      function_definitions,
				      fasta_data,
				      fasta_data_kept_functions,
				      good_functions,
				      good_roles,
				      min_reps_required,
				      recall_output_path,
				      validation_folder,
				      validation_verbose,
				      recall_min_hits,
				      recall_max_gap,
				      kmer_data_dir,
				      final_kmers,
				      n_threads))
    {
	return 1;
    }

    /*
     * If we are going to recall, make sure our output path is appropriate.
     */
    fs::path recall_calls_dir, recall_new_dir;
    if (!recall_output_path.empty())
    {
	try {
	    if (fs::is_directory(recall_output_path))
	    {
		// Ensure we have Calls and New folders.
		recall_calls_dir = recall_output_path / "Calls";
		recall_new_dir = recall_output_path / "New";
		
		if (!fs::is_directory(recall_calls_dir))
		{
		    if (!fs::create_directory(recall_calls_dir))
		    {
			std::cerr << "Error creating " << recall_calls_dir << "\n";
		    }
		}
		if (!fs::is_directory(recall_new_dir))
		{
		    if (!fs::create_directory(recall_new_dir))
		    {
			std::cerr << "Error creating " << recall_new_dir << "\n";
		    }
		}
		
	    }
	    else
	    {
		std::cerr << "Recall path " << recall_output_path << " does not exist\n";
		exit(1);
	    }
	}
	catch (fs::filesystem_error  &e)
	{
	    std::cerr << "error setting up recall output: " << e.what() << "\n";
	}
    }

    /*
     * validate our kmer_data_dir
     */
    if (!kmer_data_dir.empty())
    {
	if (!fs::is_directory(kmer_data_dir))
	{
	    if (!fs::create_directory(kmer_data_dir))
	    {
		std::cerr << "Error creating " << kmer_data_dir << "\n";
		exit(1);
	    }
	}
    }

    tbb::task_scheduler_init sched_init(n_threads);

    g_kmer_processor = new KmerProcessor(n_threads);

    fm.add_good_roles(good_roles);
    fm.add_good_functions(good_functions);

    for (auto def: function_definitions)
    {
	fm.load_id_assignments(def);
    }

    std::vector<fs::path> all_fasta_data;

    std::cerr << "load fasta\n";

    for (auto fasta: fasta_data)
    {
	fm.load_fasta_file(fasta, false);
	all_fasta_data.emplace_back(fasta);
    }

    for (auto fasta: fasta_data_kept_functions)
    {
	fm.load_fasta_file(fasta, true);
	all_fasta_data.emplace_back(fasta);
    }

    // fm.dump();
    /*
     * Process the list of functions and
     * manage the set that we want to keep.
     */
    fm.process_kept_functions(min_reps_required);
    if (!kmer_data_dir.empty())
    {
	write_function_index(kmer_data_dir, fm);
	/* Write an empty otu index and genomes file (genomes file
	 * can't be empty because kmer_search uses "-s genomes" to test
	 * for its existence). */
	{
	    fs::ofstream otu(kmer_data_dir / "otu.index");
	    otu.close();
	    fs::ofstream genomes(kmer_data_dir / "genomes");
	    genomes << "empty genomes\n";
	    genomes.close();
	}
    }
    
    /*
     * With that done, go ahead and extract kmers.
     */

    if (n_threads < 0)
    {
	for (size_t i = 0; i < all_fasta_data.size(); i++)
	{
	    load_fasta(fm, m, i, all_fasta_data[i]);
	}
    }
    else
    {
	int n = all_fasta_data.size();
	tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
			  [&fm, &m, &next_sequence_id, &all_fasta_data](const tbb::blocked_range<size_t> &r) {
			      for (size_t i = r.begin(); i != r.end(); ++i)
			      {
				  auto fasta = all_fasta_data[i];
				  std::cout << "load file " << i << " " << fasta << "\n";
				  load_fasta(fm, m, i, fasta);
			      }
			  });
    }

    /*
    g_kmer_processor->start();
    process_kmers(m);
    g_kmer_processor->stop();
    */
    par_process_kmers(m);
    std::cout << "Kept " << kept_kmers.size() << " kmers\n";
    std::cout << "distinct_signatures=" << kmer_stats.distinct_signatures << "\n";
    std::cout << "num_seqs_with_a_signature=" << kmer_stats.seqs_with_a_signature.size() << "\n";
    std::cerr << "computing weights\n";
    std::for_each(kept_kmers.begin(), kept_kmers.end(), compute_weight_of_signature);
    std::cerr << "done\n";

    if (!final_kmers.empty())
    {
	std::cerr << "writing kmers to " << final_kmers << "\n";
	fs::ofstream kf(final_kmers);
	std::for_each(kept_kmers.begin(), kept_kmers.end(), [&kf](const KeptKmer &k) {
		kf << k.kmer << "\t" << k.median_offset << "\t" << k.function_index << "\t" << k.weight << "\t" << k.otu_index << "\n";
		// kf << "\t" << k.seqs_containing_sig << "\t" << kmer_stats.seqs_with_func[k.function_index] << "\n";
	    });
	std::cerr << "done\n";
    }

    std::cerr << "writing hashtable to " << kmer_data_dir << " ...\n";
    KmerGuts *kguts = write_hashtable(kmer_data_dir, kept_kmers);
    std::cerr << "writing hashtable to " << kmer_data_dir << " done\n";
    
    /*
     * Ick.
     */
    fs::path fidx = kmer_data_dir / "function.index";
    // kguts->kmersH->function_array = kguts->load_functions(fidx.c_str(), &kguts->kmersH->function_count);
    kguts->kmersH->function_array = fm.create_kg_function_array(&kguts->kmersH->function_count);
    kguts->kmersH->otu_array = kguts->load_otus("/dev/null", &kguts->kmersH->otu_count);

    if (!recall_output_path.empty())
    {
	kguts->min_hits = recall_min_hits;
	kguts->max_gap = recall_max_gap;

	KmerGuts *shared_kguts = kguts;
	int n = all_fasta_data.size();

	if (n_threads >= 1)
	{
	    std::cerr << "starting recall of " << n << " genomes with " << n_threads << " threads\n";
	    tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
			      [&fm, &shared_kguts, recall_calls_dir, recall_new_dir, all_fasta_data, kmer_data_dir,
			       recall_min_hits, recall_max_gap, fidx](const tbb::blocked_range<size_t> &r) {

				  //KmerGuts *kguts = new KmerGuts(dir, shared_kguts->image_);
				  KmerGuts kguts(*shared_kguts);
				  //boost::thread_specific_ptr<KmerGuts> kguts;
				  //kguts.reset(new KmerGuts(*shared_kguts));

				  //kguts->min_hits = recall_min_hits;
				  //kguts->max_gap = recall_max_gap;
				  //kguts->kmersH->function_array = kguts->load_functions(fidx.c_str(), &kguts->kmersH->function_count);
				  //kguts->kmersH->otu_array = kguts->load_otus("/dev/null", &kguts->kmersH->otu_count);
				  
				  std::cout << "Process block:";
				  for(size_t i=r.begin(); i!=r.end(); ++i)
				      std::cout << " " << i;
				  std::cout << " on thread " << std::this_thread::get_id() << "\n";
				  
				  for(size_t i=r.begin(); i!=r.end(); ++i)
				  {
				      fs::path file = all_fasta_data[i];
				      recall_fasta(fm, file, &kguts, recall_calls_dir, recall_new_dir);
				  }
			      });
	}
	else
	{
	    std::cerr << "starting recall of " << n << " genomes with single thread\n";
	    for (auto file: all_fasta_data)
	    {
		//KmerGuts k(*kguts);
		//KmerGuts &k = *kguts;
		recall_fasta(fm, file, kguts, recall_calls_dir, recall_new_dir);
	    }
	}
    }

    if (!validation_folder.empty())
    {
	kguts->min_hits = recall_min_hits;
	kguts->max_gap = recall_max_gap;

	KmerGuts *shared_kguts = kguts;

	/*
	 * The validation folder will contain a subdirectory
	 * anno that holds the 2-column table of correct annotations (file per genome), and
	 * seq that holds the fasta sequence data (file per genome)
	 *
	 * For each anno file, read the 2-column table and populate the FunctionMap
	 * we're (mis)using for the id => correct function hash.
	 *
	 * Read the contents of the seq directory to obtain a list of files.
	 * 
	 * In parallel, recall and compute scores.
	 *  
	 */

	std::vector<fs::path> seq_files;
	FunctionMap correct_funcs;
	
	for (auto dit: fs::directory_iterator(validation_folder / "seq"))
	{
	    if (fs::is_regular_file(dit.path()))
		seq_files.push_back(dit.path());
	}
	for (auto dit: fs::directory_iterator(validation_folder / "anno"))
	{
	    if (fs::is_regular_file(dit.path()))
	    {
		correct_funcs.load_id_assignments(dit.path());
	    }
	}

	std::cerr << "starting call of " << seq_files.size() << " genomes with " << n_threads << " threads\n";
	tbb::parallel_for(tbb::blocked_range<size_t>(0, seq_files.size()),
			  [&fm, &shared_kguts, &seq_files, recall_min_hits, recall_max_gap, &correct_funcs, validation_verbose]
			  (const tbb::blocked_range<size_t> &r) {

			      KmerGuts kguts(*shared_kguts);
			      
			      std::cout << "Process block:";
			      for(size_t i=r.begin(); i!=r.end(); ++i)
				  std::cout << " " << i;
			      std::cout << " on thread " << std::this_thread::get_id() << "\n";
			      
			      for(size_t i=r.begin(); i!=r.end(); ++i)
			      {
				  fs::path file = seq_files[i];
				  std::cerr << "Would compute " << file << std::endl;
				  validate_fasta(fm, file, &kguts, correct_funcs, validation_verbose);
			      }
			  });
    }

    std::cerr << "all done\n";

    free(kguts->kmer_image_for_loading_);
    free(kguts->kmersH->function_array);
    free(kguts->kmersH->otu_array);
    delete kguts->kmersH;
    delete kguts;

    show_ps();
    rejected_stream.close();
    kept_stream.close();
    kept_function_stream.close();
    return 0;
}
