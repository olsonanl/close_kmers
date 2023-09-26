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
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <boost/thread/tss.hpp>
#include <boost/thread/thread.hpp>

#include "nudb_kmer_db.h"
#include "kmer_nudb.h"
#include "seed_utils.h"
#include "operators.h"
#include "fasta_parser.h"

#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/compat/thread"
#include "tbb/concurrent_hash_map.h"

using namespace seed_utils;

#define DEFINE_GLOBALS
#include "global.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include "function_map.h"

const int K = KMER_SIZE;
const int MaxSequencesPerFile = 100000;

// Kmer data type
typedef std::array<char, K> Kmer;

// Hash function on kmers
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
    FunctionIndex func_index;
    OTUIndex otu_index;
    unsigned short offset;
    unsigned int seq_id;
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
	    h = (h * 17) ^ (unsigned int) s;
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
    tbb::concurrent_unordered_map<FunctionIndex, int> seqs_with_func;
    tbb::concurrent_unordered_set<unsigned int> seqs_with_a_signature;
};

//
// Whack some globals in here. We're maintaining enough state this
// needs to be wrapped up in an object.
//

KmerStatistics kmer_stats;

tbb::mutex io_mutex;
std::ofstream rejected_stream;
std::ofstream kept_stream;

/*
 * Object to keep state for kmers we are keeping. We hang onto
 * some statistics in order to compute weights later on.
 */
struct KeptKmer
{
    Kmer kmer;
    unsigned short median_offset;		// Median offset from the end of the protein
    FunctionIndex function_index;
//    FunctionIndex function_index2;	// Secondary function in the case of ambiguity
    OTUIndex  otu_index;
    unsigned int seqs_containing_sig;	// Count of sequences containing this kmer
    unsigned int seqs_containing_function; // Count of sequences with the kmer that have the function
    float weight;
};

inline std::ostream &operator<<(std::ostream &os, const KeptKmer &k)
{
    os << k.kmer << " " << k.function_index << " " <<  k.seqs_containing_sig << " " << k.seqs_containing_function;
    return os;
}


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
    std::map<FunctionIndex, int> func_count;
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

std::set<unsigned char> ok_prot = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
			   'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y'};

void load_sequence(FunctionMap &fm, KmerAttributeMap &m, unsigned int &next_sequence_id,
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
    
    unsigned int seq_id = next_sequence_id++;

    FunctionIndex function_index = fm.lookup_index(func);

    if (function_index == UndefinedFunction)
    {
	function_index = fm.lookup_index("hypothetical protein");
	if (function_index == UndefinedFunction)
	{
	    std::cerr << "No function defined for hypothetical protein\n";
	    exit(1);
	}
    }

    // std::cout << "Load seq " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
    if (function_index == UndefinedFunction)
	return;

    kmer_stats.seqs_with_func[function_index]++;

    for (auto it = seq.begin(); it < seq.end() - K + 1; it++)
    {
        unsigned short n = (unsigned short) std::distance(it, seq.end());
	Kmer kmer;
	// Kmer::iterator kiter = kmer.begin();
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
	    m.insert({kmer, { function_index, UndefinedOTU, n, seq_id }});
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

//void process_set(Kmer &kmer, std::map<FunctionIndex, int> &func_count, int count, std::vector<KmerAttributes> &set)

void process_set(KmerSet &set)
{
    FunctionIndex best_func_1 = UndefinedFunction, best_func_2 = UndefinedFunction;
    int best_count_1 = -1, best_count_2 = -1;
    
    /*
    auto elt = std::max_element(set.func_count.begin(), set.func_count.end(),
				[](auto a, auto b) { return a.second < b.second; });

    FunctionIndex best_func = elt->first;
    int best_count = elt->second;
    */

    // we want the top two elements by value in the map; don't know how
    // with standard STL without copying to a vector, but if we're copying
    // anyway we can just search for them.
    for (auto x: set.func_count)
    {
	if (best_func_1 == UndefinedFunction)
	{
	    best_func_1 = x.first;
	    best_count_1 = x.second;
	}
	else if (x.second > best_count_1)
	{
	    best_func_2 = best_func_1;
	    best_count_2 = best_count_1;

	    best_func_1 = x.first;
	    best_count_1 = x.second;
	}
	else if (x.second > best_count_2)
	{
	    best_func_2 = x.first;
	    best_count_2 = x.second;
	}
    }

    float thresh = float(set.count) * 0.8f;
    int best_count = best_count_1;
    FunctionIndex best_func = best_func_1;

    if (rejected_stream.is_open() || kept_stream.is_open())
    {
	tbb::mutex::scoped_lock lock(io_mutex);

	/*
	kept_stream << "Process set for " << set.kmer << " best1=" << best_func_1 << " best_count_1=" << best_count_1
		    << " best2=" << best_func_2 << " best_count_2=" << best_count_2 << " count=" << set.count << " thresh=" << thresh << "\n";
	*/

/*
	for (auto x: set.func_count)
	{
	    kept_stream << x.first << " " << x.second << "\n";
	}
*/
	   
	/*
	std::sort(set.set.begin(), set.set.end(), [&set](auto a, auto b) {
		return set.func_count[b.func_index] < set.func_count[a.func_index];
	    });

	for (auto x: set.set)
	{
	    kept_stream << "  " << x.func_index << " " << x.seq_id << "\n";
	}
	*/
	
	if ((float) best_count < thresh)
	{
	    // std::cout << "discard for " << best_count << " < " << thresh << "\n";
	    if (rejected_stream.is_open())
		rejected_stream << set.kmer <<  " best_count=" << best_count << " thresh=" << thresh << "\n";
	    if ((float) (best_count_1 + best_count_2) >= thresh)
	    {
		if (kept_stream.is_open())
		    kept_stream << "AMBIG\t" << best_func_1 << "\t"
				<< best_func_2 << "\t"
				<< best_count_1 << "\t"
				<< best_count_2 << "\n";
	    }
	    return;
	}
    }
    else
    {
	if ((float) best_count < thresh)
	    return;
    }

    unsigned int seqs_containing_func = 0;
    std::vector<unsigned short> offsets;

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
    unsigned short median_offset = offsets[offsets.size() / 2];
    // std::cout << seqs_containing_func << " " << median_offset<< "\n";

    kmer_stats.distinct_signatures++;
    kmer_stats.distinct_functions[best_func]++;

    kept_kmers.push_back(KeptKmer { set.kmer, median_offset, best_func, UndefinedOTU,
		(unsigned int) set.set.size(), seqs_containing_func });
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

#if 0
// Not used apparently
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
#endif

void compute_weight_of_signature(KeptKmer &kk)
{
    float NSF = (float) kmer_stats.seqs_with_a_signature.size();
    float KS  = (float) kmer_stats.distinct_signatures;
    // float KF  = kmer_stats.seqs_with_func.size();
    float NSi = (float) kk.seqs_containing_sig;
    float NFj = (float) kmer_stats.seqs_with_func[kk.function_index];
    float NSiFj = (float) kk.seqs_containing_function;

    kk.weight = std::log((NSiFj + 1.0f) / (NSi - NSiFj + 1.0f)) +
	std::log((NSF - NFj + KS) / (NFj + KS));

}

void write_function_index(const fs::path &dir, FunctionMap &fm)
{
    fm.write_function_index(dir);
}

struct RecallData
{
    std::string new_function;
    std::string old_function;
    float score;
    float weighted_score;
    float score_offset;
};

std::experimental::optional<RecallData> recall_sequence(FunctionMap &fm, KmerNudb *db,
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

    db->process_aa_seq(id, seq, calls_, 0, 0);

    /*
    for (auto x: *calls_)
    {
	std::cout << "  " << x << "\n";
    }
    */
	    
    FunctionIndex best_fi;
    RecallData ret;

    db->find_best_call(calls, best_fi, ret.new_function, ret.score, ret.weighted_score, ret.score_offset);

    ret.old_function = fm.lookup_function(id);

    return ret;
}

/*
 * Recall a fasta file using our just-computed signatures.
 *
 * calls_dir is where we write the new calls; new_dir is where we write
 * the changed annotationsl
 */
template<typename Caller>
void recall_fasta(FunctionMap &fm, const fs::path &file, Caller *caller,
		  fs::path calls_dir, fs::path new_dir)
{
    fs::ifstream ifstr(file);

    fs::path calls_path = calls_dir / file.leaf();
    fs::path new_path = new_dir / file.leaf();

    fs::ofstream calls_stream(calls_path);
    fs::ofstream new_stream(new_path);

    FastaParser parser;
    
    auto cb = [&fm, caller, &calls_stream, &new_stream](const std::string &id, const std::string &seq) {
	auto res = recall_sequence(fm, caller, id, seq);
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
#if 0
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
#endif


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
				  fs::path &deleted_fids_file,
				  int &min_reps_required,
				  fs::path &recall_output_path,
				  fs::path &validation_folder_path,
				  bool &validation_verbose,
				  int &recall_min_hits,
				  int &recall_max_gap,
				  fs::path &kmer_data_dir,
				  fs::path &final_kmers,
				  std::string &nudb_file,
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
	("deleted-features-file", po::value<fs::path>(&deleted_fids_file), "File containing list of deleted feature IDs")
	("kmer-data-dir", po::value<fs::path>(&kmer_data_dir), "Write kmer data files to this directory")
	("nudb-file", po::value<std::string>(&nudb_file), "Write saved kmers to this NuDB file base. Should be on a SSD drive.")
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

    // rejected_stream.open("/dev/shm/rejected.by.build_signature_kmers");
    // kept_stream.open("kept.by.build_signature_kmers");
    // kept_function_stream.open("function.kept.log");

    std::vector<fs::path> function_definitions;
    std::vector<fs::path> fasta_data;
    std::vector<fs::path> fasta_data_kept_functions;

    std::vector<std::string> good_functions;
    std::vector<std::string> good_roles;

    fs::path recall_output_path;
    fs::path final_kmers;
    fs::path deleted_fids_file;
    fs::path kmer_data_dir;

    fs::path validation_folder;
    bool validation_verbose = false;

    int min_reps_required = 3;

    int recall_min_hits;
    int recall_max_gap;
    
    int n_threads;

    std::string nudb_file;

    if (!process_command_line_options(argc, argv,
				      function_definitions,
				      fasta_data,
				      fasta_data_kept_functions,
				      good_functions,
				      good_roles,
				      deleted_fids_file,
				      min_reps_required,
				      recall_output_path,
				      validation_folder,
				      validation_verbose,
				      recall_min_hits,
				      recall_max_gap,
				      kmer_data_dir,
				      final_kmers,
				      nudb_file,
				      n_threads))
    {
	return 1;
    }

    std::set<std::string> deleted_fids;
    /*
     * Read deleted fids if present.
     */
    if (!deleted_fids_file.empty())
    {
	fs::ifstream ifstr(deleted_fids_file);
	std::string line;
	while (std::getline(ifstr, line, '\n'))
	{
	    std::cerr << "'" << line << "'\n";
	    deleted_fids.emplace(line);
	}
	
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
	fm.load_fasta_file(fasta, false, deleted_fids);
	all_fasta_data.emplace_back(fasta);
    }

    for (auto fasta: fasta_data_kept_functions)
    {
	fm.load_fasta_file(fasta, true, deleted_fids);
	all_fasta_data.emplace_back(fasta);
    }

    fs::path fidx = kmer_data_dir / "function.index";
    NuDBKmerDb<K> db(nudb_file);
    db.open();

    if (!recall_output_path.empty())
    {
	size_t n = all_fasta_data.size();

	if (n_threads >= 1)
	{
	    std::cerr << "starting recall of " << n << " genomes with " << n_threads << " threads\n";
	    tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
			      [&fm,
			       &db,
			       recall_calls_dir, recall_new_dir, all_fasta_data, kmer_data_dir,
			       recall_min_hits, recall_max_gap, fidx](const tbb::blocked_range<size_t> &r) {

				  KmerNudb kmer_nudb(db, fidx.string(), recall_min_hits, 0.0, recall_max_gap);

				  std::cout << "Process block:";
				  for(size_t i=r.begin(); i!=r.end(); ++i)
				      std::cout << " " << i;
				  std::cout << " on thread " << std::this_thread::get_id() << "\n";
				  
				  for(size_t i=r.begin(); i!=r.end(); ++i)
				  {
				      fs::path file = all_fasta_data[i];
				      recall_fasta(fm, file, &kmer_nudb, recall_calls_dir, recall_new_dir);
				  }
			      });
	}
	else
	{
	    std::cerr << "starting recall of " << n << " genomes with single thread\n";
	    for (auto file: all_fasta_data)
	    {
		KmerNudb kmer_nudb(db, fidx.string(), recall_min_hits, 0.0, recall_max_gap);
		recall_fasta(fm, file, &kmer_nudb, recall_calls_dir, recall_new_dir);
	    }
	}
    }

#if 0
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
#endif

    std::cerr << "all done\n";

/*
    free(kguts->kmer_image_for_loading_);
    delete kguts;
*/
    show_ps();
    rejected_stream.close();
    kept_stream.close();
    // kept_function_stream.close();
    return 0;
}
