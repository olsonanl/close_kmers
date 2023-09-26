
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
#include <mutex>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <boost/thread/tss.hpp>
#include <boost/thread/thread.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "nudb_kmer_db.h"
#include "seed_utils.h"
#include "operators.h"
#include "fasta_parser.h"
#include "kmer_types_generic.h"

#define TBB_PREVIEW_NUMA_SUPPORT 1

#include "tbb/global_control.h"
#include "tbb/info.h"
#include "tbb/parallel_for.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_hash_map.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_queue.h"

#include "welford.h"

using namespace seed_utils;

#define DEFINE_GLOBALS
#include "global.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace acc = boost::accumulators;

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
    unsigned int protein_length;
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

std::mutex io_mutex;
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
    unsigned short median_offset;		// Median offset from the end of the protein
    FunctionIndex function_index;
//    FunctionIndex function_index2;	// Secondary function in the case of ambiguity
    OTUIndex  otu_index;
    unsigned int seqs_containing_sig;	// Count of sequences containing this kmer
    unsigned int seqs_containing_function; // Count of sequences with the kmer that have the function
    float weight;
    unsigned short mean;
    unsigned short median;
    unsigned short var;
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
    for (auto iter = ks.func_count.begin(); iter != ks.func_count.end(); iter++)
    {
	os << iter->first << "\t" << iter->second << "\n";
    }
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

    if (false)
    {
	if (function_index == UndefinedFunction)
	{
	    function_index = fm.lookup_index("hypothetical protein");
	    if (function_index == UndefinedFunction)
	    {
		std::cerr << "No function defined for hypothetical protein\n";
		exit(1);
	    }
	}
    }

    // std::cout << "Load seq " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
//    if (function_index == 19345)
//    {
//	std::cout << "FN 19345 " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
//    }
    // std::cout << "Load seq " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
    if (function_index == UndefinedFunction)
    {
//	std::cout << "Skipping undef seq " << id << " '" << func << "' " << function_index << " id=" << seq_id <<  "\n";
    	return;
    }

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
	    /*
	    if ("PHSQRWTN" == kmer || function_index == 19345)
	    {
		std::cerr << kmer << " in " << id << " " << function_index << " " << seq_id << "\n";
	    }
	    */
	    m.insert({kmer, { function_index, UndefinedOTU, n, seq_id, static_cast<unsigned int>(seq.length())}});
	    // std::cout << kmer << " " << n << "\n";
	}
	else
	{
	    /*
	    if ("PHSQRWTN" == kmer || function_index == 19345)
	    {
		std::cerr << kmer << " rejected from " << id << " " << function_index << " " << seq_id << "\n";
	    }
	    */
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
void load_fasta(FunctionMap &fm, KmerAttributeMap &m, unsigned file_number, const fs::path &file,
		const std::set<std::string> &deleted_fids)
{
    fs::ifstream ifstr(file);

    FastaParser parser;
    
    unsigned next_sequence_id = file_number * MaxSequencesPerFile;

    parser.set_def_callback([&fm, &m, &next_sequence_id, &deleted_fids](const std::string &id, const std::string &def, const std::string &seq) {
	if (deleted_fids.find(id) == deleted_fids.end())
	{
	    load_sequence(fm, m, next_sequence_id, id, def, seq);
	}
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

    if ("PHSQRWTN" == set.kmer)
    {
	std::cerr << set;
	std::cerr << "best=" << best_func_1 << " best2=" << best_func_2 << "\n";
    }


    float thresh = float(set.count) * 0.8f;
    int best_count = best_count_1;
    FunctionIndex best_func = best_func_1;

    if (rejected_stream.is_open() || kept_stream.is_open())
    {
	std::lock_guard<std::mutex> lock(io_mutex);

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
	{
	    if ("PHSQRWTN" == set.kmer)
	    {
		std::cerr << "rejecting " << set.kmer << " best_count=" << best_count << " thresh=" << thresh << "\n";
	    }
	    return;
	}
    }

    unsigned int seqs_containing_func = 0;
    std::vector<unsigned short> offsets;

    acc::accumulator_set<unsigned short, acc::stats<acc::tag::mean,
						  acc::tag::median,
						  acc::tag::variance> > acc;

    for (auto item: set.set)
    {
	if (item.func_index == best_func)
	{
	    seqs_containing_func++;
	    acc(item.protein_length);
	}
	offsets.push_back(item.offset);
	kmer_stats.seqs_with_a_signature.insert(item.seq_id);
    }

    unsigned short mean = acc::mean(acc);
    unsigned short median = acc::median(acc);
    unsigned short var = acc::variance(acc);

    std::sort(offsets.begin(), offsets.end());
    unsigned short median_offset = offsets[offsets.size() / 2];
    // std::cout << seqs_containing_func << " " << median_offset<< "\n";

    kmer_stats.distinct_signatures++;
    kmer_stats.distinct_functions[best_func]++;

    kept_kmers.emplace_back(KeptKmer { set.kmer, median_offset, best_func, UndefinedOTU,
	    (unsigned int) set.set.size(), seqs_containing_func, 0.0, mean, median, var });
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

void write_nudb_data(const std::string &nudb_file, tbb::concurrent_vector<KeptKmer> &kmers)
{
    typedef NuDBKmerDb<8> KDB;

    KDB db(nudb_file);

    if (!db.exists())
    {
	std::cerr << "creating new db\n";
	db.create();
    }
    db.open();
    
    for (auto k: kmers)
    {
	// std::cout << k << "\n";
	nudb::error_code ec;
//	db.insert(k.kmer, { k.otu_index, k.median_offset, k.function_index, k.weight }, ec);
	db.insert(k.kmer, { k.otu_index, k.median_offset, k.function_index, k.weight,
		k.mean, k.median, k.var, static_cast<unsigned short>(k.seqs_containing_function) }, ec);
//	if (!ec)
//	    std::cerr << "insert error: " << ec.message() << "\n";
    }
}


void show_ps()
{
    std::string cmd = "ps uwww" + std::to_string(getpid());
    int rc = system(cmd.c_str());
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

    return true;
}

int main(int argc, char *argv[])
{
    KmerAttributeMap m;
    unsigned next_sequence_id = 0;

    // rejected_stream.open("/dev/shm/rejected.by.build_signature_kmers");
    // kept_stream.open("kept.by.build_signature_kmers");
    // kept_function_stream.open("function.kept.log");

    std::vector<fs::path> function_definitions;
    std::vector<fs::path> fasta_data;
    std::vector<fs::path> fasta_data_kept_functions;

    std::vector<std::string> good_functions;
    std::vector<std::string> good_roles;

    fs::path final_kmers;
    fs::path deleted_fids_file;
    fs::path kmer_data_dir;

    int min_reps_required = 3;
    
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
				      kmer_data_dir,
				      final_kmers,
				      nudb_file,
				      n_threads))
    {
	return 1;
    }

    FunctionMap fm((kmer_data_dir / "kept_functions.log").string());

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
	    deleted_fids.emplace(line);
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

    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, n_threads);

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
    
    if (0)
    {
	fm.dump();
	exit(0);
    }
    /*
     * With that done, go ahead and extract kmers.
     */

    if (n_threads < 0)
    {
	for (unsigned i = 0; i < (unsigned) all_fasta_data.size(); i++)
	{
	    load_fasta(fm, m, i, all_fasta_data[i], deleted_fids);
	}
    }
    else
    {
	size_t n = all_fasta_data.size();
	tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
			  [&fm, &m, &next_sequence_id, &all_fasta_data, &deleted_fids](const tbb::blocked_range<size_t> &r) {
			      for (size_t i = r.begin(); i != r.end(); ++i)
			      {
				  auto fasta = all_fasta_data[i];
				  // std::cout << "load file " << i << " " << fasta << "\n";
				  load_fasta(fm, m, (unsigned) i, fasta, deleted_fids);
			      }
			  });
    }

    /*
    g_kmer_processor->start();
    process_kmers(m);
    g_kmer_processor->stop();
    */
    std::cerr << "processing kmers\n";
    par_process_kmers(m);
    std::cout << "Kept " << kept_kmers.size() << " kmers\n";
    std::cout << "distinct_signatures=" << kmer_stats.distinct_signatures << "\n";
    std::cout << "num_seqs_with_a_signature=" << kmer_stats.seqs_with_a_signature.size() << "\n";
    std::cerr << "computing weights\n";
//    std::for_each(kept_kmers.begin(), kept_kmers.end(), compute_weight_of_signature);

    size_t nk = kept_kmers.size();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nk),
			  [](const tbb::blocked_range<size_t> &r) {
			      for (size_t i = r.begin(); i != r.end(); ++i)
			      {
				  compute_weight_of_signature(kept_kmers[i]);
			      }
			  });

    
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

    /*    
    std::cerr << "writing hashtable to " << kmer_data_dir << " ...\n";
    KmerGuts *kguts = write_hashtable(kmer_data_dir, kept_kmers);
    std::cerr << "writing hashtable to " << kmer_data_dir << " done\n";
    */

    if (!nudb_file.empty())
    {
	write_nudb_data(nudb_file, kept_kmers);
    }
    
    NuDBKmerDb<K> db(nudb_file);
    db.open();


    std::cerr << "all done\n";

    show_ps();
    rejected_stream.close();
    kept_stream.close();
    kept_function_stream.close();
    return 0;
}
