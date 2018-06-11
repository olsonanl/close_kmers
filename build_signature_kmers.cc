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
#include <unistd.h>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "seed_utils.h"
#include "operators.h"
#include "fasta_parser.h"
#include "kguts.h"

using namespace seed_utils;

#define DEFINE_GLOBALS
#include "global.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

const int K = 8;

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
    int seq_id;
};

typedef std::unordered_multimap<Kmer, KmerAttributes> KmerAttributeMap;

inline std::ostream &operator<<(std::ostream &os, const Kmer &k)
{
    std::copy(k.begin(), k.end(),
	      std::ostream_iterator<char>(os, ""));
    return os;
}

struct KmerStatistics
{
    int distinct_signatures = 0;
    std::unordered_map<int, int> distinct_functions;
    std::unordered_map<int, int> seqs_with_func;
    std::unordered_set<int> seqs_with_a_signature;
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

std::vector<KeptKmer> kept_kmers;

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

	    id_function_map_[id] = strip_func_comment(func);
	}
    }

    /*
     * Load assignments and sequence visibility data from a fasta file.
     *
     * If the def line has an assignment on it, updated the assignment map.
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
	
	fs::ifstream ifstr(file);

	FastaParser parser;

	std::string genome;

	parser.set_def_callback([this, &genome, &genome_regex, &figid_regex, &file, keep_function_flag]
				(const std::string &id, const std::string &def, const std::string &seq) {
		if (id.empty())
		    return 0;

		boost::smatch match;

		//
		// Need to always parse for [genome]
		//
		
		std::string func = def;
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
			else
			{
			    std::cerr << "cannot determine genome from file " << file << "\n";
			    // default it to the file, just to have a value
			    genome = file.leaf().string();
			}
		    }
		    else
		    {
			if (genome_loc.empty())
			{
			    std::cerr << "cannot determine genome from file " << file << "\n";
			    // default it to the file, just to have a value
			    genome = file.leaf().string();
			}
			else
			{
			    genome = genome_loc;
			}
		    }
		}

		/*
		 * If we're assigning a function, update the id to function
		 * map.
		 *
		 * Then look up the function for this id and add to the function_genome map.
		 */
		if (func.empty())
		{
		    func = id_function_map_[id];
		}
		else
		{
		    id_function_map_[id] = func;
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
	    std::cout << x.first << " " << x.second << "\n";
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

    void write_function_index(const std::string &dir) {
	std::ofstream of(dir + "/function.index");

	std::map<int, std::string> by_index;
	for (auto ent: function_index_map_)
	    by_index.insert(std::make_pair(ent.second, ent.first));
	for (auto ent: by_index)
	{
	    of << ent.first << "\t" << ent.second << "\n";
	}
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

void load_sequence(FunctionMap &fm, KmerAttributeMap &m, int &next_sequence_id,
		   const std::string &id, const std::string &def, const std::string &seq)
{
    if (id.empty())
	return;

    std::string func = def;
    if (func.empty())
    {
	func = fm.lookup_function(id);
    }
    /*
     * Empty means empty (and perhaps deleted feature).
     */
    
    if (func.empty())
    {
	return;
    }
    
    int seq_id = next_sequence_id++;

    int function_index = fm.lookup_index(func);
    // std::cout << "Load seq " << id << " '" << func << "' " << function_index << "\n";
    if (function_index < 0)
	return;

    kmer_stats.seqs_with_func[function_index]++;

    for (auto it = seq.begin(); it < seq.end() - K + 1; it++)
    {
        int n = std::distance(it, seq.end());
	Kmer kmer;
	Kmer::iterator kiter = kmer.begin();
	auto this_char = it;
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
void load_fasta(FunctionMap &fm, KmerAttributeMap &m, int &next_sequence_id, const fs::path &file)
{
    fs::ifstream ifstr(file);

    FastaParser parser;

    parser.set_def_callback([&fm, &m, &next_sequence_id](const std::string &id, const std::string &def, const std::string &seq) {
	    load_sequence(fm, m, next_sequence_id, id, def, seq);
	    return 0;
	});
    parser.parse(ifstr);
    parser.parse_complete();
}

void process_set(Kmer &kmer, std::map<int, int> &func_count, int count, std::vector<KmerAttributes> &set)
{
    /*
    if (count < 2)
    {
    rejected_stream << kmer << " count=" << count << "\n";
	return;
    }
    */

    auto elt = std::max_element(func_count.begin(), func_count.end(),
				[](auto a, auto b) { return a.second < b.second; });

    int best_func = elt->first;
    int best_count = elt->second;

    float thresh = float(count) * 0.8;

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

    for (auto item: set)
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

    kept_kmers.emplace_back(KeptKmer { kmer, median_offset, best_func, -1,
				       (int) set.size(), seqs_containing_func });
}


void process_kmers(KmerAttributeMap &m)
{
    std::map<int, int> func_count;
    std::vector<KmerAttributes> set;
    int count = 0;
    Kmer cur { 0 } ;
    
    for (auto ent: m)
    {
	const Kmer &kmer = ent.first;
	KmerAttributes &attr = ent.second;

	if (kmer != cur)
	{
	    if (count > 0)
		process_set(cur, func_count, count, set);
	    func_count.clear();
	    set.clear();
	    count = 0;
	    cur = kmer;
	}

	func_count[attr.func_index]++;
	count++;
	set.emplace_back(attr);
    }
    process_set(cur, func_count, count, set);
    
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
	log((NSF - NFj + KS) / NFj + KS);

}

void write_function_index(const std::string &dir, FunctionMap &fm)
{
    fm.write_function_index(dir);
}

KmerGuts *write_hashtable(const std::string &dir, std::vector<KeptKmer> &kmers)
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
	
    KmerGuts *kguts = new KmerGuts(dir, hashtable_size);

    for (auto k: kmers)
    {
	std::string kstr(k.kmer.begin(), k.kmer.end());
	kguts->insert_kmer(kstr,
			  k.function_index, k.otu_index, k.median_offset, k.weight);
    }
    kguts->save_kmer_hash_table(dir + "/kmer.table.mem_map");
    return kguts;
}

int recall_sequence(FunctionMap &fm, KmerGuts *kguts, const std::string &id, const std::string &seq)
{

    if (id.empty())
	return 0;
    typedef std::vector<KmerCall> call_vector_t;
    std::shared_ptr<call_vector_t> calls_ = std::make_shared<call_vector_t>();
    kguts->process_aa_seq(id, seq, calls_, 0, 0);
	    
    int best_fi;
    float best_score;
    std::string best_func;
    kguts->find_best_call(*calls_, best_fi, best_func, best_score);

    auto old = fm.lookup_function(id);
    bool diff = (old != best_func);

    std::cout << id << "\t" << diff << "\t" << old << "\t" << best_func << "\n";
	    
    return 0;
}

/*
 * Recall a fasta file using our just-computed signatures.
 */
void recall_fasta(FunctionMap &fm, const fs::path &file, KmerGuts *kguts)
{
    fs::ifstream ifstr(file);

    FastaParser parser;

    using namespace std::placeholders;

   parser.set_callback(std::bind(&recall_sequence, fm, kguts, _1, _2));

/*    parser.set_callback([&fm, kguts](const std::string &id, const std::string &seq) {

	    if (id.empty())
		return 0;
	    std::cout << "Process " << id << " " << seq << " with " << kguts << "\n";
	    typedef std::vector<KmerCall> call_vector_t;
	    std::shared_ptr<call_vector_t> calls_ = std::make_shared<call_vector_t>();
	    kguts->process_aa_seq(id, seq, calls_, 0, 0);
	    
	    int best_fi;
	    float best_score;
	    std::string best_func;
	    kguts->find_best_call(*calls_, best_fi, best_func, best_score);
	    std::cout << id << ": " << best_func << " " << best_score << "\n";
	    
	    return 0;
	});
*/
    parser.parse(ifstr);
    parser.parse_complete();
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
				  int &min_reps_required)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options]\nAllowed options";
    po::options_description desc(x.str());

    std::vector<std::string> definition_dirs;
    std::vector<std::string> fasta_dirs;
    std::vector<std::string> fasta_keep_dirs;
    std::vector<std::string> good_function_files;
    std::vector<std::string> good_role_files;

    desc.add_options()
	("definition-dir,D", po::value<std::vector<std::string>>(&definition_dirs), "Directory of function definition files")
	("fasta-dir,F", po::value<std::vector<std::string>>(&fasta_dirs), "Directory of fasta files of protein data")
	("fasta-keep-functions-dir,K", po::value<std::vector<std::string>>(&fasta_keep_dirs), "Directory of fasta files of protein data (keep functions defined here)")
	("good-functions", po::value<std::vector<std::string>>(&good_function_files), "File containing list of functions to be kept")
	("good-roles", po::value<std::vector<std::string>>(&good_role_files), "File containing list of roles to be kept")
	("min-reps-required", po::value<int>(&min_reps_required), "Minimum number of genomes a function must be seen in to be considered for kmers")
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

    load_strings(good_function_files, good_functions);
    load_strings(good_role_files, good_roles);
    
    return true;
}

int main(int argc, char *argv[])
{
    KmerAttributeMap m;
    FunctionMap fm;
    int next_sequence_id = 0;

    rejected_stream.open("rejected.by.build_signature_kmers");
    kept_stream.open("kept.by.build_signature_kmers");
    kept_function_stream.open("function.kept.log");

    std::vector<fs::path> function_definitions;
    std::vector<fs::path> fasta_data;
    std::vector<fs::path> fasta_data_kept_functions;

    std::vector<std::string> good_functions;
    std::vector<std::string> good_roles;

    int min_reps_required = 5;

    if (!process_command_line_options(argc, argv,
				      function_definitions,
				      fasta_data,
				      fasta_data_kept_functions,
				      good_functions,
				      good_roles,
				      min_reps_required))
    {
	return 1;
    }

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

    //fm.dump();
    /*
     * Process the list of functions and
     * manage the set that we want to keep.
     */
    fm.process_kept_functions(min_reps_required);
    std::string dir("tmp_out");
    write_function_index(dir, fm);

    /*
     * With that done, go ahead and extract kmers.
     */
       
    for (auto fasta: all_fasta_data)
    {
	load_fasta(fm, m, next_sequence_id, fasta);
    }

    process_kmers(m);
    std::cout << "Kept " << kept_kmers.size() << " kmers\n";
    std::for_each(kept_kmers.begin(), kept_kmers.end(), compute_weight_of_signature);

    KmerGuts *kguts = write_hashtable(dir, kept_kmers);

    /*
     * Ick.
     */
    std::string fidx = dir + "/function.index";
    kguts->kmersH->function_array = kguts->load_functions(fidx.c_str(), &kguts->kmersH->function_count);
    kguts->kmersH->otu_array = kguts->load_otus("/dev/null", &kguts->kmersH->otu_count);

    /*

    for (auto file: all_fasta_data)
    {
	recall_fasta(fm, file, kguts);
    }

    */

    show_ps();
    rejected_stream.close();
    kept_stream.close();
    kept_function_stream.close();
    return 0;
}
