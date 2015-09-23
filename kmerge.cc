/*
 * Merge Jim's kmer files which are of the form
 * kmer 0/1
 *
 * into a single large matrix. The columns (and input files) are defined
 * by ...
 */

#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/functional/hash.hpp>
#include <cmath>
#include "kmc_file.h"
#include "kmer_api.h"
#include "stringutil.h"

using namespace boost::filesystem;
namespace po = boost::program_options;

class Kmer
{
public:

    Kmer(const std::string &k)
	{
	    if (k.size() > 31)
	    {
		std::cerr << "invalid kmer size " << k.size() << "\n";
		exit(1);
	    }
	    strncpy(kmer_, k.c_str(), k.size());
	}

    Kmer()
	{
	    kmer_[0] = 0;
	}

    bool empty() {
	return kmer_[0] == 0;
    }
	    
    char kmer_[32];
};

inline bool operator==(const Kmer& lhs, const Kmer& rhs)
{
    return strcmp(lhs.kmer_, rhs.kmer_) == 0;
}
inline bool operator!=(const Kmer& lhs, const Kmer& rhs){return !(lhs == rhs);}

std::ostream& operator<<(std::ostream& out, const Kmer& k){
    return out << k.kmer_;
}

namespace std {
    template <> struct hash<Kmer>
    {
	size_t operator()(const Kmer &k) const
	    {
		return boost::hash_range(k.kmer_, k.kmer_ + 32);
	    }
    };
};

class KmerSet
{
public:
    typedef Kmer kmer_t;
    // typedef std::string kmer_t;
    
    typedef std::unordered_map<kmer_t, std::vector<bool> > map_t;

    map_t kmer_map_;
    std::map<std::string, int> file_to_column_;
    std::vector<bool> default_value_;

    std::map<std::vector<bool>, std::vector<kmer_t>> pattern_seen_;

    KmerSet();
    
    void add_files(const std::vector<std::string> &files, bool invert);

    void process_files(const std::vector<std::string> &files, bool invert);
    void process_file(const std::string &file, int idx, bool invert);
    void process_kmc_file(const std::string &file, int idx, bool invert);
    void add_kmer(int idx, const kmer_t &kmer, bool bval);

    void remove_duplicate_values();

    void dump(std::ostream &stream);
};

class Adaboost
{
public:
    Adaboost(KmerSet &kmers);

    void compute(int n_rounds);
    double compute_error(std::vector<bool> &kvec, std::vector<double> &pvec);
    std::vector<double> update_prob_table(double alpha, std::vector<double> &old_prob, std::vector<bool> &kvec);

    KmerSet &kmers_;
};

int main(int argc, char **argv)
{
    std::ostringstream x;

    std::string kmer_dir;
    std::string output_file;
    std::string sus_file;
    std::string res_file;
    int rounds;
    int max_files;
    x << "Usage: " << argv[0] << " [options] resistant-file susceptible-file\nAllowed options";
    po::options_description desc(x.str());
    desc.add_options()
	("help,h", "show this help message")
	("rounds,r", po::value<int>(&rounds)->default_value(10), "Number of rounds of Adaboost to run")
	("max-files", po::value<int>(&max_files)->default_value(-1), "Max number of files to process per data set")
	
	("kmer-dir,d", po::value<std::string>(&kmer_dir)->default_value("KMERS"), "Directory in which kmer files are found")
	("output-file,o", po::value<std::string>(&output_file), "Write output to this file instead of stdout")
	("resistant-file", po::value<std::string>(&res_file), "File containing list of kmer files for resistant genomes")
	("susceptible-file", po::value<std::string>(&sus_file), "File containing list of kmer files for susceptible genomes")
	;

    po::positional_options_description pd;
    pd.add("resistant-file", 1)
	.add("susceptible-file", 1)
	;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(pd).run(), vm);

    if (vm.count("help"))
    {
	std::cerr << desc << "\n";
	return 1;
    }

    try {
	po::notify(vm);
    } catch (po::required_option &e)
    {
	std::cerr << "Invalid command line: " << e.what() << "\n";
	std::cerr << desc << "\n";
	exit(1);
    }

    KmerSet kset;

    std::ifstream sus(sus_file);
    std::ifstream res(res_file);

    std::string line;

    std::vector<std::string> sus_files, res_files;
    int i = 0;
    while (std::getline(sus, line))
    {
	if (max_files >= 0 && i++ >= max_files)
	    break;
	
	if (line[0] == '/')
	    sus_files.push_back(line);
	else
	    sus_files.push_back(kmer_dir + "/" + line);
    }
    i = 0;
    while (std::getline(res, line))
    {
	if (max_files >= 0 && i++ >= max_files)
	    break;
	if (line[0] == '/')
	    res_files.push_back(line);
	else
	    res_files.push_back(kmer_dir + "/" + line);
    }

    kset.add_files(res_files, 0);
    kset.add_files(sus_files, 1);

    kset.process_files(res_files, 0);
    kset.process_files(sus_files, 1);

    /*
    
    if (vm.count("output-file"))
    {
	std::ofstream out(output_file);
	kset.dump(out);
    }
    else
    {
	kset.dump(std::cout);
    }

    */
    kset.remove_duplicate_values();
    Adaboost ada(kset);
    ada.compute(rounds);
}
    
KmerSet::KmerSet()
{
    kmer_map_.reserve(10000000);
}

void KmerSet::add_files(const std::vector<std::string> &files, bool invert)
{
    for (auto file: files)
    {
	std::cerr << "file " << file << "\n";
	if (file_to_column_.find(file) != file_to_column_.end())
	{
	    std::cerr << "error: file " << file << " is repeated\n";
	    exit(1);
	}
	int idx = file_to_column_.size();
	file_to_column_[file] = idx;
	default_value_.push_back(invert);
    }

}

void KmerSet::process_files(const std::vector<std::string> &files, bool invert)
{
    std::cerr << "Processing " << files.size() << " files\n";
    for (auto file: files)
    {
	int idx = file_to_column_[file];

	std::string kmc_pre = file + ".kmc_pre";
	if (endswith(file, ".kmc_pre") ||
	    endswith(file, ".kmc_suf"))
	{
	    process_kmc_file(file.substr(0, file.length() - 8), idx, invert);
	}
	else if (boost::filesystem::is_regular_file(kmc_pre))
	{
	    process_kmc_file(file, idx, invert);
	}
	else
	{
	    process_file(file, idx, invert);
	}
    }
}

void KmerSet::process_kmc_file(const std::string &file_base, int idx, bool invert)
{
    CKMCFile kmc_file;

    if (!kmc_file.OpenForListing(file_base))
    {
	std::cerr << "could not open kmc file " << file_base << "\n";
	exit(1);
    }

    uint32 kmer_len, mode, counter_size, lut_prefix_len, signature_len, min_count;
    uint64 max_count, total_kmers;
    
    kmc_file.Info(kmer_len, mode, counter_size, lut_prefix_len, signature_len, min_count, max_count, total_kmers);

    std::cerr << "Processing kmer file #" << idx << " " << file_base << " kmer_size: " << kmer_len << " kmer count: " << total_kmers << "\n";
    
    CKmerAPI kmer(kmer_len);
    uint32 count;
    while (kmc_file.ReadNextKmer(kmer, count))
    {
	// std::cerr << "read " << kmer.to_string() << " " << count << "\n";
	bool bval = count ? 1 : 0;
	if (invert)
	    bval = !bval;
	add_kmer(idx, kmer.to_string(), bval);
    }
}

void KmerSet::process_file(const std::string &file, int idx, bool invert)
{
    std::ifstream str(file);
    
    std::cerr << "Processing file #" << idx << " " << file << "\n";
    std::string line;
    int line_num = 0;
    while (std::getline(str, line))
    {
	line_num++;
	std::size_t pos = line.find("\t");
	if (pos == std::string::npos)
	{
	    std::cerr << "Missing tab in " << file << " line " << line_num << "\n";
	    exit(1);
	}
	kmer_t kmer(line.substr(0, pos));
	int val = std::stoi(line.substr(pos + 1));

	bool bval = val ? 1 : 0;
	if (invert)
	    bval = !bval;
	// std::cerr << kmer << "\t" << val << "\t" << bval << "\n";

	add_kmer(idx, kmer, bval);
	#if 0
	if (line_num > 100)
	    break;
	#endif
    }
}

void KmerSet::add_kmer(int idx, const kmer_t &kmer, bool bval)
{
    auto it = kmer_map_.find(kmer);
    if (it == kmer_map_.end())
    {
	auto n = kmer_map_.emplace(std::make_pair(kmer, default_value_));
	n.first->second[idx] = bval;
    }
    else
    {
	it->second[idx] = bval;
    }
}

void KmerSet::remove_duplicate_values()
{
    std::cerr << "Removing duplicate matrix patterns. Original matrix size " << kmer_map_.size() << "\n";
    int ndups = 0;
    for (auto it = kmer_map_.begin(); it != kmer_map_.end(); )
    {
	auto seen_val = pattern_seen_.find(it->second);

	if (seen_val == pattern_seen_.end())
	{
	    auto n = pattern_seen_.emplace(std::make_pair(it->second, std::vector<kmer_t>()));
	    n.first->second.push_back(it->first);
	    ++it;
	}
	else
	{
	    seen_val->second.push_back(it->first);
	    // std::cerr << "dup " << it->first << " to " << seen_val->second[0] << "\n";
	    kmer_map_.erase(it++);
	    ndups++;
	}
    }
    std::cerr << "Complete. Final matrix size " << kmer_map_.size() << " after removing " << ndups << " duplicatess\n";
}

void KmerSet::dump(std::ostream &stream)
{
    for (auto ent: kmer_map_)
    {
	stream << ent.first;
	for (auto v: ent.second)
	{
	    stream << "\t" << v;
	}
	stream << "\n";
    }
}

Adaboost::Adaboost(KmerSet &kmers) : kmers_(kmers)
{
}

void Adaboost::compute(int n_rounds)
{
    int n = kmers_.default_value_.size();

    std::vector<double> prob(n, 1.0 / n);

    for (int round = 0; round < n_rounds; round++)
    {
	KmerSet::kmer_t bestk;
	double alpha;
	double error_min = 1.0;

	for (auto kent: kmers_.kmer_map_)
	{
	    double error = compute_error(kent.second, prob);
	    // std::cerr << "computed error " << error << "\n";
	    if (error < error_min)
	    {
		error_min = error;
		alpha = fabs((0.5 * log((1 - error) / error)));
		// std::cout << "computed alpha " << alpha << " from " << error << "\n";
		bestk = kent.first;
	    }
	}

	if (bestk.empty())
	{
	    std::cerr << "No bestk found at round " << round << "\n";
	    exit(1);
	}

	auto best_ent = kmers_.kmer_map_.find(bestk);
	if (best_ent == kmers_.kmer_map_.end())
	{
	    std::cerr << "Could not find entry for " << bestk << "\n";
	    exit(1);
	}

	std::cout << error_min << "\t" << alpha;
	auto seen_val = kmers_.pattern_seen_.find(best_ent->second);
	for (auto it: seen_val->second)
	    std::cout << "\t" << it;
	std::cout << "\n";
	
	prob = update_prob_table(alpha, prob, best_ent->second);
	kmers_.kmer_map_.erase(best_ent);
    }
}

std::vector<double> Adaboost::update_prob_table(double alpha,
						std::vector<double> &old_prob,
						std::vector<bool> &kvec)
{
    double Wc = exp(-alpha);
    double Wi = exp(alpha);

    double Z;
    std::vector<double> new_prob;

    auto kent = kvec.begin();
    auto pent = old_prob.begin();
    while (pent != old_prob.end())
    {
	double unprob;
	if (*kent)
	{
	    unprob = Wc * *pent;
	}
	else
	{
	    unprob = Wi * *pent;
	}
	new_prob.push_back(unprob);
	Z += unprob;
	    
	kent++;
	pent++;
    }

    for (auto it = new_prob.begin(); it != new_prob.end(); it++)
    {
	*it /= Z;
    }
    return new_prob;
}

double Adaboost::compute_error(std::vector<bool> &kvec, std::vector<double> &pvec)
{
    double error = 0.0;
    
    auto kent = kvec.begin();
    auto pent = pvec.begin();
    while (pent != pvec.end())
    {
	if (!*kent)
	{
	    error += *pent;
	}
	kent++;
	pent++;
    }
    return error;
}
