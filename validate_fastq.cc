#include "fastq_parser.h"
#include "trans_table.h"
#include "dna_seq.h"

#include <boost/program_options.hpp>
#include <algorithm>
#include <cmath>

#include <fstream>
#include <iomanip>

namespace po = boost::program_options;

void validate(std::istream &fastq)
{

    std::map<std::string, int> id_count;
    std::vector<unsigned long> sizes;
    unsigned long total = 0;

    bool valid = true;
    std::string parse_error;
    int error_line;

    TranslationTable t(TranslationTable::make_table(11));

    FastqParser parser;
    parser.set_callback([&t, &id_count, &total, &sizes](const std::string &id, const std::string &seq) {
	    if (!id.empty())
	    {
		id_count[id]++;
		total += seq.size();
		sizes.push_back(seq.size());

		#if 0

		DNASequence dna(id, seq);
		// std::cout << seq << "\n";
		auto prots = dna.get_possible_proteins(t);
		for (auto ent: prots)
		{
		    int frame;
		    std::list<std::string> proteins;
		    std::tie(frame, proteins) = ent;
		    int i = 0;
		    for (auto prot: proteins)
		    {
			i++;
			if (prot.length() > 10)
			{
			    std::cout << ">" << id << ":" << frame << ":" << i << "\n" << prot << "\n";
			}
		    }
		}
		#endif
	    }
	    return 0;
	});

    parser.set_error_callback([&valid, &parse_error, &error_line](const std::string err, int line, const std::string id) {
	    parse_error = err;
	    error_line = line;
	    valid = false;
	    return false;
	});
    
    parser.init_parse();
    parser.parse(fastq);
    parser.parse_complete();

    if (valid)
    {
	std::cout << "valid\t1\n";
	std::cout << "n_seqs\t" << sizes.size() << "\n";

	if (!sizes.empty())
	{
	    double n = (double) sizes.size();
	    double mean = (double) total / n;
	    
	    double stddev = 0.0;
	    if (sizes.size() > 1)
	    {
		double accum = 0.0;
		std::for_each (std::begin(sizes), std::end(sizes), [&accum, mean](const int i) {
			double d = (double) i;
			accum += (d - mean) * (d - mean);
		    });
		
		stddev = sqrt(accum / (n-1.0));
	    }
	    
	    std::cout << std::fixed << std::setprecision(2);

	    std::cout << "total_size\t" << total << "\n";
	    std::cout << "mean\t" << mean << "\n";
	    std::cout << "stddev\t" << stddev << "\n";
	}
    }
    else
    {
	std::cout << "valid\t0\n";
	std::cout << "n_seqs\t" << sizes.size() << "\n";
	std::cout << "error_message\t" << parse_error << "\n";
	std::cout << "error_line\t" << error_line << "\n";
    }
}

int main(int argc, char* argv[])
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options] fastq-file\nAllowed options";
    po::options_description desc(x.str());

    std::string fastq_file;
    
    desc.add_options()
	("fastq-file",  po::value<std::string>(&fastq_file)->required(), "fastq file input")
	("help,h", "show this help message")
	;
    po::positional_options_description pd;
    pd.add("fastq-file", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(pd).run(), vm);

    if (vm.count("help"))
    {
	std::cout << desc << "\n";
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

    try {
	if (fastq_file == "-")
	{
	    validate(std::cin);
	}
	else
	{
	    std::ifstream fastq(fastq_file);

	    validate(fastq);
	}
    }
    catch (std::exception &e)
    {
	std::cerr << "Caught exception: " << e.what() << "\n";
	return -1;
    }
}
