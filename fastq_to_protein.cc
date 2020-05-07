#include "fastq_parser.h"
#include "trans_table.h"
#include "dna_seq.h"

#include <boost/program_options.hpp>
#include <algorithm>
#include <cmath>

#include <fstream>
#include <iomanip>

namespace po = boost::program_options;

void to_proteins(const std::string &fastq_file, std::ostream &output)
{
    std::ifstream fastq(fastq_file);

    std::string parse_error;
    int error_line;

    TranslationTable trans_table(TranslationTable::make_table(11));

    FastqParser parser;
    parser.set_callback([&trans_table, &output](const std::string &id, const std::string &seq) {
	    if (!id.empty())
	    {
		DNASequence dna(id, seq);
		// std::cout << seq << "\n";
		auto prots = dna.get_possible_proteins(trans_table);
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
			    output << ">" << id << ":" << frame << ":" << i << "\n" << prot << "\n";
			}
		    }
		}
	    }
	    return 0;
	});

    parser.set_error_callback([&parse_error, &error_line](const std::string err, int line, const std::string id) {
	    parse_error = err;
	    error_line = line;
	    return false;
	});
    
    parser.init_parse();
    parser.parse(fastq);
    parser.parse_complete();
}

int main(int argc, char* argv[])
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options] fastq-file\nAllowed options";
    po::options_description desc(x.str());

    std::string fastq_file;
    std::string output_file;
    
    desc.add_options()
	("output-file,o", po::value<std::string>(&output_file), "output file")
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

    std::ostream *out = &std::cout;
    std::ofstream file_out;
    if (!output_file.empty())
    {
	file_out.open(output_file);
	if (!file_out)
	{
	    std::cerr << "Error opening " << output_file << "\n";
	    exit(1);
	}
	out = &file_out;
    }

    try {
	to_proteins(fastq_file, *out);
    }
    catch (std::exception &e)
    {
	std::cerr << "Caught exception: " << e.what() << "\n";
	return -1;
    }
}
