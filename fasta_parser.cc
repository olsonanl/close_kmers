#include "fasta_parser.h"
#include <cctype>

FastaParser::FastaParser(std::function<int(const std::string &id, const std::string &seq)> on_seq) : on_seq_(on_seq)
{
    init_parse();
}

void FastaParser::init_parse()
{
    cur_state_ = s_start;
    cur_id_ = "";
    cur_seq_ = "";
}

void FastaParser::parse(std::istream &stream)
{
    char c;
    init_parse();
    while (stream.get(c))
    {
	parse_char(c);
    }
    on_seq_(cur_id_, cur_seq_);
}

void FastaParser::parse_complete()
{
    on_seq_(cur_id_, cur_seq_);
    cur_id_ = "";
    cur_seq_ = "";
}

	
