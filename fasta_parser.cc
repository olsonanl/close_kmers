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

void FastaParser::parse_char(char c)
{
    /*
    std::cout << "top " << c << " " << (0+cur_state_) << std::endl;
    std::cout << cur_id_ << "\n";
    std::cout << cur_seq_ << "\n";
    */
    std::string err;
    switch (cur_state_)
    {
    case s_start:
	if (c != '>')
	{
	    err = "Missing >";
	}
	else
	{
	    cur_state_ = s_id;
	}
	break;
	    
    case s_id:
	if (isblank(c))
	{
	    cur_state_ = s_defline;
	}
	else if (c == '\n')
	{
	    cur_state_ = s_data;
	}
	else
	{
	    cur_id_.push_back(c);
	}
	break;

    case s_defline:
	if (c == '\n')
	{
	    cur_state_ = s_data;
	}
	break;

    case s_data:
	if (c == '\n')
	{
	    cur_state_ = s_id_or_data;

	}
	else if (isalpha(c))
	{
	    cur_seq_.push_back(c);
	}
	else
	{
	    err = "Bad data character";
	}
	break;

    case s_id_or_data:
	if (c == '>')
	{
	    if (on_seq_)
		on_seq_(cur_id_, cur_seq_);
	    cur_id_ = "";
	    cur_seq_ = "";
	    cur_state_ = s_id;
	}
	else if (c == '\n')
	{
	    cur_state_ = s_id_or_data;
	}
	else if (isalpha(c))
	{
	    cur_seq_.push_back(c);
	    cur_state_ = s_data;
	}
	else
	{
	    err = "Bad data character";
	}
	break;
    }
    if (!err.empty())
    {
	std::cerr << "Error found " << err << std::endl;
	return;
    }
}
	
