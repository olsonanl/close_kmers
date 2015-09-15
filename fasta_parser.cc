#include "fasta_parser.h"
#include <cctype>

FastaParser::FastaParser(std::function<int(const std::string &id, const std::string &seq)> on_seq) : on_seq_(on_seq)
{
}

void FastaParser::parse(std::istream &stream)
{
    enum state {
	s_start,
	s_id,
	s_defline,
	s_data,
	s_id_or_data,
    };
    
    char c;
    state s = s_start;
    std::string cur_id;
    std::string cur_seq;
    while (stream.get(c))
    {
//	std::cout << "top " << c << " " << s << std::endl;
	std::string err;
	switch (s)
	{
	case s_start:
	    if (c != '>')
	    {
		err = "Missing >";
	    }
	    else
	    {
		s = s_id;
	    }
	    break;
	    
	case s_id:
	    if (isblank(c))
	    {
		s = s_defline;
	    }
	    else if (c == '\n')
	    {
		s = s_data;
	    }
	    else
	    {
		cur_id.push_back(c);
	    }
	    break;

	case s_defline:
	    if (c == '\n')
	    {
		s = s_data;
	    }
	    break;

	case s_data:
	    if (c == '\n')
	    {
		s = s_id_or_data;
	    }
	    else if (isalpha(c))
	    {
		cur_seq.push_back(c);
	    }
	    else
	    {
		err = "Bad data character";
	    }
	    break;

	case s_id_or_data:
	    if (c == '>')
	    {
		on_seq_(cur_id, cur_seq);
		cur_id = "";
		cur_seq = "";
		s = s_id;
	    }
	    else if (c == '\n')
	    {
		s = s_id_or_data;
	    }
	    else if (isalpha(c))
	    {
		cur_seq.push_back(c);
		s = s_data;
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
    on_seq_(cur_id, cur_seq);
}
	
