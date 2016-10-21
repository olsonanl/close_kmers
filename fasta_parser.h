#ifndef _fasta_parser_h
#define _fasta_parser_h

#include <iostream>
#include <functional>

class FastaParser
{
public:

    enum state {
	s_start = 0,
	s_id,
	s_defline,
	s_data,
	s_id_or_data,
    };


    FastaParser(std::function<int(const std::string &id, const std::string &seq)> on_seq);
    void parse(std::istream &stream);
    void init_parse();
    inline void parse_char(char c)
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
		err = "Bad data character '";
		err += c;
		err += "'";
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
		err = "Bad id or data character '";
		err += c;
		err += "'";
	    }
	    break;
	}
	if (!err.empty())
	{
	    std::cerr << "Error found " << err << std::endl;
	    return;
	}
    }

    void parse_complete();

private:
    state cur_state_;
    std::string cur_id_;
    std::string cur_seq_;
    
    std::function<int(const std::string &id, const std::string &seq)> on_seq_;
};


#endif
