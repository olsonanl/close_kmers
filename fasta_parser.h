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
    void parse_char(char c);
    void parse_complete();

private:
    state cur_state_;
    std::string cur_id_;
    std::string cur_seq_;
    
    std::function<int(const std::string &id, const std::string &seq)> on_seq_;
};


#endif
