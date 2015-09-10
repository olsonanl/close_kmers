#ifndef _fasta_parser_h
#define _fasta_parser_h

#include <iostream>
#include <functional>

class FastaParser
{
public:
    FastaParser(std::function<int(const std::string &id, const std::string &seq)> on_seq);
    void parse(std::istream &stream);

private:
    std::function<int(const std::string &id, const std::string &seq)> on_seq_;
};


#endif
