#ifndef _dna_seq_h
#define _dna_seq_h

#include <string>
#include <list>

class TranslationTable;

class DNASequence
{
public:
    DNASequence(const std::string &id, const std::string &seq) : id_(id), seq_(seq) {};

    const std::string &id() { return id_; }
    const std::string &seq() { return seq_; }
    const std::string &reverse_seq();

    std::string get_translated_frame(const TranslationTable &trans, int frame);
    std::list<std::pair<int, std::list<std::string> > > get_possible_proteins(const TranslationTable &trans);


private:
    const std::string &id_;
    const std::string &seq_;

    std::string reverse_seq_;

    std::string::value_type complement(std::string::value_type c) {
	switch (c)
	{
	case 'a':
	    return 't';
	case 'A':
	    return 'T';
	    
	case 'c':
	    return 'g';
	case 'C':
	    return 'G';
	    
	case 'g':
	    return 'c';
	case 'G':
	    return 'C';
	    
	case 't':
	case 'u':
	    return 'a';
	case 'T':
	case 'U':
	    return 'A';
	
	case 'm':
	    return 'k';
	case 'M':
	    return 'K';
	    
	case 'r':
	    return 'y';
	case 'R':
	    return 'Y';
	    
	case 'w':
	    return 'w';
	case 'W':
	    return 'W';
	    
	case 's':
	    return 'S';
	case 'S':
	    return 'S';
	    
	case 'y':
	    return 'r';
	case 'Y':
	    return 'R';
	    
	case 'k':
	    return 'm';
	case 'K':
	    return 'M';
	    
	case 'b':
	    return 'v';
	case 'B':
	    return 'V';
	    
	case 'd':
	    return 'h';
	case 'D':
	    return 'H';
	    
	case 'h':
	    return 'd';
	case 'H':
	    return 'D';
	    
	case 'v':
	    return 'b';
	case 'V':
	    return 'B';
	    
	case 'n':
	    return 'n';
	case 'N':
	    return 'N';
	    
	default:
	    return c;
	}
    }
};

#endif /* _dna_seq_h */
