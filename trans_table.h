#ifndef _trans_table_h
#define _trans_table_h

/*
    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M------**--*----M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
*/


#include <string>
#include <cstdint>
#include <iterator> 
#include <vector>
#include <iostream>

class TranslationTable
{
public:

    struct raw_code
    {
	std::string aas;
	std::string starts;
	std::string base1;
	std::string base2;
	std::string base3;
    };

    static TranslationTable make_table(int code);
    TranslationTable(const raw_code &raw);
//    TranslationTable(const TranslationTable &) {}

    std::string translate(std::string::const_iterator beg, std::string::const_iterator end) const;

private:

    std::vector<std::string::value_type> aa_table_;
    std::vector<std::string::value_type> start_table_;

    void parse_raw_code(const raw_code &raw);

    uint8_t encode_char(std::string::value_type c) const
    {
	switch (c)
	{
	case 'a':
	case 'A':
	    return 0;
	
	case 'c':
	case 'C':
	    return 1;
	
	case 'g':
	case 'G':
	    return 2;
	
	case 't':
	case 'u':
	case 'T':
	case 'U':
	    return 3;
	
	default:
	    return 4;
	}
    }

    uint8_t encode_triple(std::string::value_type c1, std::string::value_type c2, std::string::value_type c3) const
    {
	uint8_t e1 = encode_char(c1);
	uint8_t e2 = encode_char(c2);
	uint8_t e3 = encode_char(c3);
	uint8_t offset = 64;

	if (e1 < 4 && e2 < 4 && e3 < 4)
	    offset = e1 * 16 + e2 * 4 + e3;

	return offset;
    }


};

#endif /* _trans_table_h */
