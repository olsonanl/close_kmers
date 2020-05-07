#include "trans_table.h"

#include <iostream>
#include <tuple>
#include <stdexcept>
#include <cassert>

static TranslationTable::raw_code table_11
{
        "    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	"  Starts = ---M------**--*----M------------MMMM---------------M------------",
	"  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
	"  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
	"  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
};

TranslationTable TranslationTable::make_table(int code)
{
    if (code == 11)
    {
	return TranslationTable(table_11);
    }
    else
    {
	throw std::runtime_error("invalid genetic code " + code);
    }
}

TranslationTable::TranslationTable(const raw_code &raw) :
    aa_table_(65),
    start_table_(65)
{
    parse_raw_code(raw);
}

void TranslationTable::parse_raw_code(const raw_code &raw)
{
    auto aa_pos = raw.aas.find(" = ") + 3;
    auto starts_pos = raw.starts.find(" = ") + 3;
    auto base1_pos = raw.base1.find(" = ") + 3;
    auto base2_pos = raw.base2.find(" = ") + 3;
    auto base3_pos = raw.base3.find(" = ") + 3;

    assert(aa_pos == starts_pos);
    assert(aa_pos == base1_pos);
    assert(aa_pos == base2_pos);
    assert(aa_pos == base3_pos);

    for (auto pos = aa_pos; pos < raw.aas.length(); pos++)
    {
	auto aa = raw.aas[pos];
	auto start = raw.starts[pos];
	auto base1 = raw.base1[pos];
	auto base2 = raw.base2[pos];
	auto base3 = raw.base3[pos];

	uint8_t offset = encode_triple(base1, base2, base3);
	aa_table_[offset] = aa;
	start_table_[offset] = start;
    }
    aa_table_[64] = 'X';
    start_table_[64] = '-';
}

std::string TranslationTable::translate(std::string::const_iterator beg, std::string::const_iterator end) const
{
    std::string ret;
    for (auto it = beg; it != end;)
    {
	if (std::distance(it,end) >= 3)
	{
	    auto c1 = *it++;
	    auto c2 = *it++;
	    auto c3 = *it++;
	    std::string::value_type aa;
	    uint8_t offset = encode_triple(c1, c2, c3);
	    aa = aa_table_[offset];
	    ret.push_back(aa);
	}
	else
	    break;
    }
    return ret;
}
