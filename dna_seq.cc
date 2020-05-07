#include "dna_seq.h"
#include "trans_table.h"

#include <boost/algorithm/string.hpp>

#include <stdexcept>
#include <utility>

std::list<std::pair<int, std::list<std::string> > > DNASequence::get_possible_proteins(const TranslationTable &trans)
{
    std::list<std::pair<int, std::list<std::string> > > ret;

    for (int frame: { 1, 2, 3, -1, -2, -3 })
    {
	std::string p = get_translated_frame(trans, frame);
	std::list<std::string> l;
	boost::split(l, p, boost::is_any_of("*"), boost::token_compress_on);
	auto item = std::make_pair(frame, l);
	ret.push_back(item);
    }

    return ret;
}

std::string DNASequence::get_translated_frame(const TranslationTable &trans, int frame)
{
    const std::string &mseq = frame < 0 ? reverse_seq() : seq();

    if (frame < -3 || frame == 0 || frame > 3)
	throw std::runtime_error("Invalid frame " + std::to_string(frame) + " requested in get_translated_Frame");
    
    size_t offset = abs(frame) - 1;

    // std::cout << "for frame " << frame << ": " << mseq << "\n";
    auto beg = mseq.begin() + offset;
    return trans.translate(beg, mseq.end());
}

const std::string &DNASequence::reverse_seq()
{
    if (reverse_seq_.empty())
    {
	for (auto it = seq_.rbegin(); it != seq_.rend(); it++)
	    reverse_seq_.push_back(complement(*it));
    }
    return reverse_seq_;
}
