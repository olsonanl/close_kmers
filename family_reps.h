#ifndef _FAMILY_REPS_H
#define _FAMILY_REPS_H

/*
 * Family reps database.
 *
 * In-memory database mapping from local family ID
 * to a set of pegs that are the representative proteins
 * in that family. We also store location information for each protein.
 */

#include <boost/filesystem.hpp>
#include <string>

#include <map>
#include <vector>

#include <cstdint>

class FamilyReps
{
public:

    FamilyReps();
    void load_reps_file(const boost::filesystem::path &file);
    void load_reps_directory(const boost::filesystem::path &dir);
    
    struct RepData
    {
	std::string feature_id;
	std::string contig;
	uint32_t contig_length;
	uint32_t start;
	uint32_t end;
	char strand;
    };

    std::map<std::string, std::vector<RepData>> reps_;
};


#endif // _FAMILY_REPS_H
