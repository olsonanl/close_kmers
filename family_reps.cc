#include "family_reps.h"

#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <stdexcept>
#include <iostream>

using namespace boost::filesystem;

FamilyReps::FamilyReps()
{
}

void FamilyReps::load_reps_file(const path &file)
{
    ifstream f(file);

    if (f.fail())
    {
	throw std::runtime_error("load_reps_file " + file.string() + " failed");
    }

    std::string line;
    std::vector<std::string> cols;

    // skip headers
    std::getline(f, line);
    int line_number = 1;

    try {
	
	while (std::getline(f, line))
	{
	    line_number++;
	    boost::split(cols, line, boost::is_any_of("\t"));
	    //std::cout << "have " << cols[2] << "\n";

	    if (cols.size() < 10)
	    {
		std::cerr << "Short line " << line_number << " in " << file.string() << "\n";
		continue;
	    }
	    
	    reps_[cols[3]].emplace_back(RepData{ cols[2],
			cols[5], 
			cols[9].length() ? uint32_t(std::stoul(cols[9])) : 0,
			uint32_t(std::stoul(cols[6])),
			uint32_t(std::stoul(cols[7])),
			cols[8][0]});
	}
    }
    catch (std::exception &e)
    {
	std::cerr << "Error loading " << file.string() << " at line " << line_number << ": " << e.what() << "\n";
    }
	    
    /*
    for (auto x: reps_)
    {
	std::cout << x.first << ": \n";
	for (auto y: x.second)
	{
	    std::cout << y.feature_id << " " << y.strand << " " << y.start << "\n";
	}
    }
    */
}


void FamilyReps::load_reps_directory(const path &dir)
{
    if (!is_directory(dir))
	throw std::runtime_error("load_reps_directory: " + dir.string() + " is not a directory");

    for (auto di = directory_iterator(dir); di != directory_iterator(); di++)
    {
	std::cout << di->path().string() << "\n";
	load_reps_file(di->path());
    }
}

