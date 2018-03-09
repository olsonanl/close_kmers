#include "family_reps.h"
#include <stdexcept>

using namespace boost::filesystem;

int main(int argc, char **argv)
{
    if (argc != 2)
	throw std::runtime_error("bad args");

    FamilyReps reps;

    path p(argv[1]);

    if (is_directory(p))
    {
	reps.load_reps_directory(p);
    }
    else
    {
	reps.load_reps_file(p);
    }
	
}
