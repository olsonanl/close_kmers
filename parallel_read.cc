#include <boost/thread/thread.hpp>

template <class Callback>
int read_block(const std::string &file, size_t block, size_t chunk_size, Callback cb);


template <class Callback>
void parallel_read(const std::string &file, int n_threads, Callback cb)
{
    std::ifstream is;
    is.open(file);
    if (is.fail())
    {
	std::cerr << "failed opening\n";
	exit(1);
    }

    is.seekg(0, is.end);
    size_t fsize = is.tellg();

    size_t chunk_size = (size_t) ceil((float) fsize / (float) (n_threads - 1));
    std::cerr << "size=" << fsize << " chunk=" << chunk_size << "\n";

    boost::thread_group pool;

    for (size_t block = 0; block < (size_t) n_threads; block++)
    {
	pool.create_thread([block, file, chunk_size, cb]() {
		read_block(file, block, chunk_size, cb);
	    });
    }
    pool.join_all();
}

template <class Callback>
int read_block(const std::string &file, size_t block, size_t chunk_size, Callback cb)
{
    size_t start = block * chunk_size;
    size_t end = start + chunk_size - 1;
    size_t nlines = 0;

    std::ifstream is;
    is.open(file);
    if (is.fail())
    {
	std::cerr << "failed opening\n";
	exit(1);
    }

    is.seekg(start, is.beg);
    
    std::cerr << "start block " << block << " start=" << start << " end=" << end << "\n";
    
    std::string line;

    size_t where = start;
    
    if (block > 0)
    {
	std::getline(is, line);
	where += line.size() + 1;
	//std::cerr << "discard='" << line << "'\n";
    }
    
    if (where > (end + 1))
    {
	//std::cerr << "finish after discard!\n";
	return 0;
    }
    
    while (std::getline(is, line))
    {
	where += line.size() + 1;
	cb(line);
	//std::cerr << "Read '" << line << "' now at " << where << "\n";
	nlines++;
	if (is.eof() || where > (end+1))
	{
	    //std::cerr << "finishing block at " << where << "\n";
	    break;
	}
    }
    std::cerr << "block " << block << " " << nlines << " lines\n";
    return (int) nlines;
}
