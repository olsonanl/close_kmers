#include "nudb_kmer_db.h"
#include "kmer_nudb.h"
#include <memory>

#include <stdexcept>

int main(int argc, char **argv)
{
    std::string path = "/dev/shm/out9.dat";

    auto cb = [](
	void const* key,        // A pointer to the item key
	std::size_t key_size,   // The size of the key (always the same)
	void const* data,       // A pointer to the item data
	std::size_t data_size,  // The size of the item data
	nudb::error_code& ec          // Indicates an error (out parameter)
	) {
	auto *kmer = static_cast<const NuDBKmerDb<8>::key_type *>(key);
	auto *kdata = static_cast<const NuDBKmerDb<8>::KData *>(data);
	double dev = std::sqrt(static_cast<double>(kdata->var));
	std::cout << *kmer << " "
		  << kdata->function_index << " " 
		  << kdata->n_proteins << " " 
		  << kdata->mean << " " 
		  << kdata->median << " " 
		  << dev << "\n" ;
    };
    nudb::error_code ec;
    nudb::visit(path, cb, [](auto, auto) {}, ec);
}

