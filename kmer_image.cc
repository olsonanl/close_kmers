#include "kmer_image.h"

#include <iostream>
#include <sys/file.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>

#include "global.h"
#include <boost/program_options.hpp>

KmerImage::KmerImage(const std::string &data_dir)
    : image_(0)
    , image_size_(0)
    , data_dir_(data_dir)
{
    attach();
}

void KmerImage::attach()
{
    if (image_)
	detach();

    image_ = map_image_file(data_dir_, image_size_);
}

void KmerImage::detach()
{
    if (image_)
    {
	munmap(image_, image_size_);
	image_ = 0;
    }
}


kmer_memory_image_t *KmerImage::map_image_file(const std::string &data_dir, size_t &size)
{
    std::string fileM = data_dir + "/kmer.table.mem_map";
    int fd;
    std::cout << "mmap " << fileM << "\n";
    if ((fd = open(fileM.c_str(), O_RDONLY)) < 0) {
	perror("open");
	exit(1);
    }

    /*
     * Set up for creating memory image from file. Start by determining file size
     * on disk with a stat() call.
     */
    struct stat sbuf;
    if (stat(fileM.c_str(), &sbuf) == -1) {
	fprintf(stderr, "stat %s failed: %s\n", fileM.c_str(), strerror(errno));
	exit(1);
    }
    unsigned long long file_size = (unsigned long long) sbuf.st_size;
    size = file_size;
    
    /* 
     * Memory map.
     */
    int flags = MAP_SHARED;
    
#ifdef MAP_POPULATE
    bool do_populate = false;
    if (g_parameters->count("no-populate-mmap"))
    {
	do_populate = ! ( (*g_parameters)["no-populate-mmap"].as<bool>() );
    }

    if (do_populate)
	flags |= MAP_POPULATE;
#endif
    
    kmer_memory_image_t *image = (kmer_memory_image_t *) mmap((caddr_t)0, file_size, PROT_READ, flags, fd, 0);

    if (image == (kmer_memory_image_t *)(-1)) {
	fprintf(stderr, "mmap of kmer_table %s failed: %s\n", fileM.c_str(), strerror(errno));
	exit(1);
    }

    /* Validate overall file size vs the entry size and number of entries */
    if (file_size != ((sizeof(sig_kmer_t) * image->num_sigs) + sizeof(kmer_memory_image_t))) {
	fprintf(stderr, "Version mismatch for file %s: file size does not match\n", fileM.c_str());
	exit(1);
    }
    
    /* 
     * Our image is mapped. Validate against the current version of this code.
     */
    if (image->version != (long long) KMER_VERSION) {
	fprintf(stderr, "Version mismatch for file %s: file has %lld code has %lld\n", 
		fileM.c_str(), image->version, (long long) KMER_VERSION);
	exit(1);
    }
    
    if (image->entry_size != (unsigned long long) sizeof(sig_kmer_t)) {
	fprintf(stderr, "Version mismatch for file %s: file has entry size %lld code has %lld\n",
		fileM.c_str(), image->entry_size, (unsigned long long) sizeof(sig_kmer_t));
	exit(1);
    }
    fprintf(stderr, "Set size_hash=%lld from file size %lld\n", image->num_sigs, file_size);
    return image;
}

