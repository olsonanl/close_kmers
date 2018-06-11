#ifndef _kmer_image_h
#define _kmer_image_h

#include <string>

#define KMER_VERSION 1

/*
 * num_sigs is the number of buckets in the hash table.
 */
typedef struct kmer_memory_image {
    unsigned long long num_sigs;
    unsigned long long entry_size;
    long long  version;
} kmer_memory_image_t;

typedef struct sig_kmer {
    unsigned long long  which_kmer;
    int  otu_index;
    unsigned short  avg_from_end;
    int  function_index;
    float function_wt;
} sig_kmer_t;

class KmerImage
{
public:

    KmerImage(const std::string &data_dir);

    void attach();
    void detach();
    kmer_memory_image_t *map_image_file(const std::string &data_dir, size_t &size);

    kmer_memory_image_t *image() { return image_; }
    kmer_memory_image_t *image_;
    size_t image_size_;
    std::string data_dir_;
};


#endif
