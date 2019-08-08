
#include "kguts.h"
#include "kmer_encoder.h"

constexpr char KmerEncoder::prot_alpha[20];

KmerEncoder::KmerEncoder()
{
    for (int i = 0; i < 255; i++)
	aa_to_offset_[i] = 20;
    for (int i = 0; i < 20; i++)
	aa_to_offset_[(int) prot_alpha[i]] = i;
}
