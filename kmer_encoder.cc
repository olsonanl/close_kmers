
#include "kguts.h"
#include "kmer_encoder.h"

constexpr char KmerEncoder::prot_alpha[20];
constexpr char KmerEncoder::genetic_code[64];

KmerEncoder::KmerEncoder()
{
    for (int i = 0; i < 255; i++)
	aa_to_offset_[i] = 20;
    for (unsigned char i = 0; i < 20; i++)
	aa_to_offset_[(int) prot_alpha[i]] = i;
}
