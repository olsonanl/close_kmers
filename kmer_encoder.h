#ifndef _kmer_encoder_h
#define _kmer_encoder_h

#include <cstdint>
#include <vector>

#include "kmer_params.h"

class KmerEncoder
{
public:
    KmerEncoder();

    inline unsigned char to_amino_acid_off(uint8_t c) {
	return aa_to_offset_[c];
    };

    unsigned long long encoded_kmer(unsigned char *p) {
//    std::cout << "encode '" << p << "'\n";
	unsigned long long encodedK = *p;
	int i;
	for (i=1; (i <= KMER_SIZE-1); i++) {
	    encodedK = (encodedK * 20) + *(p+i);
	}

	if (encodedK > MAX_ENCODED) {
	    fprintf(stderr,"bad encoding - input must have included invalid characters\n");
	    for (i=0; (i < KMER_SIZE); i++) {
		fprintf(stderr,"%d ",*(p+i));
	    }
	    fprintf(stderr,"\n");
	    exit(2);
	}
	return encodedK;
    }

    unsigned long long encoded_aa_kmer(const char *p) {
	unsigned char aa_off[KMER_SIZE];
	int j;
	for (j=0; (j < KMER_SIZE); j++) {
	    char prot_c = *(p+j);
	    aa_off[j] = to_amino_acid_off(prot_c);
	    /*
	     * If we have an invalid char, return an invalid encoding.
	     */
	    if (aa_off[j] >= 20)
		return MAX_ENCODED + 1;
	}
	return encoded_kmer(aa_off);
    }

    template <class T, size_t N>
	unsigned long long encoded_aa_kmer(std::array<T, N> &p) {
	unsigned char aa_off[N+1];
	unsigned char *ap = aa_off;
	for (auto prot_c: p)
	{
	    *ap = to_amino_acid_off(prot_c);
	    /*
	     * If we have an invalid char, return an invalid encoding.
	     */
	    if (*ap >= 20)
		return MAX_ENCODED + 1;
	    ap++;
	}
	*ap = 0;
	return encoded_kmer(aa_off);
    }

    void decoded_kmer(unsigned long long encodedK,char *decoded) {
  
	int i;
	*(decoded+KMER_SIZE) = '\0';
	unsigned long long x = encodedK;
	
	for (i=KMER_SIZE-1; (i >= 0); i--) {
	    *(decoded+i) = prot_alpha[x % 20];
	    x = x / 20;
	}
    }

    uint8_t aa_to_offset_[256];
    static constexpr char prot_alpha[20] = {
	'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' 
    };

};


#endif // _kmer_encoder_h
