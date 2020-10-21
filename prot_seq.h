#ifndef _prot_seq_h
#define _prot_seq_h

#include <string>
#include <cstring>

class ProteinSequence
{
public:
ProteinSequence(const std::string &id, const std::string &seq) : id_(id), seq_(seq) {};

    const std::string &id() { return id_; }
    const std::string &seq() { return seq_; }

    template <int N>
    struct KmerIterator
    {
	const char *ptr_;
	char buf[N+1];
	typedef KmerIterator self_type;
	
	KmerIterator(const char *ptr)  : ptr_(ptr) {};
	
	self_type operator++() { self_type i = *this; ptr_++; return i; }
	self_type operator++(int junk) { ptr_++; return *this; }
	const char * operator*() {
	    strncpy(buf, ptr_, N);
	    buf[N] = 0;
	    return buf;
	}
	bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
	bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }
	
	
    };

    template <int N>
    KmerIterator<N> kmer_begin() { return KmerIterator<N>(seq_.c_str()); }
    template <int N>
    KmerIterator<N> kmer_end() { return KmerIterator<N>(seq_.c_str() + seq_.size() - N + 1); }
    
   
private:
    const std::string id_;
    const std::string seq_;
};

#endif /* _prot_seq_h */
