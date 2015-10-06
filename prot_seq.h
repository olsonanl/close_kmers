#ifndef _prot_seq_h
#define _prot_seq_h

#include <string>

class ProteinSequence
{
public:
ProteinSequence(const std::string &id, const std::string &seq) : id_(id), seq_(seq) {};

    const std::string &id() { return id_; }
    const std::string &seq() { return seq_; }
private:
    const std::string id_;
    const std::string seq_;
};

#endif /* _prot_seq_h */
