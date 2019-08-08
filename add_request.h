#ifndef _add_request_h
#define _add_request_h

#include "compute_request.h"
#include "krequest2.h"
#include "kmer.h"
#include "fasta_parser.h"
#include <memory>

class AddRequest : public ComputeRequest
{
public:
    AddRequest(std::shared_ptr<KmerRequest2> owner, std::shared_ptr<KmerPegMapping> mapping, size_t content_length, bool chunked = false);
    ~AddRequest() { std::cerr << "Destroy addrequest " << this << "\n"; }
    void run();

private:
    bool silent_;
    std::shared_ptr<KmerPegMapping> mapping_;
    size_t content_length_;
    bool chunked_;
    std::shared_ptr<KmerRequest2> owner_;
    FastaParser parser_;

    typedef std::vector<std::pair<std::string, std::string> > work_list_t;
    std::shared_ptr<work_list_t> current_work_;

    int on_parsed_seq(const std::string &id, const std::string &seq);
    void on_data(boost::system::error_code err, size_t bytes);

};

#endif
