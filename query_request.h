#ifndef _query_request_h
#define _query_request_h

#include "compute_request.h"
#include "krequest2.h"
#include "kmer.h"
#include "fasta_parser.h"
#include <memory>

class QueryRequest : public ComputeRequest
{
public:
    QueryRequest(std::shared_ptr<KmerRequest2> owner, int content_length, bool chunked = false);
    ~QueryRequest() { std::cerr << "Destroy queryrequest " << this << "\n"; }
    void run();

private:
    int content_length_;
    bool chunked_;
    std::shared_ptr<KmerRequest2> owner_;
    FastaParser parser_;
    bool header_written_;

    typedef std::vector<std::pair<std::string, std::string> > work_list_t;
    std::shared_ptr<work_list_t> current_work_;

    int on_parsed_seq(const std::string &id, const std::string &seq);
    void on_data(boost::system::error_code err, size_t bytes);

};

#endif
