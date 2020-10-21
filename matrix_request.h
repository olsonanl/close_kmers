#ifndef _matrix_request_h
#define _matrix_request_h

#include "compute_request.h"
#include "krequest2.h"
#include "kmer.h"
#include "fasta_parser.h"
#include <memory>
#include "prot_seq.h"

class MatrixRequest : public ComputeRequest
{
public:
    MatrixRequest(std::shared_ptr<KmerRequest2> owner, std::shared_ptr<KmerPegMapping> mapping, size_t content_length, bool chunked = false);
    ~MatrixRequest() { std::cerr << "Destroy matrixrequest " << this << "\n"; }
    void run();

private:
    std::shared_ptr<KmerPegMapping> mapping_;
    size_t content_length_;
    bool chunked_;
    std::shared_ptr<KmerRequest2> owner_;
    FastaParser parser_;

    std::map<KmerPegMapping::encoded_id_t, size_t> matrix_proteins_;
    std::map<std::pair<KmerPegMapping::encoded_id_t, KmerPegMapping::encoded_id_t>, unsigned long> distance_;

    typedef std::vector<ProteinSequence> work_list_t;
    std::shared_ptr<work_list_t> current_work_;

    int on_parsed_seq(const std::string &id, const std::string &seq);

    void on_hit(KmerPegMapping::encoded_id_t id, const KmerGuts::hit_in_sequence_t &kmer);
    void on_data(boost::system::error_code err, size_t bytes);

    void process_results();
};

#endif
