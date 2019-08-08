#ifndef _fq_process_request_h
#define _fq_process_request_h

#include "compute_request.h"
#include "krequest2.h"
#include "kmer.h"
#include "fastq_parser.h"
#include <memory>
#include "prot_seq.h"
#include "zlib_support.h"
#include "dna_seq.h"
#include "trans_table.h"
#include "family_mapper.h"
#include <iostream>

class FqProcessRequest : public ComputeRequest
{
public:
    FqProcessRequest(std::shared_ptr<KmerRequest2> owner, std::shared_ptr<KmerPegMapping> mapping, bool family_mode, size_t content_length, bool chunked = false);
    ~FqProcessRequest() { std::cerr << "Destroy fq_processrequest " << this << "\n"; }
    void run();

private:
    int on_parsed_seq(const std::string &id, const std::string &seq,  FamilyMapper &mapper, std::ostream &os);
    void on_first_data(boost::system::error_code err, size_t bytes);
    void on_compressed_data(boost::system::error_code err, size_t bytes);
    void on_uncompressed_data(boost::system::error_code err, size_t bytes);
    void process_compressed_block();
    void process_block(bool finished);
    void process_decompressed_block(std::shared_ptr<boost::asio::streambuf> buf, zlib_support::ReturnStatus status);

    void process_data(boost::asio::streambuf &buf, bool finished, std::function<void()> cb);
    void finish_data();


    bool family_mode_;
    std::shared_ptr<KmerPegMapping> mapping_;
    size_t content_length_;
    bool chunked_;
    std::shared_ptr<KmerRequest2> owner_;
    FastqParser parser_;

    std::shared_ptr<zlib_support::GzipDecoder> gzip_decoder_;

    TranslationTable trans_table_;
    bool header_written_;
};

#endif
