#ifndef _threadpool_h
#define _threadpool_h

#include <memory>
#include <boost/asio.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/tss.hpp>
#include "kguts.h"

class ThreadPool : public std::enable_shared_from_this<ThreadPool>
{
public:

    ThreadPool(const std::string &kmer_dir);
    ~ThreadPool();
    
    void start(int n_threads);
    void stop();

    template <typename CompletionHandler>
	BOOST_ASIO_INITFN_RESULT_TYPE(CompletionHandler, void ())
	post(BOOST_ASIO_MOVE_ARG(CompletionHandler) handler) {
	io_service_.post(handler);
    };

    std::string kmer_dir_;
    std::unique_ptr<boost::asio::io_service::work> work_;

    KmerGuts::kmer_memory_image_t *image_;
    boost::thread_specific_ptr<KmerGuts> kguts_;
    boost::thread_specific_ptr<int> thread_index_;

    boost::asio::io_service io_service_;
    boost::thread_group thread_pool_;
};

#endif /* _threadpool_h */
