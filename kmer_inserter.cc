#include "kmer_inserter.h"

KmerInserter::KmerInserter(int n_workers, std::shared_ptr<KmerPegMapping> mapping) :
    n_workers_(n_workers),
    mapping_(mapping)
{
    for (int i = 0; i < n_workers; i++)
    {
	queues_.emplace_back();
    }
}

void KmerInserter::start()
{
    for (int i = 0; i < n_workers(); i++)
    {
	thread_pool_.create_thread([this, i]() {
		thread_main(i);
	    });
    }
}

void KmerInserter::stop()
{
    for (int i = 0; i < n_workers(); i++)
    {
	std::cerr << "Push end work to " << i << "\n";
	std::shared_ptr<WorkElement> endwork = std::make_shared<WorkElement>();
	endwork->shutdown = true;
	push_work(i, endwork);
    }
    std::cerr << "Await KI threads\n";
    thread_pool_.join_all();
}

void KmerInserter::thread_main(int modulus)
{
    WorkQueue &my_queue = queues_[modulus];

    while (1)
    {
	std::shared_ptr<WorkElement> work;
	my_queue.pop(work);
	if (work->shutdown)
	{
	    std::cerr << "Shutting down KmerInserter " << modulus << "\n";
	    break;
	}
//	boost::lock_guard<boost::mutex> guard(mtx_);

//	std::cerr << "KI " << modulus << " process work item of length " << work->work.size() << "\n";
	for (auto item: work->work)
	{
//	    std::cerr << "   " << item.second << " " << item.first << "\n";
	    mapping_->add_fam_mapping(item.second, item.first);
	}
    }
}
