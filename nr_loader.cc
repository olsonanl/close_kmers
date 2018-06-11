#include "nr_loader.h"

#include <boost/filesystem.hpp>

#include "kmer.h"
#include "threadpool.h"
#include "fasta_parser.h"

/*
 * Preload a families NR file.
 *
 * We do this by parsing out the file into chunks and feeding those chunks to
 * the threadpool to parse and insert.
 */

NRLoader::NRLoader(NRLoadState &load_state, const std::string &file,
		   std::shared_ptr<KmerPegMapping> root_mapping,
		   std::shared_ptr<ThreadPool> thread_pool,
		   size_t n_files, bool family_mode,
		   KmerInserter &inserter
    ) :
    family_mode_(family_mode),
    inserter_(inserter),
    load_state_(load_state),
    my_load_state_(file),
    file_(file),
    root_mapping_(root_mapping),
    thread_pool_(thread_pool),
    n_files_(n_files),
    chunks_started_(0),
    chunks_finished_(0)
{
    std::cerr << "Create NRLoader " << family_mode << "\n";
}

void NRLoader::start()
{
    thread_pool_->post([this]() {
		    load_families();
		    std::cout << "load complete on " << file_ << "\n";
		    load_state_.pending_dec();
		});
}

void NRLoader::load_families()
{
    std::cerr << "Begin load of " << file_ << "\n";
    cur_work_ = std::make_shared<seq_list_t>();

    size_t fsize = boost::filesystem::file_size(file_);

    max_size_ = fsize / thread_pool_->size() / int(ceil(10.0 / (float) n_files_));
    if (max_size_ < 1000000)
	max_size_ = 1000000;
    cur_size_ = 0;

    {
	boost::lock_guard<boost::mutex> lock(dbg_mut_);
	std::cerr << "tp size=" << thread_pool_->size() << " max_size_=" << max_size_ << "\n";
    }

    FastaParser parser;
    parser.set_callback(std::bind(&NRLoader::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2));

    std::ifstream inp(file_);
    parser.parse(inp);
    parser.parse_complete();

    if (cur_size_)
    {
	auto sent_work = cur_work_;
	cur_work_ = 0;
	cur_size_ = 0;

	my_load_state_.pending_inc();
	std::cerr << "post final work of size " << sent_work->size() << "\n";
	thread_pool_->post(std::bind(&NRLoader::thread_load, this, sent_work, chunks_started_++));
    }

    std::cerr << file_ << " loader awaiting completion of tasks\n";
    my_load_state_.pending_wait();
    std::cerr << file_ << " loader completed\n";
}

int NRLoader::on_parsed_seq(const std::string &id, const std::string &seq)
{
    cur_work_->push_back(std::make_pair(id, seq));
    cur_size_ += seq.size();
    if (cur_size_ >= max_size_)
    {
	auto sent_work = cur_work_;
	cur_work_ = std::make_shared<seq_list_t>();
	cur_size_ = 0;

	my_load_state_.pending_inc();
	std::cerr << file_ << " post work of size " << sent_work->size() << "\n";
	thread_pool_->post(std::bind(&NRLoader::thread_load, this, sent_work, chunks_started_++));
    }

    return 0;
}

void NRLoader::on_hit(KmerGuts::hit_in_sequence_t hit, KmerPegMapping::encoded_id_t &enc_id, size_t seq_len)
{
    if (family_mode_)
    {
	root_mapping_->add_fam_mapping(enc_id, hit.hit.which_kmer);
    }
    else
    {
	root_mapping_->add_mapping(enc_id, hit.hit.which_kmer);
    }
}

void NRLoader::on_hit_fam(KmerGuts::hit_in_sequence_t hit, KmerPegMapping::encoded_family_id_t &enc_id,
			  size_t seq_len)
{
    root_mapping_->add_fam_mapping(enc_id, hit.hit.which_kmer);
}


/*
 * Process a bucket of work in the worker thread.
 *
 * If we are in family load mode, we collect a set of
 * kmer=>family_id mappings to be inserted. In order to save synchronization overhead, we batch
 * these up by encoded-kmer modulo n-insert-workers and at the end of a sequence push the
 * work to the appropriate insert-worker input queue.
 */
void NRLoader::thread_load(std::shared_ptr<seq_list_t> sent_work, int count)
{
    KmerGuts *kguts = thread_pool_->kguts_.get();
    try {
	for (auto seq_entry: *sent_work)
	{
	    std::string &id = seq_entry.first;
	    std::string &seq = seq_entry.second;
	    
	    KmerPegMapping::encoded_id_t enc_id = root_mapping_->encode_id(id);

	    std::function<void(KmerGuts::hit_in_sequence_t)> hit_cb;

	    std::vector<std::shared_ptr<KmerInserter::WorkElement>> insert_work_list;

	    if (family_mode_)
	    {
		for (int i = 0; i < inserter_.n_workers(); i++)
		    insert_work_list.emplace_back(std::make_shared<KmerInserter::WorkElement>());
		
		auto fam_id_iter = root_mapping_->peg_to_family_.find(enc_id);
		if (fam_id_iter == root_mapping_->peg_to_family_.end())
		{
		    std::string dec = root_mapping_->decode_id(enc_id);
		    std::cerr << "NO FAM FOR id='" << id << "' enc_id='" << enc_id << "' dec='" << dec << "'\n";
		    my_load_state_.pending_dec();
		    return;
		}
		KmerPegMapping::encoded_family_id_t fam_id = fam_id_iter->second;
		hit_cb = [this, insert_work_list, fam_id, enc_id](KmerGuts::hit_in_sequence_t hit)
		{
		    int modulus = (int) (hit.hit.which_kmer % (unsigned long long) inserter_.n_workers());
		    insert_work_list[modulus]->work.emplace_back(std::make_pair(hit.hit.which_kmer, fam_id));
		};
		// hit_cb = std::bind(&NRLoader::on_hit_fam, this, std::placeholders::_1, fam_id, seq.length(), insert_work_list);
	    }
	    else
	    {		   
		hit_cb = std::bind(&NRLoader::on_hit, this, std::placeholders::_1, enc_id, seq.length());
	    }

	    kguts->process_aa_seq(id, seq, 0, hit_cb, 0);
	    if (family_mode_)
	    {
		for (int i = 0; i < inserter_.n_workers(); i++)
		{
		    auto w = insert_work_list[i];

		    if (w->work.size() > 0)
		    {
			inserter_.push_work(i, w);
		    }
		}
	    }
		
	}
    }
    catch (std::exception &e)
    {
	std::cerr << "initial load exception " << e.what() << "\n";
    }
    catch (...)
    {
	std::cerr << "initial load default exception\n";
    }
    {
	boost::lock_guard<boost::mutex> lock(dbg_mut_);
	chunks_finished_++;
 	std::cerr << file_ << " finish item " << count << " chunks_finished=" << chunks_finished_ << " chunks_started=" << chunks_started_ << "\n";
    }
    my_load_state_.pending_dec();
}
