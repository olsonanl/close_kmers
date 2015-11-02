
#include "kserver.h"

#include <boost/asio.hpp>
#include <boost/filesystem.hpp>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include "global.h"

KmerRequestServer::KmerRequestServer(boost::asio::io_service& io_service,
				     const std::string &port,
				     const std::string &port_file,
				     std::shared_ptr<ThreadPool> thread_pool
    ) :
    io_service_(io_service),
    acceptor_(io_service_),
    port_(port),
    port_file_(port_file),
    signals_(io_service_),
    thread_pool_(thread_pool),
    load_state_("master")
{
    
    mapping_map_ = std::make_shared<std::map<std::string, std::shared_ptr<KmerPegMapping> > >();

    auto x = mapping_map_->insert(std::make_pair("", std::make_shared<KmerPegMapping>()));
    // std::cerr << "insert created new: " << x.second << " key='" << x.first->first << "'\n";
    auto root_mapping = x.first->second;
    
    root_mapping->load_genome_map((*g_parameters)["kmer-data-dir"].as<std::string>() + "/genomes");

    /*
     * If we are preloading a families file, start that off in the background
     * using the thread pool.
     */
    if (g_parameters->count("families-file"))
    {
	std::string ff = (*g_parameters)["families-file"].as<std::string>();
	load_state_.pending_inc();
	thread_pool_->post([this,root_mapping, ff]() {
		std::cerr << "Loading families from " << ff << "...\n";
		root_mapping->load_families(ff);
		std::cerr << "Loading families from " << ff << "... done\n";
		load_state_.pending_dec();
	    });
    }

    if (g_parameters->count("reserve-mapping"))
    {
	int count = (*g_parameters)["reserve-mapping"].as<int>();
	std::cerr << "Reserving " << count << " bytes in mapping table\n";
	root_mapping->reserve_mapping_space(count);
    }

    std::vector<NRLoader *> active_loaders;

    if (g_parameters->count("families-nr"))
    {
	auto files = (*g_parameters)["families-nr"].as<std::vector<std::string> >();
	for (auto file: files)
	{
	    load_state_.pending_inc();
	    NRLoader *loader = new NRLoader(load_state_, file, root_mapping, thread_pool_, files.size());
	    active_loaders.push_back(loader);
	    loader->start();
	}
    }
	
    std::cerr << "wait for threads to finish\n";

    load_state_.pending_wait();
    
    std::cerr << "done waiting\n";

    for (auto v: active_loaders)
    {
	// std::cerr << "remove loader for " << v->file_ << "\n";
	delete v;
    }
}

/*
 * Need to split startup code from constructor due to use of enable_shared_from_this().
 */
void KmerRequestServer::startup()
{

    /*
     * Set up for clean signal handling / termination
     */
    signals_.add(SIGINT);
    signals_.add(SIGTERM);
    signals_.add(SIGQUIT);
    do_await_stop();

    /*
     * Set up listener
     */

    boost::asio::ip::tcp::resolver resolver(io_service_);
    boost::asio::ip::tcp::endpoint endpoint = *resolver.resolve({"0.0.0.0", port_});
    acceptor_.open(endpoint.protocol());
    acceptor_.set_option(boost::asio::ip::tcp::acceptor::reuse_address(true));
    acceptor_.bind(endpoint);
    acceptor_.listen();
    std::cout << "Listening on " << acceptor_.local_endpoint() << "\n";
    if (!port_file_.empty())
    {
	std::ofstream out(port_file_);
	out << acceptor_.local_endpoint().port() << "\n";
	out.close();
    }
	    
    do_accept2();
}

void KmerRequestServer::do_accept2()
{
    std::shared_ptr<KmerRequest2> r = std::make_shared<KmerRequest2>(shared_from_this(), io_service_, mapping_map_, thread_pool_);
    // std::cerr << "create " << r << "\n";
    acceptor_.async_accept(r->socket(),
			   boost::bind(&KmerRequestServer::on_accept2, this,
				       boost::asio::placeholders::error, r));
    //std::cerr << "leaving do_accept2 r use count=" << r.use_count() << "\n";
}

void KmerRequestServer::on_accept2(boost::system::error_code ec, std::shared_ptr<KmerRequest2> r)
{
    // std::cerr << "on-accept2 use=" << r.use_count() << "\n";
    g_timer.start();
    // Check whether the server was stopped by a signal before this
    // completion handler had a chance to run.
    if (!acceptor_.is_open())
    {
	std::cout << "not open\n";
	return;
    }
    
    if (!ec)
    {
	/*
	 * Connection has come in.
	 * Begin parsing the request line and headers.
	 */
	
	active_.insert(r);
	// std::cerr << "start read " << r << "\n";
	r->do_read();
    }
    
    do_accept2();
}

void KmerRequestServer::deactivate(std::shared_ptr<KmerRequest2> x)
{
    active_.erase(x);
}

void KmerRequestServer::do_await_stop()
{
    signals_.async_wait(
	[this](boost::system::error_code ec, int signo)
	{
	    std::cout << "Exiting with signal " << signo << "\n";
	    acceptor_.close();
	});
}

/*
 * Preload a families NR file.
 *
 * We do this by parsing out the file into chunks and feeding those chunks to
 * the threadpool to parse and insert.
 */

NRLoader::NRLoader(NRLoadState &load_state, const std::string &file,
		   std::shared_ptr<KmerPegMapping> root_mapping,
		   std::shared_ptr<ThreadPool> thread_pool,
		   int n_files) :
    load_state_(load_state),
    file_(file),
    root_mapping_(root_mapping),
    thread_pool_(thread_pool),
    chunks_started_(0),
    chunks_finished_(0),
    my_load_state_(file),
    n_files_(n_files)
{
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
    // std::cerr << "Begin load of " << file_ << "\n";
    cur_work_ = std::make_shared<KmerRequestServer::seq_list_t>();

    size_t fsize = boost::filesystem::file_size(file_);

    max_size_ = fsize / thread_pool_->size() / int(ceil(10.0 / (float) n_files_));
    if (max_size_ < 1000000)
	max_size_ = 1000000;
    cur_size_ = 0;

    FastaParser parser(std::bind(&NRLoader::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2));

    std::ifstream inp(file_);
    parser.parse(inp);
    parser.parse_complete();

    if (cur_size_)
    {
	auto sent_work = cur_work_;
	cur_work_ = 0;
	cur_size_ = 0;

	my_load_state_.pending_inc();
	//std::cerr << "post final work of size " << sent_work->size() << "\n";
	thread_pool_->post(std::bind(&NRLoader::thread_load, this, sent_work, chunks_started_++));
    }

    // std::cerr << file_ << " loader awaiting completion of tasks\n";
    my_load_state_.pending_wait();
    // std::cerr << file_ << " loader completed\n";
}

int NRLoader::on_parsed_seq(const std::string &id, const std::string &seq)
{
    cur_work_->push_back(std::make_pair(id, seq));
    cur_size_ += seq.size();
    if (cur_size_ >= max_size_)
    {
	auto sent_work = cur_work_;
	cur_work_ = std::make_shared<KmerRequestServer::seq_list_t>();
	cur_size_ = 0;

	my_load_state_.pending_inc();
	// std::cerr << file_ << " post work of size " << sent_work->size() << "\n";
	thread_pool_->post(std::bind(&NRLoader::thread_load, this, sent_work, chunks_started_++));
    }

}

void NRLoader::on_hit(KmerGuts::hit_in_sequence_t hit, KmerPegMapping::encoded_id_t &enc_id)
{
    root_mapping_->add_mapping(enc_id, hit.hit.which_kmer);
}

void NRLoader::thread_load(std::shared_ptr<KmerRequestServer::seq_list_t> sent_work, int count)
{
    KmerGuts *kguts = thread_pool_->kguts_.get();
    try {
	for (auto seq_entry: *sent_work)
	{
	    std::string &id = seq_entry.first;
	    std::string &seq = seq_entry.second;
	    
	    KmerPegMapping::encoded_id_t enc_id = root_mapping_->encode_id(id);
	    
	    auto hit_cb = std::bind(&NRLoader::on_hit, this, std::placeholders::_1, enc_id);

	    kguts->process_aa_seq(id, seq, 0, hit_cb, 0);
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
// 	std::cerr << file_ << " finish item " << count << " chunks_finished=" << chunks_finished_ << " chunks_started=" << chunks_started_ << "\n";
    }
    my_load_state_.pending_dec();
}

