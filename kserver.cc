
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
    thread_pool_(thread_pool)
{
    mapping_map_ = std::make_shared<std::map<std::string, std::shared_ptr<KmerPegMapping> > >();

    auto x = mapping_map_->insert(std::make_pair("", std::make_shared<KmerPegMapping>()));
    std::cerr << "insert created new: " << x.second << " key='" << x.first->first << "'\n";
    auto root_mapping = x.first->second;
    
    root_mapping->load_genome_map((*g_parameters)["kmer-data-dir"].as<std::string>() + "/genomes");

    ld_pending_count_= 0;

    /*
     * If we are preloading a families file, start that off in the background
     * using the thread pool.
     */
    if (g_parameters->count("families-file"))
    {
	std::string ff = (*g_parameters)["families-file"].as<std::string>();
	{
	    boost::lock_guard<boost::mutex> lock(ld_mut_);
	    ld_pending_count_++;
	}
	thread_pool_->post([this,root_mapping, ff]() {
		std::cerr << "Loading families from " << ff << "\n";
		root_mapping->load_families(ff);
		{
		    boost::lock_guard<boost::mutex> lock(ld_mut_);
		    std::cerr << "finishing family file load pending=" << ld_pending_count_ << "\n";
		    ld_pending_count_--;
		}
		ld_cond_.notify_one();
	    });
    }

    if (g_parameters->count("reserve-mapping"))
    {
	int count = (*g_parameters)["reserve-mapping"].as<int>();
	std::cerr << "Reserving " << count << " bytes in mapping table\n";
	root_mapping->reserve_mapping_space(count);
    }

    if (g_parameters->count("families-nr"))
    {
	load_families_nr(root_mapping, (*g_parameters)["families-nr"].as<std::string>());
    }
	
    std::cerr << "wait for threads to finish\n";

    boost::unique_lock<boost::mutex> lock(ld_mut_);
    while (ld_pending_count_ > 0)
    {
	std::cerr << "pending=" << ld_pending_count_ << "\n";
	ld_cond_.wait(lock);
    }
    
    std::cerr << "done waiting\n";
    ld_cur_work_ = 0;
}

/*
 * Preload a families NR file.
 *
 * We do this by parsing out the file into chunks and feeding those chunks to
 * the threadpool to parse and insert.
 */
void KmerRequestServer::load_families_nr(std::shared_ptr<KmerPegMapping> &mapping,
					 const std::string &file)
{
    ld_root_mapping_ = mapping;

    ld_cur_work_ = std::make_shared<seq_list_t>();

    size_t fsize = boost::filesystem::file_size(file);

    ld_max_size_ = fsize / thread_pool_->size() / 10;
    if (ld_max_size_ < 1000000)
	ld_max_size_ = 1000000;
    ld_cur_size_ = 0;

    FastaParser parser(std::bind(&KmerRequestServer::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2));

    // FastaParser parser([this,&cur_work,&mapping,&cur_size,max_size,&mut,&cond,&pending_count](const std::string &id, const std::string &seq)

    std::ifstream inp(file);
    parser.parse(inp);
    parser.parse_complete();
}

int KmerRequestServer::on_parsed_seq(const std::string &id, const std::string &seq)
{
    ld_cur_work_->push_back(std::make_pair(id, seq));
    ld_cur_size_ += seq.size();
    if (ld_cur_size_ >= ld_max_size_)
    {
	auto sent_work = ld_cur_work_;
	ld_cur_work_ = std::make_shared<seq_list_t>();
	ld_cur_size_ = 0;
	
	{
	    boost::lock_guard<boost::mutex> lock(ld_mut_);
	    ld_pending_count_++;
	}
	
	std::cerr << "post work of size " << sent_work->size() << "\n";
	thread_pool_->post(std::bind(&KmerRequestServer::thread_load, this, sent_work));
    }

}

void KmerRequestServer::on_hit(sig_kmer_t &hit, KmerPegMapping::encoded_id_t &enc_id)
{
    ld_root_mapping_->add_mapping(enc_id, hit.which_kmer);
}

void KmerRequestServer::thread_load(std::shared_ptr<seq_list_t> sent_work)
{
    KmerGuts *kguts = thread_pool_->kguts_.get();
    try {
	for (auto seq_entry: *sent_work)
	{
	    std::string &id = seq_entry.first;
	    std::string &seq = seq_entry.second;
	    
	    KmerPegMapping::encoded_id_t enc_id = ld_root_mapping_->encode_id(id);
	    
	    auto hit_cb = std::bind(&KmerRequestServer::on_hit, this, std::placeholders::_1, enc_id);

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
	boost::lock_guard<boost::mutex> lock(ld_mut_);
	std::cerr << "finishing pending=" << ld_pending_count_ << "\n";
	ld_pending_count_--;
    }
    ld_cond_.notify_one();
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
