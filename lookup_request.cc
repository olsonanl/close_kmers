#include "lookup_request.h"
#include "kserver.h"

#include <string>
#include <boost/bind.hpp>
#include <boost/asio/buffer.hpp>
#include <boost/asio/buffers_iterator.hpp>
#include "global.h"

template<class T>
struct less_second : std::binary_function<T,T,bool>
{
    inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
};  

inline std::string make_string(boost::asio::streambuf& streambuf)
{
    return {buffers_begin(streambuf.data()),
	    buffers_end(streambuf.data())};
}

LookupRequest::LookupRequest(std::shared_ptr<KmerRequest2> owner, std::shared_ptr<KmerPegMapping> mapping,
			     bool family_mode, int content_length, bool chunked) :
    family_mode_(family_mode),
    mapping_(mapping),
    content_length_(content_length),
    chunked_(chunked),
    owner_(owner),
    kmer_hit_threshold_(3),
    header_written_(false),
    find_best_match_(false),
    allow_ambiguous_functions_(false),
    best_match_method_(best_by_accumulation),
    target_genus_(owner->parameters()["target_genus"]),
    target_genus_id_(0)
{
    parser_.set_callback(std::bind(&LookupRequest::on_parsed_seq, this, std::placeholders::_1, std::placeholders::_2));
    try {
	kmer_hit_threshold_ = std::stoi(owner->parameters()["kmer_hit_threhsold"]);
    } catch (const std::invalid_argument& ia) {}
    try {
	find_best_match_ = std::stoi(owner->parameters()["find_best_match"]);
    } catch (const std::invalid_argument& ia) {}
    try {
	allow_ambiguous_functions_ = std::stoi(owner->parameters()["allow_ambiguous_functions"]);
    } catch (const std::invalid_argument& ia) {}

    const std::string &m = owner->parameters()["best_match_method"];
    if (m == "count")
	best_match_method_ = best_by_count;
    else if (m == "weighted_count")
	best_match_method_ = best_by_weighted_count;
    else if (m == "accumulation")
	best_match_method_ = best_by_accumulation;
    else if (!m.empty())
	std::cerr << "Warning: best_match_method '" << m << "' not recognized\n";
    
    try {
	std::string tg(mapping_->lookup_genus(target_genus_));
	if (!tg.empty())
	{
	    target_genus_id_ = std::stoul(tg);
	}
    } catch (const std::invalid_argument& ia) {}

    // std::cerr << "created LookupRequest kmer_hit_threhsold=" << kmer_hit_threshold_ << " find_best_match=" << find_best_match_ << " target_genus=" << target_genus_ << " target_genus_id_=" << target_genus_id_ << " allow_ambiguous_functions_=" << allow_ambiguous_functions_ << "\n";
}

void LookupRequest::run()
{
    if (owner_->request().size() > 0)
    {
	boost::system::error_code err;
	on_data(err, 0);
    }
    else
    {
	boost::asio::async_read(owner_->socket(), owner_->request(),
				boost::asio::transfer_at_least(content_length_),
				boost::bind(&LookupRequest::on_data, this,
					    boost::asio::placeholders::error,
					    boost::asio::placeholders::bytes_transferred));
    }
}

void LookupRequest::on_data(boost::system::error_code err, size_t bytes)
{
    if (!err || err == boost::asio::error::eof)
    {
 	// std::string r = make_string(owner_->request());
	// std::cerr << "bytes=" << bytes << " content_length_=" << content_length_ << " err=" << err <<  "\n";
	// std::cerr << "Read buffer contains: "  << r << std::endl;
	
	current_work_ = std::make_shared<work_list_t>();

	boost::asio::streambuf::const_buffers_type bufs = owner_->request().data();
	int n = 0;
	for (auto x = boost::asio::buffers_begin(bufs); x != boost::asio::buffers_end(bufs); x++)
	{
	    parser_.parse_char(*x);
	    n++;
	}

	content_length_ -= n;

	// std::cerr << "Processed " << n << " remaining " << content_length_ << "\n";
	if (err == boost::asio::error::eof || content_length_ == 0)
	{
	    parser_.parse_complete();
	}

	owner_->request().consume(n);

	/*
	 * We have parsed out this packet.
	 * We post the data to the thread pool for processing, and the
	 * results are to be posted back here so we can push them to the
	 * outgoing socket and queue the next read.
	 */

	std::shared_ptr<work_list_t> cur = current_work_;
	current_work_ = 0;
	owner_->thread_pool()->post([this, cur, err]() {
		KmerGuts *kguts = owner_->thread_pool()->kguts_.get();
		kguts->set_parameters(owner_->parameters());

		auto sbuf = std::make_shared<boost::asio::streambuf>();
		{
		    std::ostream os(sbuf.get());

		    if (!header_written_)
		    {
			owner_->write_header(os, 200, "OK");
			os << "\n";
			header_written_ = true;
		    }

		    for (auto work_item: *cur)
		    {
			auto id = work_item.id();
			auto seq = work_item.seq();

			seq_score_.clear();

			/*
			 * We accumulate potential function calls if we are going to attempt
			 * to return best family match.
			 */
			typedef std::vector<KmerCall> call_vector_t;
			std::shared_ptr<call_vector_t> calls = 0;
			if (find_best_match_ && family_mode_)
			    calls = std::make_shared<call_vector_t>();

			// std::cerr << "Lookup " << id << " " << seq << "\n";
			kguts->process_aa_seq(id, seq, calls,
					      std::bind(&LookupRequest::on_hit, this, std::placeholders::_1),
					      0);

			/*
			 * Process results for this sequence.
			 *
			 * Our putative family matches have been collected into hit_count and hit_total;
			 * these are maps from either encoded family id (if we are in family mode)
			 * or protein id (if we are not in family mode).
			 *
			 * If we are in find_best_match mode, we determine the called function
			 * for this protein using KmerGuts::find_best_call(). We hereby declare that
			 * find_best_match is only respected in family mode.
			 *
			 * We need to find the best match for both global and local families. For the
			 * local family match we must have been provided a genus name in the query
			 * via the "target_genus" parameter.
			 *
			 * Since the family results we get are pgf/plf pairs, in order to vote on the
			 * best family we must accumulate counts for the pgfs. We can select the best
			 * plf by choosing the plf with the highest score that matches the target genus.
			 *
			 * This requires one scan over the matched families, plus overhead for maintaining
			 * the pgf->score mapping plus a scan over pgfs to determine the best value.
			 *
			 * If we are not in find_best_match mode, we will return all potential family
			 * matches, sorted by score. We do this by converting the hit count map into
			 * a vector and sorting the vector via the less_second operator (element.second
			 * on the map entries is the count value).
			 *
			 */

			if (find_best_match_ && family_mode_)
			{
			    int best_call_fi;
			    float best_call_score;
			    std::string best_call_function;
			    kguts->find_best_call(*calls, best_call_fi, best_call_function, best_call_score);
			    // std::cout << "Best call for " << id << ": " << best_call_fi << " " << best_call_function << " " << best_call_score << "\n";

			    /*
			     * The KmerGuts best-function caller may return a function of the form
			     *   function1 ?? function2
			     * This means there was an ambiguous call between the two functions.
			     *
			     * We allow, if allow_ambiguous_functions_ is set, family matching
			     * to occur with the best of the two ambiguous functions.
			     *
			     * If we don't allow ambiguous function matching, we we set
			     * the best function to "hypothetical protein" to match the
			     * typical behavior. 
			     */

			    std::string ambig_function;
			    bool do_ambig_test = false;
			    if (best_call_function.empty())
				best_call_function = "hypothetical protein";
			    else
			    {
				size_t where = best_call_function.find(" ?? ");
				if (where != std::string::npos)
				{
				    if (allow_ambiguous_functions_)
				    {
					ambig_function = best_call_function.substr(where + 4);
					best_call_function = best_call_function.substr(0, where);
					do_ambig_test = true;
					// std::cerr << "allow ambiguous lookup '" << best_call_function << "' '" << ambig_function << "'\n";
				    }
				    else
				    {
					best_call_function = "hypothetical protein";
				    }
				}
			    }
			    
			    struct top_score
			    {
				float score;
				std::string fam;
				std::string function;
			    };

			    top_score best_lf({score: 0.0 });
			    top_score best_gf({score: 0.0 });

			    std::unordered_map<std::string, float> pgf_rollup, pgf_rollup_ambig;

				
			    for (auto hit_ent: seq_score_)
			    {
				KmerPegMapping::encoded_id_t eid = hit_ent.first;
				const sequence_accumulated_score_t &score_ent = hit_ent.second;

				if (score_ent.hit_total < kmer_hit_threshold_)
				    continue;
				
				KmerPegMapping::family_id_to_family_map_t::iterator fent = mapping_->family_data_.find(eid);
				if (fent == mapping_->family_data_.end())
				    continue;

				const KmerPegMapping::family_data_t &fam_data = fent->second;
				if (do_ambig_test)
				{
				    if (fam_data.function == best_call_function)
				    {
					pgf_rollup[fam_data.pgf] += score_ent.weighted_total;
				    }
				    else if (fam_data.function == ambig_function)
				    {
					pgf_rollup_ambig[fam_data.pgf] += score_ent.weighted_total;
				    }
				    else
				    {
					continue;
				    }
				}
				else
				{
				    if (fam_data.function == best_call_function)
					pgf_rollup[fam_data.pgf] += score_ent.weighted_total;
				    else
					continue;
				}


				// std::cerr << score << " " << best_lf.score << " " << fam_data.genus_id << " " << target_genus_id_ << "\n";
				if (score_ent.weighted_total > best_lf.score && fam_data.genus_id == target_genus_id_)
				{
				    best_lf.score = score_ent.weighted_total;
				    best_lf.fam = fam_data.plf;
				    best_lf.function = fam_data.function;
				}				    
			    }

			    std::unordered_map<std::string, float> *matching_rollup = &pgf_rollup;
			    if (do_ambig_test && best_lf.function == ambig_function)
				matching_rollup = &pgf_rollup_ambig;
			    for (auto pgf_ent: *matching_rollup)
			    {
				const std::string &pgf = pgf_ent.first;
				const float &score = pgf_ent.second;

				if (score > best_gf.score)
				{
				    best_gf.score = score;
				    best_gf.fam = pgf;
				}
			    }

			    /*
			     * We have found our best scores. Write output.
			     */
			    os << id << "\t" << best_gf.fam << "\t" << best_gf.score << "\t" << best_lf.fam << "\t" << best_lf.score << "\t" << (do_ambig_test ? best_lf.function : best_call_function) << "\t" << best_call_score << "\n";
			}
			else
			{
			    typedef std::pair<KmerPegMapping::encoded_id_t, sequence_accumulated_score_t> data_t;
			    
			    std::vector<data_t> vec;
			    for (auto it: seq_score_)
			    {
				vec.push_back(it);
			    }
			    
			    std::sort(vec.begin(), vec.end(), [](const data_t &lhs, const data_t &rhs) {
				    return lhs.second.weighted_total > rhs.second.weighted_total;
				});
			    
			    os << id << "\n";
			    for (auto it: vec)
			    {
				auto eid = it.first;
				const sequence_accumulated_score_t &score_ent = it.second;
				
				if (score_ent.hit_total < kmer_hit_threshold_)
				    break;
				
				if (family_mode_)
				{
				    /*
				     * To report we map the id back to its family and report data from there.
				     */
				    unsigned int score = score_ent.hit_count;
				    unsigned int total = score_ent.hit_total;
				    float weighted = score_ent.weighted_total;
				    auto fent = mapping_->family_data_[eid];
				    float scaled = (float) score / (float) fent.total_size;
				    os << score << "\t" << total << "\t" << weighted << "\t" << fent.pgf << "\t" << fent.plf << "\t" << fent.total_size << "\t" << fent.count << "\t" << scaled << "\t" << fent.function << "\n";
				}
				else
				{
				    std::string peg = mapping_->decode_id(eid);
				    os << peg << "\t" << score_ent.hit_count;
				    // os << eid << "\t" << peg << "\t" << score;
				    
				    auto fhit = mapping_->peg_to_family_.find(eid);
				    if (fhit != mapping_->peg_to_family_.end())
				    {
					auto fam = mapping_->family_data_[fhit->second];
					os << "\t" << fam.pgf << "\t" << fam.plf << "\t" << fam.function;
				    }
				    os << "\n";
				}
			    }
			    os << "//\n";
			}

		    }
		    os.flush();
		}
		owner_->io_service().post([this, sbuf, err](){
			/*
			 * Back in the main thread here. We can write our response.
			 */
			// std::cout << "post response in " << pthread_self() << "\n";

			// std::cerr << "write results size " << sbuf->size() << "\n";
			boost::asio::async_write(owner_->socket(), boost::asio::buffer(sbuf->data()),
						 [this, err, sbuf](const boost::system::error_code &err2, const long unsigned int &bytes2){
						     // std::cerr << "write done in " << pthread_self() << " err=" << err.message() << " content_length=" << content_length_ <<"\n";
						     // std::cerr << "   err2=" << err2.message() << " bytes2=" << bytes2 << "\n";
						     if (err == boost::asio::error::eof || content_length_ == 0)
						     {
							 process_results();
						     }
						     else
						     {
							 boost::asio::async_read(owner_->socket(), owner_->request(),
										 boost::asio::transfer_at_least(content_length_),
										 boost::bind(&LookupRequest::on_data, this,
											     boost::asio::placeholders::error,
											     boost::asio::placeholders::bytes_transferred));
						     }
						 });
		    });
	    });
    }
    else
    {
	std::cerr << "ERROR in lookup_request: " << err << "\n";
    }
	
}

int LookupRequest::on_parsed_seq(const std::string &id, const std::string &seq)
{
    current_work_->push_back(ProteinSequence(id, seq));
    return 0;
}

void LookupRequest::on_hit(KmerGuts::hit_in_sequence_t kmer)
{
    if (family_mode_)
    {
	/*
	 * In family mode we look at the mapping's kmer_to_family_id map
	 * to determine the families that map to that kmer. We roll up the
	 * hit counts for each of those families using the value in the map.
	 */
	KmerPegMapping::family_map_type_t::iterator ki = mapping_->kmer_to_family_id_.find(kmer.hit.which_kmer);
	if (ki != mapping_->kmer_to_family_id_.end())
	{
	    float weight = 1.0f / (float) ki->second.size();
	    for (auto ent : ki->second)
	    {
		// auto fent = mapping_->id_to_family_[ent.first];
		// std::cout << "got ent " << ent.first << " " << fent.second << " with count " << ent.second << "\n";
		sequence_accumulated_score_t &s = seq_score_[ent.first];
		s.hit_count += ent.second;
		s.hit_total++;
		s.weighted_total += weight;
	    }
	}
    }
    else
    {
	auto ki = mapping_->kmer_to_id_.find(kmer.hit.which_kmer);
	if (ki != mapping_->kmer_to_id_.end())
	{
	    // std::cout << "got mapping for " << kmer.hit.which_kmer << "\n";
	    for (auto eid: ki->second)
	    {
		seq_score_[eid].hit_count++;
	    }
	}
    }
}

void LookupRequest::process_results()
{
    owner_->socket().close();
    owner_->exit_request();
}
