#include "family_mapper.h"

FamilyMapper::FamilyMapper(KmerGuts *kguts,
			   std::shared_ptr<KmerPegMapping> mapping) :
    find_best_match_(false),
    family_mode_(true),
    kmer_hit_threshold_(3),
    allow_ambiguous_functions_(false),
    mapping_(mapping),
    target_genus_id_(0),
    find_reps_(false),
    kguts_(kguts)
    
{
}

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

void FamilyMapper::ingest_protein(const std::string &id, const std::string &seq)
{
    seq_score_.clear();

    find_best_match_ = family_mode_ = true;
    /*
     * We accumulate potential function calls if we are going to attempt
     * to return best family match.
     */
    calls_ = 0;
    if (find_best_match_ && family_mode_)
	calls_ = std::make_shared<call_vector_t>();

    // std::cerr << "Lookup " << id << " " << seq << "\n";
    kguts_->process_aa_seq(id, seq, calls_,
			  std::bind(&FamilyMapper::on_hit, this, std::placeholders::_1),
			  0);
}

FamilyMapper::best_match_t FamilyMapper::find_best_family_match(const std::string &id, const std::string &seq)
{
    FunctionIndex best_call_fi;
    float best_call_score, best_weighted_score, best_score_offset;
    std::string best_call_function;

    find_best_match_ = family_mode_ = true;

    ingest_protein(id, seq);

    //
    // find_best_call 
    // It may return a 2-way ambiguous call as "F1 ? F2" 
    //
    /*
     * The KmerGuts find_best_call routine collapses a set of hits to a single best call.
     *
     * It may return a 2-way ambiguous call:
     *
     *   F1 ?? F2
     *
     * if F1 and F2 are too close to call, but are both above the
     * threshold for the next best call.
     *
     * We allow, if allow_ambiguous_functions_ is set, family matching
     * to occur with the best of the two ambiguous functions.
     *
     * If we don't allow ambiguous function matching, we we set
     * the best function to "hypothetical protein" to match the
     * typical behavior. 
     */

    
    kguts_->find_best_call(*calls_, best_call_fi, best_call_function, best_call_score, best_weighted_score, best_score_offset);

    // std::cout << "Best call for " << id << ": " << best_call_fi << " " << best_call_function << " " << best_call_score << "\n";
	
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
	// if (score_ent.weighted_total > best_lf.score && fam_data.genus_id == target_genus_id_)
	if (score_ent.weighted_total > best_lf.score)
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
    // os << id << "\t" << best_gf.fam << "\t" << best_gf.score << "\t" << best_lf.fam << "\t" << best_lf.score << "\t" << (do_ambig_test ? best_lf.function : best_call_function) << "\t" << best_call_score << "\n";

    return { best_gf.fam, best_gf.score, best_lf.fam, best_lf.score, (do_ambig_test ? best_lf.function : best_call_function), best_call_score };
}

void FamilyMapper::find_all_matches(std::ostream &os, const std::string &id, const std::string &seq)
{
    find_best_match_ = false;
    family_mode_ = true;

    ingest_protein(id, seq);
    
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
/*
  if (find_reps_)
  {
  auto reps = owner_->server()->family_reps();
  if (reps)
  {
  auto rep_it = reps->reps_.find(fent.plf);
  if (rep_it != reps->reps_.end())
  {
  for (auto rep: rep_it->second)
  {
  os << rep.feature_id << "\t" << rep.contig << "\t" << rep.contig_length << "\t" << rep.start << "\t"
  << rep.end << "\t" << rep.strand << "\n";
  }
  }
  }
  os << "///\n";
  }
*/
	    }
	    else
	    {
		std::string peg = mapping_->decode_id(eid);
		os << peg << "\t" << score_ent.hit_count;
		
		auto fhit = mapping_->peg_to_family_.find(eid);
		if (fhit != mapping_->peg_to_family_.end())
		{
		    auto fam = mapping_->family_data_[fhit->second];
		    os << "\t" << fam.pgf << "\t" << fam.plf << "\t" << fam.function << "\n";
		    
		}
		else
		{
		    os << "\n";
		}
	    }
	}
    os << "//\n";
}

void FamilyMapper::on_hit(const KmerGuts::hit_in_sequence_t &kmer)
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
	    KmerPegMapping::family_counts_t &counts = ki->second;
	    float weight = 1.0f / (float) counts.size();
	    for (KmerPegMapping::encoded_family_id_t ent : ki->second)
	    {
		//auto fent = mapping_->id_to_family_[ent.first];
		//std::cout << "got ent " << ent.first << " " << fent.second << " with count " << ent.second << "\n";

		// std::cout << "got ent " << ki->first << " " << ent << "\n";
		//sequence_accumulated_score_t &s = seq_score_[ent.first];
		sequence_accumulated_score_t &s = seq_score_[ent];
		s.increment(ent, weight);
	    }
	}
	else
	{
	    // std::cout << "No mapping for k er.hit.which_kmer\n";
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
