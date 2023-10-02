#ifndef GEESE_MEAT_LIKELIHOOD_HPP
#define GEESE_MEAT_LIKELIHOOD_HPP 1

#include "geese-bones.hpp"

inline double Geese::likelihood(
    const std::vector< double > & par,
    bool as_log,
    bool use_reduced_sequence,
    size_t ncores,
    bool no_update_normalizing_constant
) {

    // Checking whether the model is initialized
    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    double ll = 0.0;

    // Updating normalizing constants
    if (!no_update_normalizing_constant)
        model->update_normalizing_constants(par0);

    // Following the prunning sequence
    std::vector< size_t > * preseq;

    if (use_reduced_sequence)
    {

        preseq = &this->reduced_sequence;

    }
    else
    {   

        preseq = &this->sequence;

    }

    // The first time it is called, it need to generate the corresponding
    // hashes of the columns so it is fast to access then (saves time
    // hashing and looking in the map.)
    auto arrays2support = model->get_arrays2support();

    for (auto& i : *preseq)
    {

        // We cannot compute probability at the leaf, we need to continue
        if (this->nodes[i].is_leaf())
            continue;

        // Since we are using this a lot...
        Node & node = nodes[i];

        // Iterating through states
        for (size_t s = 0u; s < states.size(); ++s)
        {

            // Starting the prob
            size_t array_id = node.narray[s];

            // Retrieving the sets of arrays
            const std::vector< PhyloArray > & psets =
                *(model->get_pset(array_id));

            const std::vector<double> & psets_stats =
                *(model->get_pset_stats(array_id));

            std::vector< std::vector< size_t > > & locations = pset_loc[
                arrays2support->operator[](array_id)
                ];
            
            // Summation over all possible values of X
            const auto & node_offspring = node.offspring;
            std::vector< double > totprob_n(psets.size(), 0.0);
            size_t nstate = 0u;
            #if defined(_OPENMP) || defined(__OPENMP)
            #pragma omp parallel for num_threads(ncores) \
                shared(locations, psets, psets_stats, totprob_n) \
                firstprivate(nfunctions, par0, node_offspring, array_id)
            #endif
            for (size_t n = 0u; n < psets.size(); ++n) // x = psets->begin(); x != psets->end(); ++x)
            {

                // Retrieving the pset
                const auto & x = psets[n];

                if (!x.is_dense())
                    throw std::logic_error("This is only supported for dense arrays.");

                const std::vector< size_t > & location_x = locations[n];

                // Extracting the possible values of each offspring
                double off_mult = 1.0;

                for (auto o = 0u; o < x.ncol(); ++o)
                {

                    // Setting the node
                    const Node * n_off = node_offspring[o];
                    
                    // In the case that the offspring is a leaf, then we need to
                    // check whether the state makes sense.
                    if (n_off->is_leaf())
                    {
                        for (auto f = 0u; f < nfunctions; ++f)
                        {
                            if (n_off->annotations[f] != 9u)
                            {

                                if (x(f, o) != n_off->annotations[f])
                                {

                                    off_mult = -1.0;
                                    break;

                                }
                                
                            }

                        }

                        // Going out
                        if (off_mult < 0)
                            break;
                
                        continue;

                    }

                    // Retrieving the location to the respective set of probabilities
                    off_mult *= n_off->subtree_prob[location_x[o]];

                }

                // Is this state valid?
                if (off_mult < 0.0)
                {

                    ++nstate;
                    continue;
                    
                }

                // Multiplying by P(x|x_n), the transition probability
                std::vector< double > temp_stats;
                temp_stats.reserve(par0.size());
                for (auto p = 0u; p < par0.size(); ++p)
                    temp_stats.push_back(psets_stats[par0.size() * nstate + p]);

                nstate++;

                // Use try catch in the following line
                try {
                    off_mult *= model->likelihood(
                        par0,
                        temp_stats,
                        array_id,
                        false,
                        true
                    );
                } catch (std::exception & e) {

                    auto err = std::string(e.what());

                    std::string state_str = "";
                    for (const auto & ss : states[s])
                        state_str += std::to_string(ss) + " ";

                    err = "Error computing the likelihood at node " +
                        std::to_string(node.id) + " with state " + state_str +
                        ". Error message:\n" +
                        err;

                    throw std::runtime_error(err);
                    
                }

                // Adding to the total probabilities
                totprob_n[n] = off_mult;

            }

            // Setting the probability at the node
            node.subtree_prob[s] = 0.0;
            for (size_t n = 0u; n < psets.size(); ++n)
                node.subtree_prob[s] += totprob_n[n];

        }

        // All probabilities should be completed at this point
        if (node.parent == nullptr)
        {

            for (size_t s = 0u; s < states.size(); ++s)
            {

                double tmpll = 1.0;

                for (auto k = 0u; k < nfunctions; ++k)
                {

                    tmpll *= states[s][k] ? par_root[k] : (1 - par_root[k]);

                }

                ll += tmpll * node.subtree_prob[s];

            }
        }

    }

    // In the case that the sequence is empty, then it means
    // that we are looking at a completely unnanotated tree,
    // thus the likelihood should be one
    if (preseq->size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 1.0;


    return as_log ? std::log(ll) : ll;

}
#endif