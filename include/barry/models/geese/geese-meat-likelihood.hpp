#ifndef GEESE_MEAT_LIKELIHOOD_HPP
#define GEESE_MEAT_LIKELIHOOD_HPP 1

#include "geese-bones.hpp"

inline void pset_loop(
    size_t n,
    size_t s,
    size_t nfunctions,
    const size_t node_id,
    double norm_const_i,
    std::vector< double > & totprob_n,
    const std::vector< double > & par0,
    const std::vector<std::vector<bool>> & states,
    const std::vector< PhyloArray > & psets,
    const std::vector<double> & psets_stats,
    const std::vector< std::vector< size_t > > & locations,
    const std::vector<geese::Node *> & node_offspring
) 
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
        return;

    // Use try catch in the following line
    try {

        off_mult *= barry::likelihood_(
            &psets_stats[par0.size() * n],
            par0,
            norm_const_i,
            par0.size(),
            false
        );

    } catch (std::exception & e) {

        auto err = std::string(e.what());

        std::string state_str = "";
        for (const auto & ss : states[s])
            state_str += std::to_string(ss) + " ";

        err = "Error computing the likelihood at node " +
            std::to_string(node_id) + " with state " + state_str +
            ". Error message:\n" +
            err;

        throw std::runtime_error(err);
        
    }

    // Adding to the total probabilities
    totprob_n[n] = off_mult;

}

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
        model->update_normalizing_constants(par0, ncores);

    // Following the prunning sequence
    const std::vector< size_t > & preseq = use_reduced_sequence ?
        this->reduced_sequence : this->sequence;

    // The first time it is called, it need to generate the corresponding
    // hashes of the columns so it is fast to access then (saves time
    // hashing and looking in the map.)
    const auto & arrays2support = *(model->get_arrays2support());
    const auto & normconst = model->get_normalizing_constants();

    for (auto& i : preseq)
    {

        // We cannot compute probability at the leaf, we need to continue
        if (this->nodes[i].is_leaf())
            continue;

        // Since we are using this a lot...
        Node & node = nodes[i];
        const size_t node_id  = node.id;

        // Iterating through states
        for (size_t s = 0u; s < states.size(); ++s)
        {

            // Starting the prob
            size_t array_id = node.narray[s];
            size_t support_id = arrays2support[array_id];
            double norm_const_i = normconst[support_id];

            // Retrieving the sets of arrays
            const std::vector< PhyloArray > & psets =
                *(model->get_pset(array_id));

            const std::vector<double> & psets_stats =
                *(model->get_pset_stats(array_id));

            std::vector< std::vector< size_t > > & locations = pset_loc[support_id];

            // Making sure parallelization makes sense
            if (psets.size() < 128)
                ncores = 1u;
            
            // Summation over all possible values of X
            const auto & node_offspring = node.offspring;
            std::vector< double > totprob_n(psets.size(), 0.0);
            #if defined(_OPENMP) || defined(__OPENMP)
            if (ncores > 1u)
            {
                #pragma omp parallel for num_threads(ncores) \
                    shared(\
                        locations, psets, psets_stats, totprob_n, node, states,\
                        par0, node_offspring, nfunctions, array_id, norm_const_i, \
                        s, node_id) default(none)
                for (size_t n = 0u; n < psets.size(); ++n) 
                {
                    pset_loop(
                        n, s, nfunctions, node_id, norm_const_i, totprob_n,
                        par0, states, psets, psets_stats, locations, 
                        node_offspring
                    );
                }
            } else {
                for (size_t n = 0u; n < psets.size(); ++n) 
                {
                    pset_loop(
                        n, s, nfunctions, node_id, norm_const_i, totprob_n,
                        par0, states, psets, psets_stats, locations, 
                        node_offspring
                    );
                }
            }
            #else
            for (size_t n = 0u; n < psets.size(); ++n) 
            {
                pset_loop(
                    n, s, nfunctions, node_id, norm_const_i, totprob_n,
                    par0, states, psets, psets_stats, locations, 
                    node_offspring
                );
            }
            #endif
            

            // Setting the probability at the node
            node.subtree_prob[s] = 0.0;
            auto & nsp = node.subtree_prob[s];
            #if defined(_OPENMP) || defined(__OPENMP)
            #pragma omp simd reduction(+:nsp)
            #endif
            for (size_t n = 0u; n < psets.size(); ++n)
                nsp += totprob_n[n];


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
    if (preseq.size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 1.0;


    return as_log ? std::log(ll) : ll;

}
#endif