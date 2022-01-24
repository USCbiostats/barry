#ifndef GEESE_MEAT_LIKELIHOOD_HPP
#define GEESE_MEAT_LIKELIHOOD_HPP 1

#include "geese-bones.hpp"

inline double Geese::likelihood(
    const std::vector< double > & par,
    bool as_log,
    bool use_reduced_sequence
) {

    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    std::vector< unsigned int > tmpstate(nfunctions);

    double ll = 0.0;

    Node * n_off;

    // Following the prunning sequence
    std::vector< unsigned int > * preseq;

    if (use_reduced_sequence)
    {

        preseq = &this->reduced_sequence;

    }
    else
    {   

        preseq = &this->sequence;

    }


    for (auto& i : *preseq)
    {

        // We cannot compute probability at the leaf, we need to continue
        if (this->nodes[i].is_leaf())
            continue;

        // Since we are using this a lot...
        Node & node = nodes[i];

        // Iterating through states
        for (unsigned int s = 0u; s < states.size(); ++s)
        {

            // Starting the prob
            double totprob = 0.0;

            // Retrieving the sets of arrays
            const std::vector< phylocounters::PhyloArray > * psets =
                model->get_pset(node.narray[s]);

            const std::vector< std::vector<double> > * psets_stats =
                model->get_pset_stats(node.narray[s]);

            // Summation over all possible values of X
            unsigned int nstate = 0u;

            for (auto x = psets->begin(); x != psets->end(); ++x)
            {

                // Extracting the possible values of each offspring
                double off_mult = 1.0;

                for (auto o = 0u; o < x->ncol(); ++o)
                {

                    // Setting the node
                    n_off = node.offspring[o];

                    // First, getting what is the corresponding state
                    std::fill(tmpstate.begin(), tmpstate.end(), 0u);

                    x->get_col_vec(&tmpstate, o, false);

                    // In the case that the offspring is a leaf, then we need to
                    // check whether the state makes sense.
                    if (n_off->is_leaf())
                    {

                        if (vec_diff(tmpstate, n_off->annotations))
                        {

                            off_mult = -1.0;
                            break;

                        } 

                        continue;

                    }

                    // Retrieving the location to the respective set of probabilities
                    unsigned int loc = map_to_nodes[tmpstate];

                    off_mult *= node.offspring[o]->subtree_prob[loc];

                }

                // Is this state valid?
                if (off_mult < 0.0)
                {

                    ++nstate;
                    continue;
                    
                }

                // Multiplying by P(x|x_n), the transition probability
                off_mult *= model->likelihood(
                    par0,
                    psets_stats->at(nstate++),
                    node.narray[s]
                );

                // Adding to the total probabilities
                totprob += off_mult;

            }

            // Setting the probability at the node
            node.subtree_prob[s] = totprob;

        }

        // All probabilities should be completed at this point
        if (node.parent == nullptr)
        {

            for (unsigned int s = 0u; s < states.size(); ++s)
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