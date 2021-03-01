#include "aphylomodel-bones.hpp"

#ifndef APHYLOMODEL_MEAT_LIKELIHOOD_HPP
#define APHYLOMODEL_MEAT_LIKELIHOOD_HPP

double APhyloModel::likelihood(const std::vector< double > & par) {

    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root) {
        p = std::exp(p)/(std::exp(p) + 1);
    }

    std::vector< unsigned int > tmpstate(nfunctions);

    double ll = 0.0;

    // Following the prunning sequence
    for (auto& i : this->sequence) {

        // We cannot compute probability at the leaf, we need to continue
        if (this->nodes[i].is_leaf())
            continue;

        // If it is located right before a leaf, then the computations
        // are a bit different
        if (this->nodes[i].offspring.at(0u)->is_leaf()) {

            // Iterating through the different parent states
            for (unsigned int s = 0u; s < states.size(); ++s) {

                // Update the normalizing constants
                double numer = model_const.get_norm_const(par0, nodes[i].idx_cons[s]);
                double denom = model_full.get_norm_const(par0, nodes[i].idx_full[s]);

                // Computing the probability at "leaf" level
                nodes[i].probabilities[s] = numer/denom;

            }

        } else {

            // Iterating through states
            for (unsigned int s = 0u; s < states.size(); ++s) {

                // Starting the prob
                double totprob = 0.0;

                // Retrieving the sets of arrays
                const std::vector< phylocounters::PhyloArray > * psets = model_full.get_pset(
                    nodes[i].idx_full[s]
                );

                const std::vector< std::vector<double> > * psets_stats = model_full.get_stats(
                    nodes[i].idx_full[s]
                );

                // Summation over all possible values of X
                unsigned int nstate = 0u;
                for (auto x = psets->begin(); x != psets->end(); ++x) {

                    // Extracting the possible values of each offspring
                    double off_mult = 1.0;
                    for (auto o = 0u; o < x->ncol(); ++o) {

                        // First, getting what is the corresponding state
                        std::fill(tmpstate.begin(), tmpstate.end(), 0u);
                        x->get_col_vec(&tmpstate, o, false);

                        // Retrieving the location to the respective set of probabilities
                        unsigned int loc = nodes[i].offspring[o]->idx_full[map_to_nodes[tmpstate]];
                        off_mult *= nodes[i].offspring[o]->probabilities[loc];

                    }

                    // Multiplying by P(x|x_n)
                    off_mult *= model_full.likelihood(
                        par0,
                        psets_stats->at(nstate++),
                        nodes[i].idx_full[s]
                    );

                    // Adding to the total probabilities
                    totprob += off_mult;

                }

                // Setting the probability at the node
                nodes[i].probabilities[s] = totprob;

            }

            // All probabilities should be completed at this point
            if (nodes[i].parent == nullptr) {
                for (unsigned int s = 0u; s < states.size(); ++s) {
                    double tmpll = 1.0;
                    for (auto k = 0u; k < nfunctions; ++k) {
                        tmpll *= states[s][k] ? par_root[k] : (1 - par_root[k]);
                    }

                    ll += tmpll * nodes[i].probabilities[s];

                }
            }



        }


    }

    return ll;

}
#endif