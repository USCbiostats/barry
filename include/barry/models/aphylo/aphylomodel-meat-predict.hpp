#include "aphylomodel-bones.hpp"

#ifndef APHYLOMODEL_MEAT_PREDICT_HPP
#define APHYLOMODEL_MEAT_PREDICT_HPP 1


inline std::vector< std::vector<double> > APhyloModel::predict(
    const std::vector< double > & par
    ) {

    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root) {
        p = std::exp(p)/(std::exp(p) + 1);
    }

    // Making room 
    std::vector< std::vector<double> > res(nnodes());

    // Inverse sequence
    std::vector< unsigned int > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    // Generating probabilities at the root-level (root state)
    std::vector< double > rootp(states.size(), 1.0);
    for (unsigned int i = 0u; i < rootp.size(); ++i) {
        for (unsigned int j = 0u; j < nfuns(); ++j) {
            rootp[i] *= states[i][j] ? par_root[j] : (1.0 - par_root[j]);
        }
    }

    // Step 1: Computing the probability at the root node
    std::vector< double > tmp_prob(nfuns(), 0.0);
    unsigned int root_id = preorder[0u];
    Node * tmp_node = &nodes[root_id];
    double tmp_likelihood = likelihood(par);
    for (unsigned int s = 0u; s < states.size(); ++s) {

        // Overall state probability P(x_s | D)
        tmp_node->probability[s] = tmp_node->subtree_prob[s] * rootp[s] /
            tmp_likelihood;

        // Marginalizing the probabilities P(x_sf | D)
        for (unsigned int f = 0u; f < nfuns(); ++f) {
            if (states[s][f])
                tmp_prob[f] += tmp_node->probability[s];
        }
        

    }

    // Storing the final prob
    res[nodes[preorder[0u]].id] = tmp_prob;
    
        // Going in the opposite direction
    for (auto& i : preorder) {

        // We just started from the root
        if (nodes[i].parent == nullptr)
            continue;

        std::vector< double > zerovec(nfuns(), 0.0);
        res[nodes[i].id] = zerovec;

        // Iterating through the parent state
        for (unsigned int s = 0u; s < states.size(); ++s) {

            // Resetting the prob vec
            std::fill(tmp_prob.begin(), tmp_prob.end(), 0.0);

            // All the ways in which we can go from x_pk to x_nk, we first
            // need to find its location in the model (loc):
            unsigned int loc = nodes[i].parent->narray[s];
            
            // Iterating through the possible staes of the parent
            for (unsigned int sp = 0u; sp < states.size(); ++sp) {

                // tmp_prob[]

            }

            // Retrieving the corresponding arrays and stats that will be
            // use to marginalize
            auto p_arrays = model.get_pset(loc);
            auto p_stats  = model.get_stats(loc);
            for (unsigned int a = 0u; a < p_arrays->size(); ++a) {

                // Iterating through offspring
                for (unsigned int o = 0u; o < p_arrays->size(); ++o) {
                
                    // Iterating through function
                    for (unsigned int f = 0u; f < nfuns(); ++f) {

                        if (p_arrays->at(a).get_cell(f, o, false) != 0u) {
                            tmp_prob[f] += model.likelihood(par0, p_stats->at(a), loc);
                        }

                    }
                }

            }

            for (unsigned int f = 0u; f < nfuns(); ++f) {
                res[nodes[i].id][f] += nodes[i].parent->subtree_prob[s] * 
                    nodes[i].parent->subtree_prob[s];
            }
            
            

        }



    }

    return res;

}

#endif