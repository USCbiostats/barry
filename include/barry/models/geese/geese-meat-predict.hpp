#include "geese-bones.hpp"

#ifndef GEESE_MEAT_PREDICT_HPP
#define GEESE_MEAT_PREDICT_HPP 1


inline std::vector< std::vector<double> > Geese::predict(
    const std::vector< double > & par,
    std::vector< std::vector< double > > * res_prob
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
    tmp_node->probability.resize(states.size(), 0.0);
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

        Node & node = nodes[i];

        // We just started from the root
        if (node.parent == nullptr)
            continue;

        // Reserving space
        node.probability.resize(states.size(), 0.0);

        std::vector< double > zerovec(nfuns(), 0.0);
        res[node.id] = zerovec;

        // Need to identify what is the position of the node with respect to
        // its siblings
        unsigned int n_pos = 0u;
        for (unsigned int n = 0u; n < node.parent->offspring.size(); ++n) {

            if (node.parent->offspring.at(n)->id == node.id)
                break;
            
            ++n_pos;
        }

        // Iterating through the offspring state P(x_n^p | D)
        std::fill(node.probability.begin(), node.probability.end(), 0.0);
        for (unsigned int s = 0u; s < states.size(); ++s) {         

            // Iterating throught the parent state
            for (unsigned int s_p = 0u; s_p < states.size(); ++s_p) {

                // All the ways in which we can go from x_pk to x_nk, we first
                // need to find its location in the model (loc):
                unsigned int loc = node.parent->narray[s_p];
                
                // Retrieving the corresponding arrays and stats that will be
                // use to marginalize
                auto p_arrays = support->get_pset(loc);
                auto p_stats  = support->get_stats(loc);

                double prob = 0.0;

                // Iterating through all the cases in which x_p -> x_nk^p
                for (unsigned int a = 0u; a < p_arrays->size(); ++a) {
                
                    // Should we include this?
                    bool includeit = true;
                    for (unsigned int k = 0u; k < nfuns(); ++k)
                        if (p_arrays->at(a).get_cell(k, n_pos, false) != static_cast<unsigned int>(states[s][k])) {
                            
                            // If it does not match, then jump to the next state
                            includeit = false;
                            break;
                            
                        }
                    

                    // If not to be included, then we go to the next
                    if (!includeit)
                        continue;

                    // Computing the likelihood 
                    prob += support->likelihood(par0, p_stats->at(a), loc);
                    
                }

                // Finalizing
                node.probability[s] += node.parent->probability[s_p] * prob;
            }
            
            

        }

        // Computing marginal probabilities. For this we need to integrate out
        // function by function.
        std::fill(tmp_prob.begin(), tmp_prob.end(), 0.0);
        for (unsigned int s = 0u; s < states.size(); ++s) {
        
            for (unsigned int k = 0u; k < nfuns(); ++k) {

                // If the state is true, then include it, otherwise, don't
                if (states[s][k])
                    tmp_prob[k] += node.probability[s];

            }
                        
        }

        res[node.id] = tmp_prob;

    }

    // If the user requires the probability matrix per state
    if (res_prob != nullptr) {
        res_prob->resize(nnodes());
        for (auto& i : sequence) {
            res_prob->at(nodes[i].id) = nodes[i].probability;
        }
    }
        

    return res;

}

#endif