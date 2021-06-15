// #include "geese-bones.hpp"

#ifndef GEESE_MEAT_PREDICT_HPP
#define GEESE_MEAT_PREDICT_HPP 1

inline std::vector< std::vector<double> > Geese::predict_backend(
    const std::vector< double > & par,
    bool use_reduced_sequence,
    const std::vector< uint > & preorder
)
{

    // Splitting the probabilities
    std::vector< double > par_terms(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    // Generating probabilities at the root-level (root state)
    std::vector< double > rootp(this->states.size(), 1.0);
    for (unsigned int i = 0u; i < rootp.size(); ++i)
    {

        for (unsigned int j = 0u; j < nfuns(); ++j)
            rootp[i] *= states[i][j] ? par_root[j] : (1.0 - par_root[j]);
        
    }

    // Making room 
    std::vector< std::vector<double> > res(nnodes());

    // Step 1: Computing the probability at the root node
    std::vector< double > tmp_prob(nfuns(), 0.0);
    unsigned int root_id = preorder[0u];
    Node * tmp_node      = &nodes[root_id];
    tmp_node->probability.resize(states.size(), 0.0);
    double tmp_likelihood = likelihood(par, false, use_reduced_sequence);

    for (unsigned int s = 0u; s < states.size(); ++s) {

        // Overall state probability P(x_s | D)
        tmp_node->probability[s] = tmp_node->subtree_prob[s] * rootp[s] /
            tmp_likelihood;

        // Marginalizing the probabilities P(x_sf | D)
        for (unsigned int f = 0u; f < nfuns(); ++f) {

            // Since the probability, the expected value, is for
            // observing an x = 1, then we need to make sure that we
            // are multiplying by the corresponding state
            if (states[s][f])
                tmp_prob[f] += tmp_node->probability[s];
        }
        

    }

    // Storing the final prob
    res[nodes[preorder[0u]].ord] = tmp_prob;
    
        // Going in the opposite direction
    for (auto& i : preorder)
    {

        Node & node = nodes[i];

        // We just started from the root
        if (node.parent == nullptr)
            continue;

        // Reserving space
        node.probability.resize(states.size(), 0.0);

        std::vector< double > zerovec(nfuns(), 0.0);
        res[node.ord] = zerovec;

        // Need to identify what is the position of the node with respect to
        // its siblings
        unsigned int n_pos = 0u;
        for (unsigned int n = 0u; n < node.parent->offspring.size(); ++n)
        {

            if (node.parent->offspring[n]->id == node.id)
                break;
            
            ++n_pos;
            
        }

        // Iterating through the offspring state P(x_n^p | D)
        std::fill(node.probability.begin(), node.probability.end(), 0.0);
        for (unsigned int s = 0u; s < states.size(); ++s)
        {         

            // Iterating throught the parent state
            for (unsigned int s_p = 0u; s_p < states.size(); ++s_p)
            {

                // All the ways in which we can go from x_pk to x_nk, we first
                // need to find its location in the model (loc):
                unsigned int loc = node.parent->narray[s_p];
                
                // Retrieving the corresponding arrays and stats that will be
                // use to marginalize
                auto p_arrays = model->get_pset(loc);
                auto p_stats  = model->get_stats(loc);

                double prob = 0.0;

                // Iterating through all the cases in which x_p -> x_nk^p
                for (unsigned int a = 0u; a < p_arrays->size(); ++a)
                {
                
                    const auto & A = p_arrays->operator[](a);

                    // Should we include this?
                    bool includeit = true;
                    for (unsigned int k = 0u; k < nfuns(); ++k)
                        if (A(k, n_pos, false) != static_cast<unsigned int>(states[s][k]))
                        {
                            
                            // If it does not match, then jump to the next state
                            includeit = false;
                            break;
                            
                        }
                    

                    // If not to be included, then we go to the next
                    if (!includeit)
                        continue;

                    // Computing the likelihood 
                    prob += model->likelihood(par_terms, p_stats->at(a), loc);
                    
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

        res[node.ord] = tmp_prob;

    }
        
    return res;

}

inline std::vector< std::vector<double> > Geese::predict(
    const std::vector< double > & par,
    std::vector< std::vector< double > > * res_prob,
    bool leave_one_out,
    bool only_annotated,
    bool use_reduced_sequence
)
{

    INITIALIZED()

    // Inverse sequence
    std::vector< unsigned int > preorder;
    if (only_annotated)
        preorder = this->reduced_sequence;
    else
        preorder = this->sequence;

    std::reverse(preorder.begin(), preorder.end());

    // Full prediction (first run, right now I am doing this
    // twice. Need to fix in the future)
    std::vector< std::vector<double> > res = predict_backend(
        par, use_reduced_sequence, preorder
        );

    // If the user requires the probability matrix per state
    if (res_prob != nullptr)
    {

        res_prob->resize(nnodes());
        for (auto& i : sequence)
            res_prob->at(nodes[i].ord) = nodes[i].probability;

    }


    // In this case, we need to update the predictions, mostly of the annotated
    // leaf nodes. Because of 
    if (leave_one_out)
    {

        std::vector< unsigned int > default_empty(nfuns(), 9u);
        for (auto& n : nodes)
        {

            if (n.second.is_leaf()) {

                Node & ntmp = n.second;

                // We only make the changes if it is not all missing
                bool use_it = false;
                for (auto& n_state : ntmp.annotations)
                    if (n_state != 9u)
                    {

                        use_it = true;
                        break;

                    }
                

                if (!use_it)
                    continue;

                // Recording the original annotation
                auto old_ann = ntmp.annotations;

                // Removing the entire gene
                update_annotations(ntmp.id, default_empty);

                // Making the prediction
                res[ntmp.ord] = (
                    predict_backend(par, use_reduced_sequence, preorder)
                )[ntmp.ord];

                // Restoring the gene
                update_annotations(ntmp.id, old_ann);

                if (res_prob != nullptr)
                    res_prob->operator[](ntmp.ord) = ntmp.probability;


            }

        }

    }

    
    return res;

}

#endif