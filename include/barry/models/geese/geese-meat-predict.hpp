// #include "geese-bones.hpp"

#ifndef GEESE_MEAT_PREDICT_HPP
#define GEESE_MEAT_PREDICT_HPP 1

inline std::vector< std::vector<double> > Geese::predict_backend(
    const std::vector< double > & par,
    bool use_reduced_sequence,
    const std::vector< size_t > & preorder
)
{

    // Splitting the probabilities
    std::vector< double > par_terms(par.begin(), par.end() - nfuns());
    std::vector< double > par_root(par.end() - nfuns(), par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    // Generating probabilities at the root-level (root state)
    std::vector< double > rootp(this->states.size(), 1.0);
    for (size_t s = 0u; s < rootp.size(); ++s)
    {

        for (size_t f = 0u; f < nfuns(); ++f)
            rootp[s] *= states[s][f] ? par_root[f] : (1.0 - par_root[f]);
        
    }

    // Making room 
    std::vector< std::vector<double> > res(
        nnodes(), std::vector<double>(nfuns())
        );

    // Step 1: Computing the probability at the root node
    std::vector< double > tmp_prob(nfuns(), 0.0); 
    size_t root_id = preorder[0u];
    Node * tmp_node = &nodes[root_id];
    tmp_node->probability.resize(states.size(), 0.0);
    double tmp_likelihood = likelihood(par, false, use_reduced_sequence);

    if (!std::isfinite(tmp_likelihood))
    {
        throw std::runtime_error("Likelihood is not finite");
    }

    for (size_t s = 0u; s < states.size(); ++s)
    {

        // Overall state probability P(x_s | D)
        tmp_node->probability[s] = tmp_node->subtree_prob[s] * rootp[s] /
            tmp_likelihood;

        if (!std::isfinite(tmp_node->probability[s]))
        {
            throw std::runtime_error("Probability is not finite");
        }
            
            

        

        // Marginalizing the probabilities P(x_sf | D)
        for (size_t f = 0u; f < nfuns(); ++f)
        {

            // Since the probability, the expected value, is for
            // observing an x = 1, then we need to make sure that we
            // are multiplying by the corresponding state
            if (states[s][f])
                tmp_prob[f] += tmp_node->probability[s];

        }
        

    }

    // Storing the final prob
    res[nodes[root_id].ord] = tmp_prob;
    
    // Retrieving the powersets probabilities
    const auto & pset_probs = *(model->get_pset_probs());
    const auto & arrays2support = *(model->get_arrays2support());
    const auto & pset_locations = *(model->get_pset_locations());

    // This will start from the root node and go down
    for (auto & i : preorder)
    {

        // Leafs have nothing to do here
        Node & parent = nodes[i];
        if (parent.is_leaf())
            continue;

        // Creating space.
        std::vector< std::vector< double > > everything_below;
        std::vector< std::vector< double > > everything_above;

        everything_below.reserve(states.size());
        everything_above.reserve(states.size());
        
        // All combinations of the the parent states
        // So psets[s] = combinations of offspring given state s.
        //    psets[s][i] = The ith combination of offspring given state s.
        std::vector< std::vector< PhyloArray > > psets;
        psets.reserve(states.size());

        // Making space for the offspring
        for (auto & off : parent.offspring)
        {
            off->probability.resize(states.size(), 0.0);
            std::fill(off->probability.begin(), off->probability.end(), 0.0);
        }

        // Iterating through the parent states
        for (size_t s = 0u; s < states.size(); ++s)
        {
            std::vector< double > below;
            std::vector< double > above;
            std::vector< PhyloArray > pset;

            // Array id in the support and support id of the array
            size_t array_id   = parent.narray[s];
            size_t support_id = arrays2support[array_id];

            // Retrieving powerset of stats and arrays
            const auto & pset_arrays   = model->get_pset(array_id);

            // Going over all possible combinations given parent is state s
            for (size_t p = 0u; p < pset_arrays->size(); ++p)
            {

                // Corresponding graph and target stats
                const PhyloArray & array_p = pset_arrays->at(p);

                // Adding to the map, we only do this during the first run,
                // afterwards, we need to actually look for the array.
                bool in_the_set = true; /// < True if the array belongs to the set
                
                // Everything below just need to be computed only once
                // and thus, if already added, no need to go through all of this!
                double everything_below_p = 1.0;
                for (size_t off = 0u; off < parent.offspring.size(); ++off)
                {

                    // Below leafs, the everything below is 1.
                    if (parent.offspring[off]->is_leaf())
                    {

                        // But we can only includ it if the current state actually
                        // matches the leaf data (otherwise the prob is 0)
                        const auto & off_ann = parent.offspring[off]->annotations;
                        for (size_t f = 0u; f < nfuns(); ++f)
                        {

                            if ((off_ann[f] != 9u) && (off_ann[f] != array_p(f, off)))
                            {
                                in_the_set = false;
                                break;
                            }
                                
                        }

                        if (!in_the_set)
                            break;

                        continue;

                    } else {

                        // Getting the offspring state, and how it maps, only
                        // if it is not an offspring
                        const auto & off_state = array_p.get_col_vec(off);
                        size_t loc = this->map_to_state_id.find(off_state)->second;

                        everything_below_p *= parent.offspring[off]->subtree_prob[loc];

                    }

                }

                // If it is not in the set, then continue to the next array
                if (!in_the_set)
                    continue;

                pset.push_back(array_p); // Generating a copy
                
                // - With focal node, conditioning on it beening status s.
                // - But the offspring probabilities are the central ones here.
                // - So the saved values are for computing P(x_offspring | Data)
                below.push_back(everything_below_p);

                // The first run, we only need to grow the list
                above.push_back(
                    pset_probs[pset_locations[support_id]]
                    // model->likelihood(
                    //     par_terms, target_p, parent.narray[s], false
                    // )
                    *  parent.probability[s] / parent.subtree_prob[s]
                );


            } // end for psets

            // Storing the psets
            psets.push_back(std::move(pset));
            everything_below.push_back(std::move(below));
            everything_above.push_back(std::move(above));

            
        } // end for states

        // Marginalizing at the state level for each offspring
        for (size_t s = 0u; s < states.size(); ++s)
        {

            for (size_t p = 0u; p < everything_above[s].size(); ++p)
            {

                // p-th pset
                const auto & pset_p = psets[s][p];

                // Updating the probability (it is the product)
                everything_above[s][p] *= everything_below[s][p];

                for (size_t off = 0u; off < parent.offspring.size(); ++off)
                {

                    // Figuring out the state of the offspring
                    auto cvec = pset_p.get_col_vec(off);
                    size_t off_s = this->map_to_state_id[cvec];
                    parent.offspring[off]->probability[off_s] += everything_above[s][p];

                    // We integrate over the offspring itsefl
                    for (size_t f = 0u; f < nfuns(); ++f)
                    {
                        if (cvec[f] == 1u)
                            res[parent.offspring[off]->ord][f] += everything_above[s][p];
                    }
                    


                }

            }
        }

        // Finally, we can marginalize the values at the 
        // gene function level.
        for (const auto & off : parent.offspring)
        {

            // Checking that probabilities add up to one
            for (size_t f = 0u; f < nfuns(); ++f)
            {
                if ((res[off->ord][f] > 1.00001) || (res[off->ord][f] < -.0000009))
                {
                    auto msg = "[geese] Out-of-range probability for node.id " +
                        std::to_string(off->id) + " for function " +
                        std::to_string(f) + ": " +
                        std::to_string(res[off->ord][f]);

                    throw std::logic_error(msg);
                    
                } 

                if (res[off->ord][f] > 1.0)
                    res[off->ord][f] = 1.0;
                else if (res[off->ord][f] < 0.0)
                    res[off->ord][f] = 0.0;

            }

        }

    } // end for over preorder
        
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
    std::vector< size_t > preorder;
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

    std::vector< double > par0(par.begin(), par.end() - nfuns());
    model->update_pset_probs(par0, 1u);

    // In this case, we need to update the predictions, mostly of the annotated
    // leaf nodes. Because of 
    if (leave_one_out)
    {

        std::vector< size_t > default_empty(nfuns(), 9u);
        for (auto& n : nodes)
        {

            if (n.second.is_leaf())
            {

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