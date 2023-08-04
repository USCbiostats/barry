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
    Node * tmp_node      = &nodes[root_id];
    tmp_node->probability.resize(states.size(), 0.0);
    double tmp_likelihood = likelihood(par, false, use_reduced_sequence);

    for (size_t s = 0u; s < states.size(); ++s)
    {

        // This is mostly relevant when the root is connected
        // to a leaf node, where the array may not be valid 
        // according to the dynamic rules, i.e., gains/loses.
        // If 2 genes have a function, the root state is zero,
        // and the rule is no more than one gain, then the
        // state is not valid.
        if (!tmp_node->arrays_valid[s])
        {
            tmp_node->probability[s] = 0.0;
            continue;
        }

        // Overall state probability P(x_s | D)
        tmp_node->probability[s] = tmp_node->subtree_prob[s] * rootp[s] /
            tmp_likelihood;

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
    res[nodes[preorder[0u]].ord] = tmp_prob;
    size_t n_pars = par_terms.size();

    for (auto & i : preorder)
    {

        // Leafs have nothing to do here
        Node & parent = nodes[i];
        if (parent.is_leaf())
            continue;

        // Creating space.
        std::vector< std::vector< double > > everything_below(states.size());
        std::vector< std::vector< double > > everything_above(states.size());
        std::vector< std::vector< phylocounters::PhyloArray > > psets(states.size());

        // Making space for the offspring
        for (auto & off : parent.offspring)
        {
            off->probability.resize(states.size(), 0.0);
            std::fill(off->probability.begin(), off->probability.end(), 0.0);
        }

        // Iterating through the parent states
        for (size_t s = 0u; s < states.size(); ++s)
        {

            // Only if it is valid: See the same block applied to the root
            if (!parent.arrays_valid[s])
            {
                parent.probability[s] = 0.0;
                continue;
            }

            // Retrieving powerset of stats and arrays
            // it is not const since we will flip the states back and forth
            // to generate the key
            const auto & pset_arrays = model->get_pset(parent.narray[s]);
            const std::vector<double> * pset_target = model->get_pset_stats(parent.narray[s]);

            for (size_t p = 0u; p < pset_arrays->size(); ++p)
            {

                // Corresponding graph and target stats
                const phylocounters::PhyloArray & array_p = pset_arrays->at(p);
                std::vector<double> target_p(n_pars, 0.0);
                for (size_t par_i = 0u; par_i < target_p.size(); ++par_i)
                    target_p[par_i] = pset_target->operator[](p * n_pars + par_i);

                phylocounters::PhyloArray tmp_array(nfuns(), array_p.ncol());
                tmp_array += array_p;

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
                        size_t loc = this->map_to_nodes[off_state];

                        everything_below_p *= parent.offspring[off]->subtree_prob[loc];

                    }

                }

                // If it is not in the set, then continue to the next array
                if (!in_the_set)
                    continue;

                psets[s].push_back(array_p); // Generating a copy
                everything_below[s].push_back(everything_below_p);

                // The first run, we only need to grow the list
                everything_above[s].push_back(
                    model->likelihood(
                        par_terms, target_p, parent.narray[s], false
                    ) *  parent.probability[s] / parent.subtree_prob[s]
                );


            } // end for psets
            
        } // end for states

        // Marginalizing at the state level
        for (size_t s = 0u; s < states.size(); ++s)
        {

            // Only if it is valid
            if (!parent.arrays_valid[s])
                continue;

            for (size_t p = 0u; p < everything_above[s].size(); ++p)
            {

                // p-th pset
                const auto & pset_p = psets[s][p];

                // Updating the probability (it is the product)
                everything_above[s][p] *= everything_below[s][p];

                for (size_t off = 0u; off < parent.offspring.size(); ++off)
                {

                    // Figuring out the state of the offspring
                    size_t off_s = this->map_to_nodes[pset_p.get_col_vec(off)];
                    parent.offspring[off]->probability[off_s] += everything_above[s][p];


                }

            }
        }

        // Finally, we can marginalize the values at the 
        // gene function level.
        for (const auto & off : parent.offspring)
        {
            for (size_t s = 0u; s < states.size(); ++s)
            {

                // Only if it is valid
                if (!parent.arrays_valid[s])
                    continue;

                for (size_t f = 0u; f < nfuns(); ++f)
                    if (states[s][f]) 
                        res[off->ord][f] += off->probability[s];

            }

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


    // In this case, we need to update the predictions, mostly of the annotated
    // leaf nodes. Because of 
    if (leave_one_out)
    {

        std::vector< size_t > default_empty(nfuns(), 9u);
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