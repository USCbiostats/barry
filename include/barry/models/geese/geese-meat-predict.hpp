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
    std::vector< double > par_terms(par.begin(), par.end() - nfuns());
    std::vector< double > par_root(par.end() - nfuns(), par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    // Generating probabilities at the root-level (root state)
    std::vector< double > rootp(this->states.size(), 1.0);
    for (unsigned int s = 0u; s < rootp.size(); ++s)
    {

        for (unsigned int f = 0u; f < nfuns(); ++f)
            rootp[s] *= states[s][f] ? par_root[f] : (1.0 - par_root[f]);
        
    }

    // Making room 
    std::vector< std::vector<double> > res(
        nnodes(), std::vector<double>(nfuns()));

    // Step 1: Computing the probability at the root node
    std::vector< double > tmp_prob(nfuns(), 0.0);
    unsigned int root_id = preorder[0u];
    Node * tmp_node      = &nodes[root_id];
    tmp_node->probability.resize(states.size(), 0.0);
    double tmp_likelihood = likelihood(par, false, use_reduced_sequence);

    for (unsigned int s = 0u; s < states.size(); ++s)
    {

        // Overall state probability P(x_s | D)
        tmp_node->probability[s] = tmp_node->subtree_prob[s] * rootp[s] /
            tmp_likelihood;

        // Marginalizing the probabilities P(x_sf | D)
        for (unsigned int f = 0u; f < nfuns(); ++f)
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

    // In wat comes next, we need to store all the possible combination
    // values at the event level. We have (# Internal nodes) events, and
    // each one of them has (# offspring possible states)
    std::vector< std::vector< double > > predictions(nnodes());
    
        // Going in the opposite direction
    for (auto& i : preorder)
    {

        printf_barry("Looking at node %i\n", i);

        Node & node = nodes[i];

        // We just started from the root
        if (node.parent == nullptr)
            continue;
        else if (node.is_leaf())
            continue;

        // Reserving space
        for (const auto & off : node.offspring)
            predictions[off->ord].resize(nfuns(), 0.0);

        // We start by computing the "Everything below"
        // Since at this point only matters the state of the offspring,
        // we just grab the first of narray.
        const auto & pset = model->get_pset(node.narray[0u]);
        std::vector< double > Prob_Xoff_given_D(pset->size(), 1.0);
        for (unsigned int p = 0u; p < pset->size(); ++p)
        {

            const auto & pset_p = pset->operator[](p);

            // Everything below Xoff
            for (unsigned int off = 0u; off < node.offspring.size(); ++off)
            {

                // Below leafs, the everything below is 1.
                if (node.offspring[off]->is_leaf())
                    continue;

                // Getting the offspring state, and how it maps
                const auto & off_state = pset_p.get_col_vec(off);
                unsigned int loc = this->map_to_nodes[off_state];

                Prob_Xoff_given_D[p] *= node.offspring[off]->subtree_prob[loc];

            }
            
        }

        // We now procede to compute everything above. For this, we don't
        // need to iterate through the offspring, but rather the parent's 
        // possible states
        std::vector< double > Prob_Everything_Above(pset->size(), 0.0);
        for (unsigned int s = 0u; s < states.size(); ++s)
        {

            // Retrieving the corresponding state of the parent
            unsigned int narray_s = node.narray[s];
            const auto & pset_s   = model->get_pset(narray_s);
            const auto & stats_s  = model->get_stats(narray_s);

            for (unsigned int p = 0u; p < pset_s->size(); ++p)
            {

                double tmp_prob = 1.0;

                // Comnputing transitional prob
                tmp_prob *= model->likelihood(
                    par_terms, stats_s->at(p), narray_s
                    );

                // Times the likelihood of observing the parent
                // in that state
                tmp_prob *= node.probability[s];

                // Divided by the induced sub tree probability
                tmp_prob /= node.subtree_prob[s];

                Prob_Everything_Above[p] += tmp_prob;
                
            }

        }

        // Now, computing the final steps, i.e., multiplying
        // everything above times everything below, and marginalizing
        // the values at the gene level
        for (unsigned int p = 0u; p < pset->size(); ++p)
        {

            const auto & pset_p = pset->operator[](p);
            Prob_Xoff_given_D[p] *= Prob_Everything_Above[p];

            for (unsigned int off = 0u; off < node.offspring.size(); ++off)
            {

                // Figuring out the state of the offspring
                unsigned int off_s = this->map_to_nodes[pset_p.get_col_vec(off)];
                node.offspring[off]->probability[off_s] += Prob_Xoff_given_D[p];

            }

        }

        // Finally, we can marginalize the values at the 
        // gene function level.
        for (const auto & off : node.offspring)
        {
            for (unsigned int s = 0u; s < states.size(); ++s)
            {

                for (unsigned int f = 0u; f < nfuns(); ++f)
                    if (states[s][f])
                        res[off->ord][f] += off->probability[s];

            }
                
                
        }
        

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