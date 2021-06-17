
#ifndef GEESE_MEAT_PREDICT_EXHAUST_HPP
#define GEESE_MEAT_PREDICT_EXHAUST_HPP 1
// #include "../../barry.hpp"
// #include "geese-bones.hpp" 

/*

    // Creating the phyloarray, nfunctions x noffspring
    n.array = phylocounters::PhyloArray(nfunctions, n.offspring.size());
    std::vector< bool > tmp_state = vector_caster<bool,uint>(n.annotations);
    std::vector< double > blen(n.offspring.size(), 1.0);
    n.array.set_data(
        new phylocounters::NodeData(blen, tmp_state, n.duplication),
        true
    );

    // We initialize all with a zero since, if excluded from the pruning process,
    // We need to set it to one (as the result of the full integration).
    n.subtree_prob.resize(states.size(), 1.0);

    // Adding the data, first through functions
    for (unsigned int k = 0u; k < nfunctions; ++k) {

        // Then through the offspring
        unsigned int j = 0;
        for (auto& o : n.offspring) {

            // If leaf, then it may have an annotation
            if (o->is_leaf()) {
                if (o->annotations.at(k) != 0) {
                    n.array.insert_cell(
                        k, j, o->annotations.at(k), false, false
                        );
                }
            } else {
                // Otherwise, we fill it with a 9.
                n.array.insert_cell(
                    k, j, 9u, false, false
                    );

            }

            ++j;

        }

    }

    // We then need to set the powerset
    if (n.arrays.size() != states.size()) {
        n.arrays.resize(states.size());
        n.narray.resize(states.size());
    }
    
    for (unsigned int s = 0u; s < states.size(); ++s) {

        n.arrays[s] = phylocounters::PhyloArray(n.array, true);
        n.arrays[s].set_data(
            new phylocounters::NodeData(blen, states[s], n.duplication),
            true
        );

        // Once the array is ready, we can add it to the model
        n.narray[s] = model->add_array(n.arrays[s]);

    }

    return;
*/

inline std::vector< std::vector<double> > Geese::predict_exhaust(
    const std::vector< double > & par,
    bool only_annotated
) {

    INITIALIZED()

    // This is only worthwhile if the number of nodes is small
    if (this->nnodes() > 6)
        throw std::overflow_error("Too many nodes! Exhaust calculation of likelihood cannot be done for such cases.");

    if (this->nfuns() > 2)
        throw std::overflow_error("Too many functions! Exhaust calculation of prediction cannot be done for such cases.");

    // Processing the probabilities --------------------------------------------
    std::vector< double > par_terms(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    // Computing all combinations ----------------------------------------------
    // The base PhyloArray will store the original set of annotations.
    phylocounters::PhyloArray base(nfuns(), nnodes());
    for (auto& n : nodes)
    {

        for (unsigned int i = 0u; i < nfuns(); ++i)
            base(i, n.second.ord) = n.second.annotations[i];

    }

    phylocounters::PhyloPowerSet pset(base);//this->nfuns(), this->nnodes());
    pset.add_rule(
            rule_empty_free<phylocounters::PhyloArray,phylocounters::PhyloRuleData>
            );
    pset.calc();
    
    // Generating the sequence preorder sequence -------------------------------
    std::vector< unsigned int > preorder;
    if (only_annotated)
        preorder = this->reduced_sequence;
    else
        preorder = this->sequence;

    std::reverse(preorder.begin(), preorder.end());
    
    // Making space for the expected values
    std::vector< double > expected(nnodes() * nfuns(), 0.0);
    
    // This vector says whether the probability has to be included in 
    // the final likelihood or not.
    for (unsigned int p = 0u; p < pset.size(); ++p)
    {
        
        // ith state
        const phylocounters::PhyloArray * s = &pset[p];
        s->print("Start pset %02i\n", p);
        // Computing the likelihood of the state s        
        double current_prob = 1.0;
        for (auto & o: preorder)
        {
            // Getting the corresponding node
            Node & n = nodes[o];

            // Nothing to do at the leaf level (leafs are calculated from parents)
            if (n.is_leaf())
                continue;

            // Extracting the parent column (without checking boundaries)
            auto par_state = s->get_col_vec(n.ord, false);

            // Need to compute the root probability (if we havent before)
            if (n.parent == nullptr)
            {

                for (unsigned int f = 0u; f < nfuns(); ++f)
                    current_prob *= par_state[f] ? par_root[f] : (1.0 - par_root[f]);

            }
        
            // Generating a copy of the observed array
            // (data is copied so that we can chage the state of the parent)
            phylocounters::PhyloArray tmparray(n.array, true);

            // Updating the state of the parent
            for (unsigned int f = 0u; f < nfuns(); ++f)
                tmparray.D()->states[f] = par_state[f] == 1u;

            // Updating offspring annotations
            int loc = 0;
            for (auto & off : n.offspring) {
                
                for (unsigned int f = 0u; f < nfuns(); ++f)
                {

                    if (s->operator()(f, off->ord) == 1u)
                        tmparray(f, loc) = 1u;
                    else
                        tmparray.rm_cell(f, loc);

                }

                // Next offspring start in the next column of the array, Duh.
                ++loc;

            }
            tmparray.print("Array of node %02i\n", n.ord);
            // Computing the likelihood
            current_prob *= model->likelihood(par_terms, tmparray, -1, false);

        }
            // this->update_annotations(n.second.id, s->get_col_vec(n.second.ord));
        
        // Adding to the overall probability
        for (auto & n: nodes)
            for (unsigned int j = 0u; j < nfuns(); ++j)
                expected[n.second.ord +  j * nnodes()] += s->operator()(j, n.second.ord) * current_prob;
        
    }

    // Coercing expected to a list vector
    std::vector< std::vector< double > > res(nnodes());
    std::vector< double > zerovec(nfuns(), 0.0);
    for (auto & n: nodes)
    {
        res[n.second.ord] = zerovec;
        for (unsigned int i = 0u; i < nfuns(); ++i)
            res[n.second.ord][i] = expected[n.second.ord +  i * nnodes()];
    }

    return res;

}
#endif