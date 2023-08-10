
#ifndef GEESE_MEAT_PREDICT_EXHAUST_HPP
#define GEESE_MEAT_PREDICT_EXHAUST_HPP 1

inline std::vector< std::vector<double> > Geese::predict_exhaust(
    const std::vector< double > & par
) {

    INITIALIZED()

    // This is only worthwhile if the number of nodes is small
    if (this->nnodes() > 6)
        throw std::overflow_error("Too many nodes! Exhaust calculation of likelihood cannot be done for such cases.");

    if (this->nfuns() > 2)
        throw std::overflow_error("Too many functions! Exhaust calculation of prediction cannot be done for such cases.");


    // Generating the sequence preorder sequence -------------------------------
    std::vector< size_t > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    std::vector< std::vector< double > > res = predict_exhaust_backend(
        par, preorder
        );

    // Looping to do LOO
    std::vector< size_t > annotated_ids = this->get_annotated_nodes();
    std::vector< size_t > missing_vec(nfuns(), 9u);
    for (auto & i : annotated_ids) {

        Node & n = nodes[i];

        auto old_ann = n.annotations;
        update_annotations(i, missing_vec);

        res[n.ord] = predict_exhaust_backend(par, preorder)[n.ord];

        update_annotations(i, old_ann);

    }

    return res;

}

inline std::vector< std::vector<double> > Geese::predict_exhaust_backend(

    const std::vector< double > & par,
    const std::vector< size_t > & preorder
) {

    // Processing the probabilities --------------------------------------------
    std::vector< double > par_terms(par.begin(), par.end() - nfuns());
    std::vector< double > par_root(par.end() - nfuns(), par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    double baseline_likelihood = this->likelihood(par);

    // Computing all combinations ----------------------------------------------
    // The base PhyloArray will store the original set of annotations.
    PhyloArray base(nfuns(), nnodes());
    for (auto& n : nodes)
    {

        for (size_t f = 0u; f < nfuns(); ++f)
            base(f, n.second.ord) = n.second.annotations[f];

    }

    PhyloPowerSet pset(base);//this->nfuns(), this->nnodes());
    pset.add_rule(
            rule_empty_free<PhyloArray,PhyloRuleData>,
            PhyloRuleData()
            );
    pset.calc();
    
    // Making space for the expected values
    std::vector< double > expected(nnodes() * nfuns(), 0.0);
    
    // This vector says whether the probability has to be included in 
    // the final likelihood or not.
    for (size_t p = 0u; p < pset.size(); ++p)
    {
        
        // ith state
        const PhyloArray * s = &pset[p];
        
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

                for (size_t f = 0u; f < nfuns(); ++f)
                    current_prob *= par_state[f] ? par_root[f] : (1.0 - par_root[f]);

            }
        
            // Generating a copy of the observed array
            // (data is copied so that we can chage the state of the parent)
            PhyloArray tmparray(n.array, true);

            // Updating the state of the parent
            for (size_t f = 0u; f < nfuns(); ++f)
                tmparray.D_ptr()->states[f] = par_state[f] == 1u;

            // Updating offspring annotations
            int loc = 0;
            for (auto & off : n.offspring) {
                
                for (size_t f = 0u; f < nfuns(); ++f)
                {

                    if (s->operator()(f, off->ord) == 1u)
                        tmparray(f, loc) = 1u;
                    else
                        tmparray.rm_cell(f, loc);

                }

                // Next offspring start in the next column of the array, Duh.
                ++loc;

            }
            
            // Computing the likelihood
            current_prob *= model->likelihood(par_terms, tmparray, -1, false);

        }
            // this->update_annotations(n.second.id, s->get_col_vec(n.second.ord));
        
        // Adding to the overall probability
        for (auto & n: nodes)
            for (size_t j = 0u; j < nfuns(); ++j)
                expected[n.second.ord +  j * nnodes()] += s->operator()(j, n.second.ord) * current_prob/
                    baseline_likelihood;
        
    }

    // Coercing expected to a list vector
    std::vector< std::vector< double > > res(nnodes());
    std::vector< double > zerovec(nfuns(), 0.0);
    for (auto & n: nodes)
    {
        res[n.second.ord] = zerovec;
        for (size_t i = 0u; i < nfuns(); ++i)
            res[n.second.ord][i] = expected[n.second.ord +  i * nnodes()];
    }

    return res;

}
#endif