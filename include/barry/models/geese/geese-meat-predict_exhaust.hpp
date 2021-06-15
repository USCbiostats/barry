
#ifndef GEESE_MEAT_PREDICT_EXHAUST_HPP
#define GEESE_MEAT_PREDICT_EXHAUST_HPP 1
// #include "../../barry.hpp"
// #include "geese-bones.hpp" 

inline std::vector< std::vector<double> > Geese::predict_exhaust(
    const std::vector< double > & par,
    bool only_annotated,
    bool use_reduced_sequence
) {

    INITIALIZED()

    // This is only worthwhile if the number of nodes is small
    if (this->nnodes() > 6)
        throw std::overflow_error("Too many nodes! Exhaust calculation of likelihood cannot be done for such cases.");

    if (this->nfuns() > 2)
        throw std::overflow_error("Too many functions! Exhaust calculation of prediction cannot be done for such cases.");

    // Computing all combinations ----------------------------------------------
    // The base PhyloArray will store the original set of annotations.
    phylocounters::PhyloArray base(nfuns(), nnodes());
    for (auto& n : nodes) {
        for (unsigned int i = 0u; i < nfuns(); ++i)
            base(i, n.second.ord) = n.second.annotations[i];
    }

    phylocounters::PhyloPowerSet pset(base);//this->nfuns(), this->nnodes());
    pset.add_rule(
            rule_empty_free<phylocounters::PhyloArray,phylocounters::PhyloRuleData>
            );
    pset.calc();

    // Making space for the expected values
    std::vector< double > expected(nnodes() * nfuns(), 0.0);
    
    // This vector says whether the probability has to be included in 
    // the final likelihood or not.
    for (unsigned int p = 0u; p < pset.size(); ++p) {
        
        // ith state
        const phylocounters::PhyloArray * s = &pset[p];

        // We now need to update the annotations to match those in the
        // powerset.
        for (auto & n: nodes)
            this->update_annotations(n.second.id, s->get_col_vec(n.second.ord));

        double prob = this->likelihood(par, false, use_reduced_sequence);
        std::cout << "Prob: " << prob << std::endl;
        // Adding to the overall probability
        for (auto & n: nodes)
            for (unsigned int j = 0u; j < nfuns(); ++j)
                expected[n.second.ord +  j * nnodes()] += n.second.annotations[j] * prob;
        
    }

    // Returning to the original set of annotations
    for (auto & n: nodes)
        this->update_annotations(n.second.id, base.get_col_vec(n.second.ord));

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