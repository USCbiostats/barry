
#ifndef GEESE_MEAT_LIKELIHOOD_EXHAUST_HPP
#define GEESE_MEAT_LIKELIHOOD_EXHAUST_HPP 1
// #include "../../barry.hpp"
// #include "geese-bones.hpp" 

inline double Geese::likelihood_exhaust(
    const std::vector< double > & par
)
{

    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    // This is only worthwhile if the number of nodes is small
    if (this->nnodes() > 6)
        throw std::overflow_error("Too many nodes! Exhaust calculation of likelihood cannot be done for such cases.");

    if (this->nfuns() > 3)
        throw std::overflow_error("Too many functions! Exhaust calculation of likelihood cannot be done for such cases.");

    // Computing all combinations ----------------------------------------------
    PhyloArray base(nfuns(), nnodes());
    for (auto& n : nodes)
    {

        for (size_t i = 0u; i < nfuns(); ++i)
            base(i, n.second.ord) = n.second.annotations[i];
            
    }

    PhyloPowerSet pset(base);//this->nfuns(), this->nnodes());
    pset.add_rule(
            rule_empty_free<PhyloArray,PhyloRuleData>,
            PhyloRuleData()
            );
    pset.calc();

    // Inverse sequence
    std::vector< size_t > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    double totprob = 0.0;
    
    // This vector says whether the probability has to be included in 
    // the final likelihood or not.
    for (size_t p = 0u; p < pset.size(); ++p)
    {
        
        // ith state
        const PhyloArray * s = &pset[p];
        
        // Following the sequence
        double prob = 1.0;
        std::vector< size_t > tmpstates(this->nfuns());

        Node * node;
        for (auto& i : preorder)
        {

            node = &nodes[i];
            std::fill(tmpstates.begin(), tmpstates.end(), 0u);
            s->get_col_vec(&tmpstates, node->ord, false);

            // Root node first
            if (node->parent == nullptr)
            {               
                // Since it is the root, the first probability is computed using
                // the root only
                for (auto k = 0u; k < this->nfuns(); ++k)
                    prob *= tmpstates[k] == 1u ? par_root[k] : (1.0 - par_root[k]);

            }
            else if (node->is_leaf())
                continue;

            // Computing the transition
            PhyloArray transition(nfuns(), node->offspring.size());

            std::vector< double > bl(node->offspring.size(), 1.0);

            std::vector< bool > sl = vector_caster<bool,size_t>(tmpstates);

            transition.set_data(
                new NodeData(bl, sl, node->duplication),
                true
            );

            // Filling the array
            for (size_t a = 0u; a < nfuns(); ++a)
            {

                for (size_t o = 0u; o < node->offspring.size(); ++o)
                {

                    if (s->get_cell(a, node->offspring[o]->id) == 1u)
                        transition(a, o) = 1u;

                }

            }

            prob *= this->model->likelihood(
                par0,
                transition,
                node->narray[this->map_to_state_id[tmpstates]],
                false
                );

        }

        totprob += prob;
    }

    return totprob;

}
#endif