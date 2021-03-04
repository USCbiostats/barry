
#ifndef APHYLOMODEL_MEAT_LIKELIHOOD_EXACT_HPP
#define APHYLOMODEL_MEAT_LIKELIHOOD_EXACT_HPP 1
#include "../../barry.hpp"
#include "aphylomodel-bones.hpp" 

double APhyloModel::likelihood_exact(const std::vector< double > & par) {

    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root) {
        p = std::exp(p)/(std::exp(p) + 1);
    }

    // This is only worthwhile if the number of nodes is small
    if (this->nnodes() > 6)
        throw std::overflow_error("Too many nodes! Exact likelihood cannot be computed for such cases.");

    if (this->nfuns() > 2)
        throw std::overflow_error("Too many functions! Exact likelihood cannot be computed for such cases.");

    // Computing all combinations ----------------------------------------------
    phylocounters::PhyloArray base(nfuns(), nnodes());
    for (auto& n : nodes) {
        for (unsigned int i = 0u; i < nfuns(); ++i)
            base(i, n.second.id) = n.second.annotations[i];
    }

    phylocounters::PhyloPowerSet pset(base);//this->nfuns(), this->nnodes());
    pset.add_rule(
            rule_empty_free<phylocounters::PhyloArray,phylocounters::PhyloRuleData>
            );
    pset.calc();

    // Inverse sequence
    std::vector< unsigned int > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    double totprob = 0.0;
    
    // This vector says whether the probability has to be included in 
    // the final likelihood or not.
    for (unsigned int p = 0u; p < pset.size(); ++p) {
        
        // ith state
        const phylocounters::PhyloArray * s = &pset[p];
        
        // Following the sequence
        double prob = 1.0;
        std::vector< unsigned int > tmpstates(this->nfuns());

        Node * node;
        for (auto& i : preorder) {

            node = &nodes[i];
            std::fill(tmpstates.begin(), tmpstates.end(), 0u);
            s->get_col_vec(&tmpstates, node->id, false);

            // Root node first
            if (node->parent == nullptr) {               
                // Since it is the root, the first probability is computed using
                // the root only
                for (auto k = 0u; k < this->nfuns(); ++k) {
                    prob *= tmpstates[k] == 1u ? par_root[k] : (1.0 - par_root[k]);
                }

            } else if (node->is_leaf())
                continue;

            // Computing the transition
            phylocounters::PhyloArray transition(nfuns(), node->offspring.size());
            std::vector< double > bl(node->offspring.size(), 1.0);
            std::vector< bool > sl = caster<bool,unsigned int>(tmpstates);
            transition.set_data(
                new phylocounters::NodeData(bl, sl),
                true
            );

            // Filling the array
            for (unsigned int a = 0u; a < nfuns(); ++a) {

                for (unsigned int o = 0u; o < node->offspring.size(); ++o) {
                    if (s->get_cell(a, node->offspring[o]->id) == 1u)
                        transition(a, o) = 1u;
                }

            }

            prob *= this->model_full.likelihood(
                par0, transition,
                node->idx_full[this->map_to_nodes[tmpstates]],
                false);

        }

        totprob += prob;
    }

    return totprob;

}
#endif