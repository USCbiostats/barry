
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

    // Computing all combinations
    barry::PowerSet<> pset(this->nfuns(), this->nnodes());
    pset.calc();

    // Reserving space
    std::vector< double > probs(pset.size(), 0.0);

    // Inverse sequence
    std::vector< unsigned int > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    double totprob = 0.0;
    std::vector< double > nodeprobs(this->nnodes(), 0.0);
    for (unsigned int p = 0u; p < pset.size(); ++p) {
        
        // ith state
        barry::BArray<> s = pset[p];
        
        // Following the sequence
        double prob = 1.0;
        std::fill(nodeprobs.begin(), nodeprobs.end(), 0.0);
        std::vector< bool > tmpstates(this->nfuns());
        for (auto& i : preorder) {

            // Is it the first node?
            if (nodes[i].parent == nullptr) {

                // Getting the state, and the corresponding index
                tmpstates = s.get_col_vec(nodes[i].id);
                unsigned int m = this->map_to_nodes[tmpstates];

            }


        }
    }

    return 0.0;

}
#endif