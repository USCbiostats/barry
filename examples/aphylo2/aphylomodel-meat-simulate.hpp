#include "aphylomodel-bones.hpp"
#include <algorithm>
#include <random>

#ifndef APHYLOMODEL_MEAT_SIMULATE_HPP
#define APHYLOMODEL_MEAT_SIMULATE_HPP 1

void APhyloModel::set_seed(const unsigned int & s) {
    rengine.seed(s);
}

std::vector< std::vector< bool > > APhyloModel::simulate(
    const std::vector< double > & par
    ) {

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfuns);
    std::vector< double > par_root(par.end() - nfuns, par.end());

    // Scaling root
    for (auto& p : par_root) {
        p = std::exp(p)/(std::exp(p) + 1);
    }

    // Making room 
    std::vector< std::vector< bool > > res(nodes.size());

    // Inverse sequence
    std::vector< unsigned int > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    // Generating probabilities at the root-level (root state)
    std::vector< double > rootp(states.size(), 1.0);
    for (unsigned int i = 0u; i < rootp.size(); ++i) {
        for (unsigned int j = 0u; j < nfuns; ++j) {
            rootp[i] *= states[i][j] ? par_root[j] : (1.0 - par_root[j]);
        }
    }

    // Preparing the random number generator
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(rengine);
    unsigned int idx = 0u;
    double cumprob = rootp[idx];
    while ((idx < rootp.size()) && (cumprob < r)) {
        cumprob += rootp[++idx];
    }
    
    // We now know the state of the root
    res[nodes[preorder[0u]].id] = states[idx];

    // Going in the opposite direction
    for (auto& i : preorder) {

        // Given the state of the current node, sample the state of the
        // offspring, all based on the current state
        auto tmp = model_full.sample(nodes[i].idx_full[idx], par0);

    }

    return res;

}

#endif