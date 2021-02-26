#include "aphylomodel-bones.hpp"
// #include <algorithm>
// #include <random>

#ifndef APHYLOMODEL_MEAT_SIMULATE_HPP
#define APHYLOMODEL_MEAT_SIMULATE_HPP 1

void APhyloModel::set_seed(const unsigned int & s) {
    rengine.seed(s);
}

template<typename Ta, typename Tb>
inline std::vector< Ta > vector_caster(const std::vector< Tb > & x) {
    std::vector< Ta > ans;
    ans.reserve(x.size());
    for (auto i = x.begin(); i != x.end(); ++i)
        ans.push_back(static_cast< Ta >(*i));
    return ans;
}

std::vector< std::vector< unsigned int > > APhyloModel::simulate(
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
    std::vector< std::vector< unsigned int > > res(nodes.size());

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
    res[nodes[preorder[0u]].id] =
        vector_caster< unsigned int, bool>(states[idx]);

    // Going in the opposite direction
    for (auto& i : preorder) {

        if (nodes[i].is_leaf())
            continue;

        // Getting the state of the node      
        unsigned int n = this->map_to_nodes[res[nodes[i].id]];

        // Given the state of the current node, sample the state of the
        // offspring, all based on the current state
        auto tmp = model_full.sample(nodes[i].idx_full[n], par0);

        // Iterating through the offspring to assign the state
        unsigned int m;
        for (unsigned int j = 0u; j < nodes[i].offspring.size(); ++j) {
            res[nodes[i].offspring[j]->id] = tmp.get_col_vec(j, false);
        }

    }

    return res;

}

#endif