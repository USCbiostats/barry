#ifndef GEESE_MEAT_SIMULATE_HPP
#define GEESE_MEAT_SIMULATE_HPP 1

inline void Geese::set_seed(const size_t & s) {
    rengine->seed(s);
}

inline std::vector< std::vector< size_t > > Geese::simulate(
    const std::vector< double > & par
    ) {

    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root) {
        p = std::exp(p)/(std::exp(p) + 1);
    }

    // Making room 
    std::vector< std::vector< size_t > > res(nodes.size());

    // Inverse sequence
    std::vector< size_t > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    // Generating probabilities at the root-level (root state)
    std::vector< double > rootp(states.size(), 1.0);
    for (size_t i = 0u; i < rootp.size(); ++i)
    {

        for (size_t j = 0u; j < nfuns(); ++j)
            rootp[i] *= states[i][j] ? par_root[j] : (1.0 - par_root[j]);

    }

    // Preparing the random number generator
    std::uniform_real_distribution<> urand(0, 1);
    double r         = urand(*rengine);
    size_t idx = 0u;
    double cumprob = rootp[idx];
    while ((idx < rootp.size()) && (cumprob < r))
    {
        cumprob += rootp[++idx];
    }

    #ifdef BARRY_DEBUG
    
    // auto totprob = std::accumulate(rootp.begin(), rootp.end(), 0.0);
    // if (totprob < 0.9999999999999999 || totprob > 1.0000000000000001)
    //     throw std::runtime_error("Root probabilities do not sum to 1!"
    //         " (totprob = " + std::to_string(totprob) + ")");
    
    #endif
    
    // We now know the state of the root
    res[nodes[preorder[0u]].ord] =
        vector_caster< size_t, bool>(states[idx]);

    // Going in the opposite direction
    for (auto& i : preorder)
    {

        if (nodes[i].is_leaf())
            continue;

        const Node & n = nodes[i];

        // Getting the id of the state
        size_t lth_state = map_to_state_id[res[n.ord]];

        // Given the state of the current node, sample the state of the
        // offspring, all based on the current state
        // auto z = n.narray;
        auto tmp = model->sample(n.narray[lth_state], par0);

        // Iterating through the offspring to assign the state
        for (size_t j = 0u; j < n.offspring.size(); ++j)
            res[n.offspring[j]->ord] = tmp.get_col_vec(j, false);

    }

    return res;

}

#endif