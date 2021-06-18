#ifndef GEESE_MEAT_SIMULATE_HPP
#define GEESE_MEAT_SIMULATE_HPP 1

inline void Geese::set_seed(const unsigned int & s) {
    rengine->seed(s);
}

inline std::vector< std::vector< unsigned int > > Geese::simulate(
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
    std::vector< std::vector< unsigned int > > res(nodes.size());

    // Inverse sequence
    std::vector< unsigned int > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    // Generating probabilities at the root-level (root state)
    std::vector< double > rootp(states.size(), 1.0);
    for (unsigned int i = 0u; i < rootp.size(); ++i)
    {
        for (unsigned int j = 0u; j < nfuns(); ++j)
            rootp[i] *= states[i][j] ? par_root[j] : (1.0 - par_root[j]);
    }

    // Preparing the random number generator
    std::uniform_real_distribution<> urand(0, 1);
    double r         = urand(*rengine);
    unsigned int idx = 0u;
    double cumprob = rootp[idx];
    while ((idx < rootp.size()) && (cumprob < r))
    {
        cumprob += rootp[++idx];
    }
    
    // We now know the state of the root
    res[nodes[preorder[0u]].ord] =
        vector_caster< unsigned int, bool>(states[idx]);

    // Going in the opposite direction
    for (auto& i : preorder)
    {

        if (nodes[i].is_leaf())
            continue;

        const Node & n = nodes[i];

        // Getting the state of the node      
        unsigned int l = map_to_nodes[res[n.ord]];

        // Given the state of the current node, sample the state of the
        // offspring, all based on the current state
        auto z = n.narray;
        auto tmp = model->sample(n.narray[l], par0);

        // Iterating through the offspring to assign the state
        for (unsigned int j = 0u; j < n.offspring.size(); ++j)
            res[n.offspring[j]->ord] = tmp.get_col_vec(j, false);

    }

    return res;

}

#endif