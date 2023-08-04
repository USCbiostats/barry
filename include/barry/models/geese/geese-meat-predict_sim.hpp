// #include "geese-bones.hpp"

#ifndef GEESE_MEAT_PREDICT_SIM_HPP
#define GEESE_MEAT_PREDICT_SIM_HPP 1

inline std::vector< std::vector<double> > Geese::predict_sim(
    const std::vector< double > & par,
    bool use_reduced_sequence,
    size_t nsims
)
{

    INITIALIZED()

    // Preparing
    std::vector< std::vector< size_t > > tmp;

    std::vector< double > zerovec(nfuns(), 0.0);
    std::vector< std::vector< double > > res_vec(nnodes(), zerovec);
    std::vector< int > counts(nnodes(), 0);

    // We will iterate through this list everytime we need to check
    // whether we have all the annotations for the conditional prob.
    auto annotated = this->get_annotated_nodes();

    for (size_t i = 0u; i < nsims; ++i)
    {

        // Generating a sample
        tmp = this->simulate(par);

        for (auto j = nodes.begin(); j != nodes.end(); ++j)
        {
            // Retrieving node
            const Node & n = j->second;

            // Checking we have all matching
            bool includeit = true;
            for (auto & id : annotated)
            {

                // Same node need not to match (since we are not conditionin
                // each node on itself!)
                if (n.id == id)
                    continue;

                const auto & ord     = nodes[id].ord; 
                const auto & n_w_ann = nodes[id].annotations;
                for (size_t f = 0u; f < nfuns(); ++f)
                {
                    // No checking missings
                    if (n_w_ann[f] == 9u)
                        continue;

                    // If this is not matching, then we cannot use it!
                    if (n_w_ann[f] != tmp[ord][f])
                    {
                        includeit = false;
                        break;
                    }

                }

                if (!includeit)
                    break;
            }

            // If it passed the test, then we can use it for counting stuff
            if (!includeit)
                continue;

            for (size_t f = 0u; f < nfuns(); ++f)
                if (tmp[n.ord][f] == 1)
                    ++res_vec[n.ord][f];

            ++counts[n.ord];

        }

    }

    // Once the simulations have finalized, we can then approximate
    // probabilities
    for (size_t i = 0u; i < nnodes(); ++i)
    {
        // printf_barry("We used %i counts for node %i.\n", counts[i], i);
        for (size_t f = 0u; f < nfuns(); ++f)
            res_vec[i][f] /= (static_cast< double >(counts[i]) + 1e-10);
    }
    
    return res_vec;

}


#endif