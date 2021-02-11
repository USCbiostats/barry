#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include "../include/barry/barry.hpp"

// The same need to be locked
RULE_FUNCTION(rule_blocked) {
    if (Array->get_cell(i, j) == 9u)
        return false;
    return true;
};

using namespace phylocounters;

template <typename T>
using Vec = std::vector< T >;

template<typename T1, typename T2>
Vec< T1 > caster(const Vec< T2 > & vec) {

    Vec< T1 > ans;
    ans.reserve(vec.size());

    for (auto &i : vec) {
        ans.push_back(*i);
    }

    return ans;

}

// Hasher
inline Vec< double > tip_keygen(const PhyloArray & array) {
    
    // Baseline data: nrows and columns
    Vec< double > dat = {
        (double) array.nrow(), (double) array.ncol()
    };
    
    // State of the parent
    for (bool i : array.data->states) {
        dat.push_back(i ? 1.0 : 0.0);
    }

    // Free cells
    for (auto i = 0u; i < array.nrow(); ++i)
        for (auto j = 0u; j < array.ncol(); ++j)
            dat.push_back((double) array.get_cell(i, j, false));
    
    return dat;
}

// Hasher
inline Vec< double > tip_keygen_baseline(const PhyloArray & array) {
    
    // Baseline data: nrows and columns
    Vec< double > dat = {
        (double) array.nrow(), (double) array.ncol()
    };
    
    // State of the parent
    for (bool i : array.data->states) {
        dat.push_back(i ? 1.0 : 0.0);
    }

    // // Free cells
    // for (auto i = 0u; i < array.nrow(); ++i)
    //     for (auto j = 0u; j < array.ncol(); ++j)
    //         dat.push_back((double) array.get_cell(i, j, false));
    
    return dat;
}


int main() {

    // Defining input data
    /* Setting the following array:
     *   o1 o2   o3 o4
     * a  1  9    0  1
     * b  1  1    0  0
     * c  0  9    9  9
     */
    Vec< uint > row = {0, 1, 0, 1, 2};
    Vec< uint > col = {0, 0, 1, 1, 1};
    Vec< uint > val = {1, 1, 9, 1, 9};

    // Stargint to measure time
    auto start = std::chrono::system_clock::now();

    // Creating the node object
    PhyloArray node0(3, 2, row, col, val);

    PhyloPowerSet pset(3, 1);
    pset.calc();

    // Each node will have its own model.
    Vec< PhyloArray > nodes(0u);
    Vec< double > blengths(3, 1.0);

    // Baseline parameters
    Vec< double > params = {0.1, -.1, 0, 0, 0};

    // Iterating through the possible sets of nodes
    PhyloModel tipmodel;
    tipmodel.set_keygen(tip_keygen_baseline);

    PhyloModel tipmodel_restricted;
    tipmodel_restricted.set_keygen(tip_keygen);
    tipmodel_restricted.add_rule(rule_blocked<PhyloArray,PhyloRuleData>);

    counter_gains(&tipmodel.counters, {0,1,2});
    counter_cogain(&tipmodel.counters, 0, 1);
    counter_cogain(&tipmodel.counters, 1, 2);

    tipmodel_restricted.set_counters(&tipmodel.counters);


    unsigned int i = 0u;
    for (auto iter = pset.begin(); iter != pset.end(); ++iter) {

        Vec< bool > states(blengths.size(), false);
        for (auto iter2 = iter->get_col(0u)->begin(); iter2 != iter->get_col(0u)->end(); ++iter2) {
            states.at(iter2->first) = true;
        }

        // Generating the corresponding node data
        nodes.push_back(node0);
        nodes.at(i).set_data(
            new NodeData(blengths, states),
            true
        );

        // Running the support
        uint l = tipmodel_restricted.add_array(nodes.at(i));
        uint k = tipmodel.add_array(nodes.at(i));

        // Computing the likelihoods
        tipmodel_restricted.likelihood_total(params);
        tipmodel.likelihood_total(params);

        // Printing the data
        double norm_num = tipmodel_restricted.normalizing_constants.at(l);
        double norm_denom = tipmodel.normalizing_constants.at(k);
        
        printf("(l, k): (%i, %i); ", tipmodel_restricted.arrays2support.at(l), tipmodel.arrays2support.at(k));
        printf("Pr(): %.8f\n", norm_num/norm_denom);

    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;

    std::cout << "Total time: " << diff.count() << std::endl;
    
    // Measuring time updating the data
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(123); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-1, 1);

    std::chrono::duration<double> mdiff;

    for (auto i = 0u; i < 1000; ++i) {

        // Updating the parameter
        Vec<double> par(params.size());
        for (auto j = 0u; j < par.size(); ++j)
            par.at(j) = dis(gen);

        // Computing the likelihoods
        auto s = std::chrono::system_clock::now();
        tipmodel_restricted.likelihood_total(par);
        tipmodel.likelihood_total(par);
        auto e = std::chrono::system_clock::now();

        if (!(i % 100)) {
            double norm_num = tipmodel_restricted.normalizing_constants.at(0u);
            double norm_denom = tipmodel.normalizing_constants.at(0u);
            printf("%5i Pr(): %.8f\n", i, norm_num/norm_denom);
        }

        if (i == 0u)
            mdiff = e-s;
        else
            mdiff += (e-s);
    }

    std::cout << "Timing 1000 updates: " << mdiff.count() << std::endl;

    return 0;
}