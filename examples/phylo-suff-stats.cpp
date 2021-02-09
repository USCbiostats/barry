#include <iostream>
#include <string>
#include "../include/barry/barry.hpp"

// The same need to be locked
RULE_FUNCTION(rule_blocked) {
    for (auto iter = data->begin(); iter != data->end(); ++iter)
        if (iter->first == i & iter->second == j)
            return false;
    // if ((i == 2 && j == 1) | (i == 0 && j == 1))
    //     return false;
    return true;
};

int main() {

    // Defining input data
    std::vector< unsigned int > row = {0, 1, 0, 1, 2};
    std::vector< unsigned int > col = {0, 0, 1, 1, 1};
    std::vector< unsigned int > val = {1, 1, 9, 1, 9};

    // Creating the node object
    phylocounters::PhyloArray node0(3, 2, row, col, val);
    
    /* Setting the following array:
     *   o1 o2
     * a  1  9  
     * b  1  1
     * c  0  9
     */
    // node0(0,0) = true;
    // node0(1,0) = true;
    // node0(1,1) = true;

    phylocounters::NodeData dat({1.0,1.0,1.0}, {false, false, false});
    node0.set_data(&dat, false);
  
    // Model
    phylocounters::PhyloModel model;
    phylocounters::counter_gains(&model.counters, {0,1,2});
    phylocounters::counter_cogain(&model.counters, 0, 1);
    phylocounters::counter_cogain(&model.counters, 1, 0);

    // Need to block a few cells
    phylocounters::PhyloRuleData free_blocks = {{0,1}, {2,1}};
    phylocounters::PhyloRule rule(
        rule_blocked<phylocounters::PhyloArray, phylocounters::PhyloRuleData>,
        &free_blocks, false
        );

    model.add_rule(rule);

    // Running the support
    model.add_array(node0);
    model.print_stats(0);

    model.likelihood_total({0, 0, 0, 0, 0});
    printf("Normalizing constant: %.4f\n", model.normalizing_constants.at(0));

    // Computing the probabilities
    double norm0 = model.normalizing_constants.at(0);

    phylocounters::PhyloModel model2;
    phylocounters::counter_gains(&model2.counters, {0,1,2});
    phylocounters::counter_cogain(&model2.counters, 0, 1);
    phylocounters::counter_cogain(&model2.counters, 1, 0);
    
    model2.add_array(node0);
    model2.print_stats(0);
    
    model2.likelihood_total({0, 0, 0, 0, 0});
    double norm1 = model2.normalizing_constants.at(0);
    printf("Normalizing constant: %.4f\n", model2.normalizing_constants.at(0));

    printf("Probability of observing this data: %.8f\n", norm0/norm1);

    return 0;
}