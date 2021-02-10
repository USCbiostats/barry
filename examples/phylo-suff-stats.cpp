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

using namespace phylocounters;

int main() {

    // Defining input data
    std::vector< unsigned int > row = {0, 1, 0, 1, 2};
    std::vector< unsigned int > col = {0, 0, 1, 1, 1};
    std::vector< unsigned int > val = {1, 1, 9, 1, 9};

    // Creating the node object
    PhyloArray node0(3, 2, row, col, val);

    // Preparing the rule data
    PhyloRuleData locked_blocks = {{0,1}, {2,1}};
    PhyloRule rule(
        rule_blocked<PhyloArray, PhyloRuleData>,
        &locked_blocks, false
        );

    PhyloPowerSet pset(3, 1);
    pset.calc();

    // Each node will have its own model.
    std::vector< PhyloArray > nodes(0u);
    std::vector< PhyloModel > models(0u);
    std::vector< PhyloModel > models_trunk(0u);
    std::vector< double > blengths(3, 1.0);

    // Baseline parameters
    std::vector< double > params = {0.1, -.1, 0, 0, 0};

    // Iterating through the possible sets of nodes
    unsigned int i = 0u;
    for (auto iter = pset.begin(); iter != pset.end(); ++iter) {

        std::vector< bool > states(blengths.size(), false);
        for (auto iter2 = iter->get_col(0u)->begin(); iter2 != iter->get_col(0u)->end(); ++iter2) {
            states.at(iter2->first) = true;
        }

        // Generating the corresponding node data
        nodes.push_back(node0);
        nodes.at(i).set_data(
            new NodeData(blengths, states),
            true
        );

        // Generating the model
        models.push_back(PhyloModel());
        models_trunk.push_back(PhyloModel());

        counter_gains(&(models.at(i).counters), {0,1,2});
        counter_cogain(&(models.at(i).counters), 0, 1);
        counter_cogain(&(models.at(i).counters), 1, 2);

        counter_gains(&(models_trunk.at(i).counters), {0,1,2});
        counter_cogain(&(models_trunk.at(i).counters), 0, 1);
        counter_cogain(&(models_trunk.at(i).counters), 1, 2);

        // Only one has the rules
        models.at(i).add_rule(rule);

        // Running the support
        models.at(i).add_array(nodes.at(i));
        models_trunk.at(i).add_array(nodes.at(i));

        // Computing the likelihoods
        models.at(i).likelihood_total(params);
        models_trunk.at(i).likelihood_total(params);

        // Printing the data
        double norm_num = models.at(i).normalizing_constants.at(0u);
        double norm_denom = models_trunk.at(i++).normalizing_constants.at(0u);
        printf("Probability of observing this data: %.8f\n", norm_num/norm_denom);

    }

    
    // /* Setting the following array:
    //  *   o1 o2
    //  * a  1  9  
    //  * b  1  1
    //  * c  0  9
    //  */
    // // node0(0,0) = true;
    // // node0(1,0) = true;
    // // node0(1,1) = true;

    // NodeData dat({1.0,1.0,1.0}, {false, false, false});
    // node0.set_data(&dat, false);
  
    // // Model
    // PhyloModel model;
    // counter_gains(&model.counters, {0,1,2});
    // counter_cogain(&model.counters, 0, 1);
    // counter_cogain(&model.counters, 1, 2);

    // // Need to block a few cells
    // PhyloRuleData locked_blocks = {{0,1}, {2,1}};
    // PhyloRule rule(
    //     rule_blocked<PhyloArray, PhyloRuleData>,
    //     &locked_blocks, false
    //     );

    // model.add_rule(rule);

    // // Running the support
    // model.add_array(node0);
    // model.print_stats(0);

    // std::vector< double > params = {0, 0, 0, 0, 0};

    // model.likelihood_total(params);
    // printf("Normalizing constant: %.4f\n", model.normalizing_constants.at(0));

    // // Computing the probabilities
    // double norm0 = model.normalizing_constants.at(0);

    // PhyloModel model2;
    // counter_gains(&model2.counters, {0,1,2});
    // counter_cogain(&model2.counters, 0, 1);
    // counter_cogain(&model2.counters, 1, 0);
    
    // model2.add_array(node0);
    // model2.print_stats(0);
    
    // model2.likelihood_total(params);
    // double norm1 = model2.normalizing_constants.at(0);
    // printf("Normalizing constant: %.4f\n", model2.normalizing_constants.at(0));

    // printf("Probability of observing this data: %.8f\n", norm0/norm1);

    return 0;
}