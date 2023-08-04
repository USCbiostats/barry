#include "tests.hpp"
#include "../include/barry/models/geese.hpp"

BARRY_TEST_CASE("Geese model with rules", "[geese with rules]") {
// int main() {

    using namespace barry::counters::phylo;

    // More interesting experiment
    std::vector< std::vector<size_t> > ann = {
        {9, 9},
        {9, 9},
        {9, 0},
        {1, 1},
        {0, 9}
    };
    // ann.resize(3u, {9, 9});

    // std::vector< size_t > geneid = {0, 1, 2, 3, 4};
    // std::vector< int >  parent = {4, 3, 3, 4, -1};
    std::vector< size_t > geneid = {0, 1, 2, 3, 4};
    std::vector< int >  parent = {-1, 0, 0, 1, 1};

    std::vector< bool > duplication(geneid.size(), true);

    Geese model(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model.get_counters(), {0, 1});
    counter_maxfuns(model.get_counters(), 2, 2);
    counter_overall_changes(model.get_counters());

    rule_dyn_limit_changes(model.get_support_fun(), 0, 0, 1);

    model.init();
    model.set_seed(100);

    // Computing the likelihood

    model.print();
    std::vector< double > params = {0,0,0,0,0,0};
    auto x = model.likelihood(params, true, true);
    printf_barry("Likelihood: %f\n", x);

    // A model with a rule avoiding gains
    Geese model2(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model2.get_counters(), {0, 1});
    counter_loss(model2.get_counters(), {0, 1});

    rule_dyn_limit_changes(model2.get_support_fun(), 0, 0, 0);

    model2.init();
    model2.set_seed(100);

    // Computing the likelihood

    model2.print();
    std::vector< double > params2 = {0, 0, 0, 0, 0, 0};
    auto x2 = model2.likelihood(params2, true, true);
    printf_barry("Likelihood: %f\n", x2);

    #ifndef CATCH_CONFIG_MAIN
    return 0;
    #endif

}
