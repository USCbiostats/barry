#include "tests.hpp"
#include "../include/barry/models/geese.hpp"

BARRY_TEST_CASE("Phylo model likelihood", "[phylo likelihood]")
{

    using namespace barry::counters::phylo;

    // More interesting experiment
    std::vector< std::vector<uint> > ann = {
        {9, 9, 9},
        {9, 9, 9},
        {9, 0, 0},
        {1, 1, 1},
        {0, 9, 9}
    };
    // ann.resize(3u, {9, 9});

    // std::vector< uint > geneid = {0, 1, 2, 3, 4};
    // std::vector< int >  parent = {4, 3, 3, 4, -1};
    std::vector< uint > geneid = {0, 1, 2, 3, 4};
    std::vector< int >  parent = {-1, 0, 0, 1, 1};

    std::vector< bool > duplication(geneid.size(), true);

    Geese model(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model.get_counters(), {0, 1});
    counter_maxfuns(model.get_counters(), 2, 2);
    counter_prop_genes_changing(model.get_counters(), 1);
    counter_prop_genes_changing(model.get_counters(), 0);
    counter_k_genes_changing(model.get_counters(), 1, 0);
    counter_k_genes_changing(model.get_counters(), 1, 1);

    model.init();
    model.set_seed(100);

    // Tryingout likelihood
    std::vector<double> params = {1, 1, -1, -5, 1, 1, 1, 1, -1, -2.5};

    double ans0a = model.likelihood(params); // 0.057089548117690694
    double ans1a = model.likelihood_exhaust(params); // 0.057089548117690736
    
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::abs(ans0a - ans1a) < .0000001);
    #endif

    // Updating an annotation --------------------------------------------------
    ann[3u][1u] = 9u;
    ann[4u][1u] = 0u;
    Geese model2(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model2.get_counters(), {0, 1});
    counter_maxfuns(model2.get_counters(), 2, 2);
    counter_prop_genes_changing(model2.get_counters(), 1);
    counter_prop_genes_changing(model2.get_counters(), 0);
    counter_k_genes_changing(model2.get_counters(), 1, 0);
    counter_k_genes_changing(model2.get_counters(), 1, 1);

    model2.init();
    model2.set_seed(100);

    // Modifying the annotations
    model.update_annotations(3u, {1u, 9u, 1u});
    model.update_annotations(4u, {0u, 0u, 9u});

    double ans0b = model.likelihood(params);  // 0.056110910038269665
    double ans1b = model2.likelihood(params); // 0.056110910038269665

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::abs(ans0b - ans1b) < .0000001);
    #endif

    // Going back
    model.update_annotations(3u, {1u, 1u, 1u});
    model.update_annotations(4u, {0u, 9u, 9u});

    double ans0c = model.likelihood(params); // 0.057089548117690694

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::abs(ans0a - ans0c) < .0000001);
    #endif

    std::cout << "Test for updating annotations" << std::endl;
    std::cout << "Likelihood(baseline): " << ans0a << std::endl;
    std::cout << "Likelihood(changed) : " << ans0b << std::endl;

}