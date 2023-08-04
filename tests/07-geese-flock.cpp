#ifndef CATCH_CONFIG_MAIN
// #include "/opt/intel/oneapi/advisor/2022.0.0/include/advisor-annotate.h"
#include <chrono>
// #define INTEL_ADVISOR
#endif

#include "tests.hpp"

#include "../include/barry/models/geese.hpp"

BARRY_TEST_CASE("Flock likelihood", "[flock-likelihood]") {

    using namespace phylocounters;

    // More interesting experiment
    std::vector< std::vector<size_t> > ann = {
        {9, 9, 1, 0},
        {9, 9, 1, 0},
        {9, 0, 1, 0},
        {1, 1, 1, 0},
        {0, 9, 1, 0}
    };
    // ann.resize(3u, {9, 9});

    // std::vector< size_t > geneid = {0, 1, 2, 3, 4};
    // std::vector< int >  parent = {4, 3, 3, 4, -1};
    std::vector< size_t > geneid = {0, 1, 2, 3, 4};
    std::vector< int >  parent = {-1, 0, 0, 1, 1};

    std::vector< bool > duplication(geneid.size(), true);

    Geese model(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model.get_counters(), {0, 1}, 1);
    counter_gains(model.get_counters(), {0, 1}, 0);
    counter_maxfuns(model.get_counters(), 2, 2, 1);
    counter_maxfuns(model.get_counters(), 2, 2, 0);
    counter_prop_genes_changing(model.get_counters(), 1);
    counter_prop_genes_changing(model.get_counters(), 0);

    model.init();
    model.set_seed(100);

    // Tryingout likelihood
    std::vector<double> params = {1, 1, 1, 1, -1, .5, .5, -1, -10, -10, -10, -10};

    auto simex = model.simulate(params);

    double ans0 = std::log(model.likelihood(params)) * 2.0;
    
    // Generating a flock ------------------------------------------------------
    Flock aflock;
    aflock.add_data(ann, geneid, parent, duplication);
    aflock.add_data(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(aflock.get_counters(), {0, 1}, 1);
    counter_gains(aflock.get_counters(), {0, 1}, 0);
    counter_maxfuns(aflock.get_counters(), 2, 2, 1);
    counter_maxfuns(aflock.get_counters(), 2, 2, 0);
    counter_prop_genes_changing(aflock.get_counters(), 1);
    counter_prop_genes_changing(aflock.get_counters(), 0);

    aflock.init();

    double ans1 = aflock.likelihood_joint(params, true);
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::abs(ans0 - ans1) < .00000001);
    #endif

    // Checking the likelihood sequence -----------------------------------------
    geneid = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    parent = {-1, 0,0, 1, 1, 2, 2, 5, 5};
    ann    = {
        {9, 9, 0, 1},
        {9, 9, 0, 1},
        {9, 9, 0, 1},
        {0, 1, 0, 1},
        {1, 9, 0, 1},
        {9, 9, 0, 1},
        {0, 1, 0, 1},
        {9, 9, 0, 1},
        {9, 9, 0, 1}
    };

    std::vector<bool> duplication2(geneid.size(), true);
    
    Geese model2(ann, geneid, parent, duplication2);

    // Adding terms
    counter_gains(model2.get_counters(), {0, 1}, 1);
    counter_gains(model2.get_counters(), {0, 1}, 0);
    counter_maxfuns(model2.get_counters(), 2, 2, 1);
    counter_maxfuns(model2.get_counters(), 2, 2, 0);
    counter_prop_genes_changing(model2.get_counters(), 1);
    counter_prop_genes_changing(model2.get_counters(), 0);

    model2.init();

    ans0 = model2.likelihood(params, true, true);
    ans1 = model2.likelihood(params, true, false);

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::abs(ans0 - ans1) < .00000001);
    #endif

    // Measuring time
    auto start = std::chrono::system_clock::now();
    // for (size_t i = 0u; i < 500; ++i)
    //     model2.likelihood(params, true);
    auto end = std::chrono::system_clock::now();
    // std::chrono::duration<double> diff0 = end-start;

    std::uniform_real_distribution<> rand(-1, 1);
    std::mt19937 ren;

    // start = std::chrono::system_clock::now();
    // for (size_t i = 0u; i < 500; ++i)
    //     model2.likelihood(params, false);
    // end = std::chrono::system_clock::now();
    // std::chrono::duration<double> diff1 = end-start;

    // std::cout << "Total time (full seq): " << diff1.count() << std::endl;
    // std::cout << "Total time (trun seq): " << diff0.count() << std::endl;

    start = std::chrono::system_clock::now();
    for (size_t i = 0u; i < 5000; ++i)
    {
        params[std::floor(std::fabs(rand(ren) * params.size()))] = rand(ren) * 2;
        aflock.likelihood_joint(params, true);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff_updating_normseq = end-start;

    std::cout << "Total time (normseq): " << diff_updating_normseq.count() << std::endl;
    

}