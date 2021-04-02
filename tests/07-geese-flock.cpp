#include "../include/barry/models/geese.hpp"

TEST_CASE("Flock likelihood", "[flock-likelihood]") {

    // More interesting experiment
    std::vector< std::vector<uint> > ann = {
        {9, 9},
        {9, 9},
        {9, 0},
        {1, 1},
        {0, 9}
    };
    // ann.resize(3u, {9, 9});

    // std::vector< uint > geneid = {0, 1, 2, 3, 4};
    // std::vector< int >  parent = {4, 3, 3, 4, -1};
    std::vector< uint > geneid = {0, 1, 2, 3, 4};
    std::vector< int >  parent = {-1, 0, 0, 1, 1};

    std::vector< bool > duplication(geneid.size(), true);

    Geese model(ann, geneid, parent, duplication);

    // Adding terms
    phylocounters::counter_gains(model.counters, {0, 1});
    phylocounters::counter_maxfuns(model.counters, 2, 2);

    model.init();
    model.set_seed(100);

    // Tryingout likelihood
    std::vector<double> params = {1, 1, -1, -10, -10};

    double ans0 = std::log(model.likelihood(params)) * 2.0;
    
    // Generating a flock ------------------------------------------------------
    Flock aflock;
    aflock.add_data(ann, geneid, parent, duplication);
    aflock.add_data(ann, geneid, parent, duplication);

    // Adding terms
    phylocounters::counter_gains(aflock.counters_ptr(), {0, 1});
    phylocounters::counter_maxfuns(aflock.counters_ptr(), 2, 2);

    aflock.init();

    double ans1 = aflock.likelihood_joint(params, true);
    REQUIRE(std::abs(ans0 - ans1) < .00000001);

    // Checking the likelihood sequence -----------------------------------------
    std::vector< std::vector<uint> > ann2 = {
        {9, 9},
        {9, 9},
        {9, 9}, // Now both are missing, before {9, 0}
        {1, 1},
        {0, 9}
    };

    Geese model2(ann2, geneid, parent, duplication);

    // Adding terms
    phylocounters::counter_gains(model2.counters, {0, 1});
    phylocounters::counter_maxfuns(model2.counters, 2, 2);

    model2.init();

    ans0 = model2.likelihood(params, true);
    ans1 = model2.likelihood(params, false);

    REQUIRE(std::abs(ans0 - ans1) < .00000001);

}