#include "../include/barry/models/geese.hpp"

TEST_CASE("Geese model prediction", "[geese prediction]") {

    using namespace barry::counters::phylo;

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
    counter_gains(model.get_counters(), {0, 1});
    counter_maxfuns(model.get_counters(), 2, 2);

    model.init();
    model.set_seed(100);

    // Tryingout likelihood
    std::vector<double> params = {1, 1, -1, -10, -10};

    std::vector<std::vector<double>> ans0a = model.predict(params);
    std::vector<std::vector<double>> ans1a = model.predict_exhaust(params);
    
    // auto simres = model.simulate()

    // Casting as vectors
    std::vector< double > ans0a_vec(0u);
    for (auto & i : ans0a)
        for (auto & j: i)
            ans0a_vec.push_back(j);

    std::vector< double > ans1a_vec(0u);
    for (auto & i : ans1a)
        for (auto & j: i)
            ans1a_vec.push_back(j);

    REQUIRE_THAT(ans0a_vec, Catch::Approx(ans1a_vec).epsilon(0.001));

}
