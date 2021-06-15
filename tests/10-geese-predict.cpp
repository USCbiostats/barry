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

    std::vector<double> margin(ans0a_vec.size(), .0001);

    REQUIRE(vabsdiff(ans0a_vec, ans1a_vec) < margin);

    // // Updating an annotation --------------------------------------------------
    // ann[3u][1u] = 9u;
    // ann[4u][1u] = 0u;
    // Geese model2(ann, geneid, parent, duplication);

    // // Adding terms
    // counter_gains(model2.get_counters(), {0, 1});
    // counter_maxfuns(model2.get_counters(), 2, 2);

    // model2.init();
    // model2.set_seed(100);

    // // Modifying the annotations
    // model.update_annotations(3u, {1u, 9u});
    // model.update_annotations(4u, {0u, 0u});


    // double ans0b = model.likelihood(params);
    // double ans1b = model2.likelihood(params);

    // REQUIRE(std::abs(ans0b - ans1b) < .0000001);

    // // Going back
    // model.update_annotations(3u, {1u, 1u});
    // model.update_annotations(4u, {0u, 9u});

    // double ans0c = model.likelihood(params);

    // REQUIRE(std::abs(ans0a - ans0c) < .0000001);

    // std::cout << "Test for updating annotations" << std::endl;
    // std::cout << "Likelihood(baseline): " << ans0a << std::endl;
    // std::cout << "Likelihood(changed) : " << ans0b << std::endl;

}