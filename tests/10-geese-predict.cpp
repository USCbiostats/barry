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
    
    std::vector< uint > geneid = {0, 1, 2, 3, 4};
    std::vector< int >  parent = {-1, 0, 0, 1, 1};
    /**
     * 0__1__3
     * |  |
     * |  |__4
     * |
     * |__2
     */

    std::vector< bool > duplication(geneid.size(), true);

    Geese model(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model.get_counters(), {0, 1});
    counter_overall_changes(model.get_counters());

    model.init();
    model.set_seed(100);

    // Model parameters to test
    std::vector<double> params = {1, -1, -.5, -5, -5};

    std::vector<std::vector<double>> ans0a = model.predict_exhaust(params);
    std::vector<std::vector<double>> ans1a = model.predict(params, nullptr, true);
    std::vector<std::vector<double>> ans2a = model.predict_sim(params, false, 100000);
    
    printf_barry("predict_exhaust():\n");
    for (auto & a : ans0a) {
        for (auto & v : a) {
            printf_barry("%.6f, ", v);
        }
        printf_barry("\n");
    }
    printf_barry("Versus predict():\n");
    for (auto & a : ans1a) {
        for (auto & v : a) {
            printf_barry("%.6f, ", v);
        }
        printf_barry("\n");
    }

    printf_barry("Versus predict_sim():\n");
    for (auto & a : ans2a) {
        for (auto & v : a) {
            printf_barry("%.6f, ", v);
        }
        printf_barry("\n");
    }

    // Casting as vectors
    std::vector< double > ans0a_vec(0u);
    for (auto & i : ans0a)
        for (auto & j: i)
            ans0a_vec.push_back(j);

    std::vector< double > ans1a_vec(0u);
    for (auto & i : ans1a)
        for (auto & j: i)
            ans1a_vec.push_back(j);

    std::vector< double > ans2a_vec(0u);
    for (auto & i : ans2a)
        for (auto & j: i)
            ans2a_vec.push_back(j);

    REQUIRE_THAT(ans0a_vec, Catch::Approx(ans2a_vec).margin(0.025));
    REQUIRE_THAT(ans0a_vec, Catch::Approx(ans1a_vec).margin(0.025));

}
