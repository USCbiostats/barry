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


    // A tricky case -----------------------------------------------------------
    // This tree was causing problems in R geese.
    // PTHR11237
    // GO:0006744
    
    // Listing the data
    std::vector< std::vector< uint > > ann_R = {
        {9}, {9}, {9}, {1}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {0}, {0}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {1},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {1}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {1}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}
    };

    std::vector< uint > geneid_R = {
        59, 1, 2, 3, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
        75, 76, 77, 78, 4, 5, 6, 79, 7, 8, 80, 9, 10, 81, 82, 11, 12, 83, 13, 14,
        15, 16, 84, 17, 18, 19, 85, 86, 20, 21, 22, 23, 24, 87, 25, 26, 88, 89, 90,
        91, 27, 28, 29, 92, 93, 94, 95, 30, 96, 31, 32, 33, 34, 97, 35, 36, 37, 38,
        39, 40, 41, 98, 99, 100, 101, 102, 103, 104, 42, 43, 44, 45, 105, 106, 46,
        47, 107, 108, 48, 49, 50, 51, 52, 53, 109, 54, 55, 110, 56, 57, 58
    };

    std::vector< int > parent_R = {
        58, 59, 59, 59, 58, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
        74, 75, 76, 77, 78, 78, 77, 76, 79, 79, 75, 80, 80, 74, 81, 82, 82, 81, 83,
        83, 73, 72, 71, 84, 84, 70, 69, 85, 86, 86, 85, 68, 67, 66, 87, 87, 65, 88,
        89, 90, 91, 91, 90, 89, 92, 93, 94, 95, 95, 96, 96, 94, 93, 92, 97, 97, 97,
        88, 64, 64, 63, 62, 98, 99, 100, 101, 102, 103, 104, 104, 103, 102, 101,
        105, 106, 106, 105, 107, 108, 108, 107, 100, 99, 98, 61, 109, 109, 60, 110,
        110, -1
    };

    std::vector< bool > dpl_R = {
        false, true, true, true, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false,
        false, false, true, true, true, true, true, true, false, true, true, false,
        false, true, true, false, true, true, true, true, false, true, true, true,
        false, false, true, true, true, true, true, true, true, true, false, false,
        false, false, true, true, true, false, false, false, false, true, true,
        true, true, true, true, true, true, true, true, true, true, true, true,
        false, false, false, false, false, false, false, true, true, true, true,
        false, false, true, true, false, false, true, true, true, true, true, true,
        false, true, true, false, true, true, true
    };

    // Creating the model
    Geese model_R(ann_R, geneid_R, parent_R, dpl_R);

    /**
     # Building the model
  term_overall_changes(model2fit, duplication = TRUE)
  term_overall_changes(model2fit, duplication = FALSE)
  term_genes_changing(model2fit, duplication = TRUE)
  term_gains(model2fit, 0:(nfunctions - 1))
  term_loss(model2fit, 0:(nfunctions - 1))
  term_gains(model2fit, 0:(nfunctions - 1), FALSE)
  term_loss(model2fit, 0:(nfunctions - 1), FALSE)

  rule_limit_changes(model2fit, 0, 0, 4, TRUE)
  rule_limit_changes(model2fit, 1, 0, 4, FALSE)
     */
    counter_overall_changes(model_R.get_counters(), true);
    counter_overall_changes(model_R.get_counters(), false);
    counter_genes_changing(model_R.get_counters(), true);
    counter_gains(model_R.get_counters(), {0}, true);
    counter_loss(model_R.get_counters(), {0}, true);
    counter_gains(model_R.get_counters(), {0}, false);
    counter_loss(model_R.get_counters(), {0}, false);

    rule_dyn_limit_changes(model_R.get_support(), 0, 0, 4, true);
    rule_dyn_limit_changes(model_R.get_support(), 1, 0, 4, true);

    model_R.init();

    // These parameter estimates were obtain from MCMC in geese
    std::vector< double > params_R = {
        0, 0, -0.580299609252017, 1.40200510315799, -2.47917447288739,
        -6.5517372352121, -3.28104569986571, 0.0608426417846538
    };
    
    auto pred_R = model_R.predict(params_R, nullptr, true, false, true);


}
