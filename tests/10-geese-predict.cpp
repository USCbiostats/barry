#include "tests.hpp"

#include "../include/barry/models/geese.hpp"

BARRY_TEST_CASE("Geese model prediction", "[geese prediction]") {

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
     * 0__1 (9, 9)__3 (1, 1)
     * |  |
     * |  |__4 (0, 9)
     * |
     * |__2 (9, 0)
     */

    std::vector< bool > duplication = {true, false, true, false, true};

    Geese model(ann, geneid, parent, duplication);
    Flock model_flock;

    (void) model_flock.add_data(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model.get_counters(), {0, 1});
    counter_gains(model.get_counters(), {0, 1}, false); // speciation
    counter_overall_changes(model.get_counters());

    // Adding terms
    counter_gains(model_flock.get_counters(), {0, 1});
    counter_gains(model_flock.get_counters(), {0, 1}, false); // speciation
    counter_overall_changes(model_flock.get_counters());

    model.init();
    model_flock.init();

    // model.get_model()->print_stats(3);
    model.set_seed(100);

    // Model parameters to test
    std::vector<double> params = {1, -1, -1.5, -1.5, -.5, -5, -5};

    std::vector<std::vector<double>> ans0a = model.predict_exhaust(params);
    std::vector<std::vector<double>> ans1a = model.predict(params, nullptr, true);

    #ifdef BARRY_VALGRIND
    std::vector<std::vector<double>> ans2a = model.predict_sim(params, false, 1);
    #else
    std::vector<std::vector<double>> ans2a = model.predict_sim(params, false, 50000);
    #endif
    
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

    #ifdef CATCH_CONFIG_MAIN
        #ifndef BARRY_VALGRIND
            REQUIRE_THAT(ans0a_vec, Catch::Approx(ans2a_vec).margin(0.025));
        #endif
            REQUIRE_THAT(ans0a_vec, Catch::Approx(ans1a_vec).margin(0.025));
    #endif

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

    counter_overall_changes(model_R.get_counters(), 1u);
    counter_overall_changes(model_R.get_counters(), 0u);
    counter_genes_changing(model_R.get_counters(), 1u);
    counter_gains(model_R.get_counters(), {0}, 1u);
    counter_loss(model_R.get_counters(), {0}, 1u);
    counter_gains(model_R.get_counters(), {0}, 0u);
    counter_loss(model_R.get_counters(), {0}, 0u);

    rule_dyn_limit_changes(model_R.get_support(), 0, 0, 4, 1u);
    rule_dyn_limit_changes(model_R.get_support(), 0, 0, 4, 0u);

    model_R.init();

    // These parameter estimates were obtain from MCMC in geese
    std::vector< double > params_R = {
        0, 0, -0.580299609252017, 1.40200510315799, -2.47917447288739,
        -6.5517372352121, -3.28104569986571, 0.0608426417846538
    };
    
    model_R.set_seed(8882);
    auto pred_R = model_R.predict(params_R, nullptr, true, false, true);

    #ifdef CATCH_CONFIG_MAIN
        #ifdef BARRY_VALGRIND
            auto pred_sim_R = model_R.predict_sim(params_R, false, 1);
        #else
            auto pred_sim_R = model_R.predict_sim(params_R, false, 50000);
        #endif
    #else
        auto pred_sim_R = model_R.predict_sim(params_R, false, 1);
    #endif

    std::vector< double > ansR_vec(0u);
    for (auto & i : pred_R)
        for (auto & j: i)
            ansR_vec.push_back(j);

    std::vector< double > ansR_sim_vec(0u);
    for (auto & i : pred_sim_R)
        for (auto & j: i)
            ansR_sim_vec.push_back(j);

    printf_barry("Versus (real dat) :\n");
    for (int i = 0; i < static_cast<int>(ansR_sim_vec.size()); ++i)
    {

        printf_barry(
            "% 3i: % 7.4f - % 7.4f = % 7.4f\n",
            i, ansR_vec[i], ansR_sim_vec[i], ansR_vec[i] - ansR_sim_vec[i]
            );

    }

    #ifdef CATCH_CONFIG_MAIN
        #ifndef BARRY_VALGRIND
            REQUIRE_THAT(ansR_vec, Catch::Approx(ansR_sim_vec).margin(0.025));    
        #endif
    #endif

    // Another tricky case -----------------------------------------------------
    std::vector< std::vector< uint > > ann_R2 = {
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {0}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {1}, {9}, {9}, {1}, {9}, {9}, {9}, {9}, {9}, {9}, {1},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {0}, {9},
        {9}, {9}, {9}, {0}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {1}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {0}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9}, {9},
        {9}, {9}, {9}, {9}, {9}, {9}
    };

    std::vector< uint > geneid_R2 = {
        145, 146, 147, 148, 149, 150, 151, 152, 153, 1, 2, 154, 155, 156, 157,
        158, 159, 160, 161, 3, 4, 5, 6, 162, 7, 8, 163, 164, 165, 166, 9, 10,
        11, 167, 12, 13, 14, 168, 15, 16, 17, 18, 169, 19, 20, 21, 170, 171,
        172, 22, 23, 24, 173, 174, 25, 26, 27, 175, 176, 177, 178, 179, 180,
        181, 182, 183, 184, 28, 29, 30, 31, 185, 186, 32, 33, 187, 34, 35, 188,
        189, 190, 191, 36, 37, 38, 192, 39, 40, 193, 41, 42, 43, 44, 45, 46,
        194, 195, 47, 48, 196, 49, 50, 197, 198, 199, 200, 201, 202, 203, 204,
        51, 52, 53, 54, 205, 55, 56, 206, 207, 208, 57, 58, 59, 209, 60, 61,
        62, 210, 63, 64, 211, 212, 65, 66, 67, 213, 214, 215, 216, 217, 218,
        219, 220, 221, 68, 69, 70, 71, 222, 72, 73, 223, 224, 225, 74, 75, 76,
        226, 77, 78, 79, 80, 81, 227, 228, 82, 83, 84, 229, 230, 231, 232, 233,
        234, 235, 236, 237, 238, 85, 86, 87, 88, 239, 89, 90, 240, 241, 242, 91,
        92, 93, 243, 94, 95, 96, 97, 244, 98, 99, 100, 245, 246, 247, 101, 102,
        248, 103, 104, 105, 249, 250, 251, 252, 253, 254, 255, 106, 256, 107,
        108, 257, 258, 259, 109, 110, 111, 260, 112, 113, 114, 115, 261, 116,
        262, 117, 118, 119, 263, 264, 120, 121, 122, 265, 123, 124, 125, 126,
        266, 127, 128, 129, 130, 131, 267, 268, 269, 132, 133, 270, 134, 135,
        136, 137, 271, 138, 139, 272, 140, 141, 142, 143, 144
    };

    std::vector< int > parent_R2 = {
        144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 153, 152, 154, 155,
        156, 157, 158, 159, 160, 161, 161, 160, 159, 158, 162, 162, 157, 163,
        164, 165, 166, 166, 166, 165, 167, 167, 164, 163, 168, 168, 156, 155,
        154, 169, 169, 151, 150, 170, 171, 172, 172, 171, 170, 173, 174, 174,
        173, 149, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 184, 183,
        182, 181, 185, 186, 186, 185, 187, 187, 180, 188, 189, 190, 191, 191,
        190, 189, 192, 192, 188, 193, 193, 179, 178, 177, 176, 175, 194, 195,
        195, 194, 196, 196, 149, 197, 198, 199, 200, 201, 202, 203, 204, 204,
        203, 202, 201, 205, 205, 200, 206, 207, 208, 208, 207, 206, 209, 209,
        199, 198, 210, 210, 197, 211, 212, 212, 211, 149, 213, 214, 215, 216,
        217, 218, 219, 220, 221, 221, 220, 219, 218, 222, 222, 217, 223, 224,
        225, 225, 224, 223, 226, 226, 216, 215, 214, 213, 227, 228, 228, 227,
        149, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 238, 237, 236,
        235, 239, 239, 234, 240, 241, 242, 242, 241, 240, 243, 243, 233, 232,
        231, 244, 244, 230, 229, 245, 246, 247, 247, 246, 248, 248, 245, 149,
        249, 250, 251, 252, 253, 254, 255, 255, 256, 256, 254, 257, 258, 259,
        259, 258, 257, 260, 260, 253, 252, 251, 261, 261, 262, 262, 250, 249,
        263, 264, 264, 263, 148, 265, 265, 265, 147, 146, 266, 266, 266, 266,
        266, 145, 267, 268, 269, 269, 268, 270, 270, 270, 270, 267, 271, 271,
        144, 272, 272, 272, 272, -1
    };

    std::vector< bool > dpl_R2 = {
        false, false, false, false, true, false, false, true, false, true, true,
        false, false, false, false, false, false, false, false, true, true,
        true, true, false, true, true, false, false, false, true, true, true,
        true, true, true, true, true, false, true, true, true, true, false,
        true, true, true, true, false, false, true, true, true, false, false,
        true, true, true, false, false, false, false, false, false, false,
        false, false, false, true, true, true, true, true, false, true, true,
        false, true, true, false, false, false, true, true, true, true, true,
        true, true, false, true, true, true, true, true, true, false, false,
        true, true, true, true, true, false, false, false, false, false, false,
        false, false, true, true, true, true, false, true, true, false, false,
        true, true, true, true, false, true, true, true, false, true, true,
        false, false, true, true, true, false, false, false, false, false,
        false, false, false, false, true, true, true, true, false, true, true,
        false, false, false, true, true, true, false, true, true, true, true,
        true, false, false, true, true, true, false, false, false, false, false,
        false, false, false, false, false, true, true, true, true, false, true,
        true, false, false, false, true, true, true, false, true, true, true,
        true, false, true, true, true, false, false, true, true, true, true,
        true, true, true, false, false, false, false, false, false, false, true,
        false, true, true, false, false, false, true, true, true, false, true,
        true, true, true, false, true, true, true, true, true, false, false,
        true, true, true, true, true, true, true, true, true, true, true, true,
        true, true, false, false, false, true, true, true, true, true, true,
        true, false, true, true, true, true, true, true, true, false
    };


    // Creating the model
    Geese model_R2(ann_R2, geneid_R2, parent_R2, dpl_R2);

    counter_overall_changes(model_R2.get_counters(), 1u);
    counter_overall_changes(model_R2.get_counters(), 0u);
    counter_genes_changing(model_R2.get_counters(), 1u);
    counter_gains(model_R2.get_counters(), {0}, 1u);
    counter_loss(model_R2.get_counters(), {0}, 1u);
    counter_gains(model_R2.get_counters(), {0}, 0u);
    counter_loss(model_R2.get_counters(), {0}, 0u);

    rule_dyn_limit_changes(model_R2.get_support(), 0, 0, 4, 1u);
    rule_dyn_limit_changes(model_R2.get_support(), 1, 0, 4, 0u);

    model_R2.init();

    // These parameter estimates were obtain from MCMC in geese
    params_R = {
        0, 0, 2.05977681012341, 2.66090067970927, -0.0123964396552612,
        -7.07827462428353, -3.26318748366662, -0.0966378176432868
    };
    
    model_R2.set_seed(8882);
    auto pred_R2 = model_R2.predict(params_R, nullptr, true, false, true);

    #ifdef CATCH_CONFIG_MAIN
        #ifdef BARRY_VALGRIND
            auto pred_sim_R2 = model_R2.predict_sim(params_R, false, 1);
        #else
            auto pred_sim_R2 = model_R2.predict_sim(params_R, false, 50000);
        #endif
    #else
        auto pred_sim_R2 = model_R2.predict_sim(params_R, false, 1);
    #endif

    std::vector< double > ansR2_vec(0u);
    for (auto & i : pred_R2)
        for (auto & j: i)
            ansR2_vec.push_back(j);

    std::vector< double > ansR2_sim_vec(0u);
    for (auto & i : pred_sim_R2)
        for (auto & j: i)
            ansR2_sim_vec.push_back(j);

    #ifdef CATCH_CONFIG_MAIN
        #ifndef BARRY_VALGRIND
            REQUIRE_THAT(ansR2_vec, Catch::Approx(ansR2_sim_vec).margin(0.025));    
        #endif
    #endif

}
