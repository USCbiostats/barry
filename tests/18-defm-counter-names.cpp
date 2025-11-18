#include "tests.hpp"
#include "../include/barry/models/defm.hpp"

BARRY_TEST_CASE("DEFM counter names", "[DEFM counter names]") {

    using namespace defm;
    
    // Setup data for testing - markov order 2
    // We need id_length = number of time observations
    std::vector<int> id = {0, 0, 0, 0};
    std::vector<int> Y = {
        1, 0, 1, 0,  // Time 0
        0, 1, 0, 1,  // Time 1
        1, 0, 1, 0   // Time 2
    };
    std::vector<double> X = {1.0, 1.0, 1.0, 1.0};

    // =========================================================================
    // Test 1: counter_generic with order 2, single variable transition
    // Expected name should show transition states
    // =========================================================================
    DEFM model1(&id[0u], &Y[0u], &X[0u], 4, 4, 1, 2);
    model1.get_model().store_psets();
    model1.set_names({"A", "B", "C", "D"}, {"X1"});
    
    // Motif: A goes from 0 to 1 to 0  (positions 0, 1, 2)
    // Position calculation: y_col * (m_order + 1) + y_row
    // A_0 = 0 * 3 + 0 = 0
    // A_1 = 0 * 3 + 1 = 1
    // A_2 = 0 * 3 + 2 = 2
    counter_generic(
        model1.get_model().get_counters(),
        {0, 1, 2},
        {false, true, false},
        2, 4,
        -1,
        &model1.get_X_names(),
        &model1.get_Y_names()
    );
    
    model1.init();
    auto names1 = model1.get_counters()->get_names();
    
    #ifndef CATCH_CONFIG_MAIN
    printf("Test 1 (counter_generic, order 2): '%s'\n", names1[0].c_str());
    #endif
    
    // For order 2, the name should show the motif structure
    // The name groups prior states (t=0, t=1) together and shows final state
    // Expected format: "Motif {A-, A+}⇨{A-}" or similar
    
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(names1.size() == 1);
    REQUIRE(names1[0].find("A") != std::string::npos);
    REQUIRE(names1[0].find("Motif") != std::string::npos);
    // Check that the name contains indicators for the transitions
    REQUIRE((names1[0].find("-") != std::string::npos || 
            names1[0].find("0") != std::string::npos));
    REQUIRE((names1[0].find("+") != std::string::npos ||
            names1[0].find("1") != std::string::npos));
    // Should contain the arrow symbol for transition
    REQUIRE((names1[0].find("->") != std::string::npos ||
            names1[0].find(u8"\u21E8") != std::string::npos ||
            names1[0].find(">") != std::string::npos));
    #endif

    // =========================================================================
    // Test 2: counter_formula with 2-group format (backwards compatible)
    // For order 2, this uses format: {vars_at_0_and_1} > {vars_at_2}
    // Note: Formula parser uses y0, y1, etc. notation, not custom names
    // =========================================================================
    DEFM model2(&id[0u], &Y[0u], &X[0u], 4, 4, 1, 2);
    model2.get_model().store_psets();
    model2.set_names({"A", "B", "C", "D"}, {"X1"});
    
    // Use 2-group format: LHS has times 0-1 (explicit), RHS has time 2 (implicit)
    // Note: must use y0, y1, etc. in formula, not custom names
    counter_formula(
        model2.get_counters(),
        "{0y0_0, y0_1} > {0y0}",
        2, 4,
        &model2.get_X_names(),
        &model2.get_Y_names()
    );
    
    model2.init();
    auto names2 = model2.get_counters()->get_names();
    
    #ifndef CATCH_CONFIG_MAIN
    printf("Test 2 (counter_formula, 2-group): '%s'\n", names2[0].c_str());
    #endif
    
    #ifdef CATCH_CONFIG_MAIN
    // Names from counter_generic and counter_formula should match
    REQUIRE(names1[0] == names2[0]);
    #endif

    // =========================================================================
    // Test 3: counter_formula with 3-group format (explicit time points)
    // For order 2, this uses format: {vars_at_0} > {vars_at_1} > {vars_at_2}
    // =========================================================================
    DEFM model3(&id[0u], &Y[0u], &X[0u], 4, 4, 1, 2);
    model3.get_model().store_psets();
    model3.set_names({"A", "B", "C", "D"}, {"X1"});
    
    // Use 3-group format: explicit time for each group
    counter_formula(
        model3.get_counters(),
        "{0y0} > {y0} > {0y0}",
        2, 4,
        &model3.get_X_names(),
        &model3.get_Y_names()
    );
    
    model3.init();
    auto names3 = model3.get_counters()->get_names();
    
    #ifndef CATCH_CONFIG_MAIN
    printf("Test 3 (counter_formula, 3-group): '%s'\n", names3[0].c_str());
    #endif
    
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(names3.size() == 1);
    // Should still produce the same name as counter_generic
    // because both represent the same motif
    REQUIRE(names3[0] == names1[0]);
    #endif

    // =========================================================================
    // Test 4: More complex - multiple variables with 3-group format
    // =========================================================================
    DEFM model4(&id[0u], &Y[0u], &X[0u], 4, 4, 1, 2);
    model4.get_model().store_psets();
    model4.set_names({"A", "B", "C", "D"}, {"X1"});
    
    // Multiple variables in 3-group format
    counter_formula(
        model4.get_counters(),
        "{y0, y2} > {y0, 0y2} > {0y0, y2}",
        2, 4,
        &model4.get_X_names(),
        &model4.get_Y_names()
    );
    
    model4.init();
    auto names4 = model4.get_counters()->get_names();
    
    #ifndef CATCH_CONFIG_MAIN
    printf("Test 4 (counter_formula, 3-group, multi var): '%s'\n", names4[0].c_str());
    #endif
    
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(names4.size() == 1);
    // Should contain both A and C (custom names for y0 and y2)
    REQUIRE(names4[0].find("A") != std::string::npos);
    REQUIRE(names4[0].find("C") != std::string::npos);
    REQUIRE(names4[0].find("Motif") != std::string::npos);
    #endif

    // =========================================================================
    // Test 5: Verify counter_generic produces same name for complex motif
    // =========================================================================
    DEFM model5(&id[0u], &Y[0u], &X[0u], 4, 4, 1, 2);
    model5.get_model().store_psets();
    model5.set_names({"A", "B", "C", "D"}, {"X1"});
    
    // Same motif as model4: {y0, y2} > {y0, 0y2} > {0y0, y2}
    // Positions: A_0=0, C_0=6, A_1=1, C_1=7, A_2=2, C_2=8
    // Coords for y0: col=0, y2: col=2
    // y0_0: 0*3+0=0, y2_0: 2*3+0=6, y0_1: 0*3+1=1, y2_1: 2*3+1=7, y0_2: 0*3+2=2, y2_2: 2*3+2=8
    counter_generic(
        model5.get_model().get_counters(),
        {0, 6, 1, 7, 2, 8},
        {true, true, true, false, false, true},
        2, 4,
        -1,
        &model5.get_X_names(),
        &model5.get_Y_names()
    );
    
    model5.init();
    auto names5 = model5.get_counters()->get_names();
    
    #ifndef CATCH_CONFIG_MAIN
    printf("Test 5 (counter_generic, multi var): '%s'\n", names5[0].c_str());
    #endif
    
    #ifdef CATCH_CONFIG_MAIN
    // Names from counter_generic and counter_formula should match
    REQUIRE(names4[0] == names5[0]);
    #endif

    // =========================================================================
    // Test 6: Order 1 for comparison - ensure it still works correctly
    // =========================================================================
    std::vector<int> Y_o1 = {
        1, 0, 1,  // Time 0
        0, 1, 0   // Time 1
    };
    std::vector<int> id_o1 = {0, 0};
    std::vector<double> X_o1 = {1.0, 1.0};
    
    DEFM model6(&id_o1[0u], &Y_o1[0u], &X_o1[0u], 2, 3, 1, 1);
    model6.get_model().store_psets();
    model6.set_names({"A", "B", "C"}, {"X1"});
    
    counter_formula(
        model6.get_counters(),
        "{0y0} > {y0}",
        1, 3,
        &model6.get_X_names(),
        &model6.get_Y_names()
    );
    
    model6.init();
    auto names6 = model6.get_counters()->get_names();
    
    #ifndef CATCH_CONFIG_MAIN
    printf("Test 6 (counter_formula, order 1): '%s'\n", names6[0].c_str());
    #endif
    
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(names6.size() == 1);
    REQUIRE(names6[0].find("A") != std::string::npos);
    #endif

    // =========================================================================
    // Test 7: Test with numeric names (no custom names set)
    // =========================================================================
    DEFM model7(&id[0u], &Y[0u], &X[0u], 4, 4, 1, 2);
    model7.get_model().store_psets();
    // Don't set names - will use default y0, y1, etc.
    
    counter_formula(
        model7.get_counters(),
        "{0y0_0, y0_1} > {0y0}",
        2, 4,
        nullptr,
        nullptr
    );
    
    model7.init();
    auto names7 = model7.get_counters()->get_names();
    
    #ifndef CATCH_CONFIG_MAIN
    printf("Test 7 (counter_formula, default names): '%s'\n", names7[0].c_str());
    #endif
    
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(names7.size() == 1);
    REQUIRE(names7[0].find("y0") != std::string::npos);
    #endif

    // =========================================================================
    // Test 8: Test with covariate - name should include covariate
    // =========================================================================
    DEFM model8(&id[0u], &Y[0u], &X[0u], 4, 4, 1, 2);
    model8.get_model().store_psets();
    model8.set_names({"A", "B", "C", "D"}, {"X1"});
    
    counter_formula(
        model8.get_counters(),
        "{0y0_0, y0_1} > {0y0} x X1",
        2, 4,
        &model8.get_X_names(),
        &model8.get_Y_names()
    );
    
    model8.init();
    auto names8 = model8.get_counters()->get_names();
    
    #ifndef CATCH_CONFIG_MAIN
    printf("Test 8 (counter_formula, with covariate): '%s'\n", names8[0].c_str());
    #endif
    
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(names8.size() == 1);
    REQUIRE(names8[0].find("A") != std::string::npos);
    REQUIRE(names8[0].find("X1") != std::string::npos);
    #endif

    #ifndef CATCH_CONFIG_MAIN
    return 0;
    #endif

}
