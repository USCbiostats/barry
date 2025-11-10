#include "tests.hpp"
#include "../include/barry/models/defm.hpp"

BARRY_TEST_CASE("DEFM motif formula with order 2", "[DEFM motif formula order 2]") {

    std::vector< size_t > res_locations1;
    std::vector< size_t > res_locations2;
    std::vector< size_t > res_locations3;
    std::vector< size_t > res_locations4;

    std::vector< bool > res_sign1;
    std::vector< bool > res_sign2;
    std::vector< bool > res_sign3;
    std::vector< bool > res_sign4;

    std::string covar_name1, covar_name2, covar_name3, covar_name4;

    // Test 1: Simple order 2 transition "{0y1} > {y1} > {0y1}"
    // Variable 1 transitions from 0 to 1, and back to 0
    // With m_order=2, y_ncol=4, we have:
    // - y1_0 at position 1*(2+1)+0 = 3
    // - y1_1 at position 1*(2+1)+1 = 4  
    // - y1_2 at position 1*(2+1)+2 = 5
    defm::defm_motif_parser(
        "{0y1} > {y1} > {0y1}", res_locations1, res_sign1, 2, 4, covar_name1
    );

    std::vector< size_t > expected_loc1 = {3, 4, 5};
    std::vector< bool > expected_sign1 = {false, true, false};

    // Test 2: Order 2 transition with explicit indices
    defm::defm_motif_parser(
        "{0y1_0} > {y1_1} > {0y1_2}", res_locations2, res_sign2, 2, 4, covar_name2
    );

    std::vector< size_t > expected_loc2 = {3, 4, 5};
    std::vector< bool > expected_sign2 = {false, true, false};

    // Test 3: Multiple variables in order 2 transition
    // "{y0_0, y2_0} > {y0_1, 0y2_1} > {0y0_2, y2_2}"
    // y0_0: 0*(2+1)+0 = 0, y2_0: 2*(2+1)+0 = 6
    // y0_1: 0*(2+1)+1 = 1, 0y2_1: 2*(2+1)+1 = 7
    // 0y0_2: 0*(2+1)+2 = 2, y2_2: 2*(2+1)+2 = 8
    defm::defm_motif_parser(
        "{y0_0, y2_0} > {y0_1, 0y2_1} > {0y0_2, y2_2}",
        res_locations3, res_sign3, 2, 4, covar_name3
    );

    std::vector< size_t > expected_loc3 = {0, 6, 1, 7, 2, 8};
    std::vector< bool > expected_sign3 = {true, true, true, false, false, true};

    // Test 4: Order 2 with covariate
    defm::defm_motif_parser(
        "{0y1_0} > {y1_1} > {0y1_2} x cov1",
        res_locations4, res_sign4, 2, 4, covar_name4
    );

    std::vector< size_t > expected_loc4 = {3, 4, 5};
    std::vector< bool > expected_sign4 = {false, true, false};

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(res_locations1 == expected_loc1);
    REQUIRE(res_sign1 == expected_sign1);
    REQUIRE(res_locations2 == expected_loc2);
    REQUIRE(res_sign2 == expected_sign2);
    REQUIRE(res_locations3 == expected_loc3);
    REQUIRE(res_sign3 == expected_sign3);
    REQUIRE(res_locations4 == expected_loc4);
    REQUIRE(res_sign4 == expected_sign4);
    REQUIRE(covar_name4 == "cov1");
    #else
    vabsdiff(res_locations1, expected_loc1);
    vabsdiff(res_sign1, expected_sign1);
    vabsdiff(res_locations2, expected_loc2);
    vabsdiff(res_sign2, expected_sign2);
    vabsdiff(res_locations3, expected_loc3);
    vabsdiff(res_sign3, expected_sign3);
    vabsdiff(res_locations4, expected_loc4);
    vabsdiff(res_sign4, expected_sign4);
    #endif

    // Now test with counter_generic to ensure it works correctly
    // Create a simple DEFM model with m_order=2
    std::vector< int > id = {0, 0, 0, 0};
    std::vector< int > Y  = {
        0, 0, 1, 1,  // Time 0
        0, 0, 1, 0,  // Time 1
        0, 0, 0, 1   // Time 2
    };
    std::vector< double > X = {
        1.0, 1.0, 1.0, 1.0  // Dummy covariate
    };

    defm::DEFM model(&id[0u], &Y[0u], &X[0u], 4, 4, 1, 2);
    model.get_model().store_psets();

    // Add a counter using the parsed formula data
    defm::counter_generic(
        model.get_model().get_counters(),
        res_locations1,
        res_sign1,
        2,
        4
    );

    model.init();

    // Just verify that the model initializes without errors
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(model.get_model().nterms() == 1);
    #else
    if (model.get_model().nterms() != 1) {
        printf("Expected 1 term, got %zu\n", model.get_model().nterms());
        return 1;
    }
    #endif

    #ifndef CATCH_CONFIG_MAIN
    return 0;
    #endif

}
