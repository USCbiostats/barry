#include "tests.hpp"
#include "../include/barry/models/defm.hpp"
// #include "../defm.hpp"

BARRY_TEST_CASE("DEFM motif formula", "[DEFM motif formula]") {

    std::vector< size_t > res_locations1a;
    std::vector< size_t > res_locations2a;
    std::vector< size_t > res_locations3a;

    std::vector< bool > res_sign1a;
    std::vector< bool > res_sign2a;
    std::vector< bool > res_sign3a;

    std::vector< size_t > res_locations1b;
    std::vector< size_t > res_locations2b;
    std::vector< size_t > res_locations3b;
    
    std::vector< bool > res_sign1b;
    std::vector< bool > res_sign2b;
    std::vector< bool > res_sign3b;

    // Expected values
    std::vector< size_t > location_1 = {3, 7};
    std::vector< size_t > location_2 = {2, 6, 3, 7};
    std::vector< size_t > location_3 = {3};

    std::vector< bool > sign_1 = {true, true};
    std::vector< bool > sign_2 = {true, true, true, false};
    std::vector< bool > sign_3 = {true};

    defm::defm_motif_parser("{y1, y3}", res_locations1a, res_sign1a, 1, 4);
    defm::defm_motif_parser("{y1_0, y3} > {y1_1, 0y3_1}", res_locations2a, res_sign2a, 1, 4);
    defm::defm_motif_parser("{y1}", res_locations3a, res_sign3a, 1, 4);

    defm::defm_motif_parser("{y1_1, y3}", res_locations1b, res_sign1b, 1, 4);
    defm::defm_motif_parser("{y1_0, y3_0} > {y1, 0y3_1}", res_locations2b, res_sign2b, 1, 4);
    defm::defm_motif_parser("{y1_1}", res_locations3b, res_sign3b, 1, 4);


    #ifdef CATCH_CONFIG_MAIN
    REQUIRE_THAT(res_locations1a, Catch::Approx(location_1).epsilon(.001));
    REQUIRE_THAT(res_locations2a, Catch::Approx(location_2).epsilon(.001));
    REQUIRE_THAT(res_locations3a, Catch::Approx(location_3).epsilon(.001));
    REQUIRE_THAT(res_locations1a, Catch::Approx(res_locations1b).epsilon(.001));
    REQUIRE_THAT(res_locations2a, Catch::Approx(res_locations2b).epsilon(.001));
    REQUIRE_THAT(res_locations3a, Catch::Approx(res_locations3b).epsilon(.001));
    REQUIRE_THAT(res_sign1a, Catch::Approx(res_sign1b).epsilon(.001));
    REQUIRE_THAT(res_sign2a, Catch::Approx(res_sign2b).epsilon(.001));
    REQUIRE_THAT(res_sign3a, Catch::Approx(res_sign3b).epsilon(.001));
    REQUIRE_THAT(res_sign1a, Catch::Approx(sign_1).epsilon(.001));
    REQUIRE_THAT(res_sign2a, Catch::Approx(sign_2).epsilon(.001));
    REQUIRE_THAT(res_sign3a, Catch::Approx(sign_3).epsilon(.001));
    #else
    vabsdiff(res_locations1a, location_1);
    vabsdiff(res_locations2a, location_2);
    vabsdiff(res_locations3a, location_3);

    vabsdiff(res_locations1a, res_locations1b);
    vabsdiff(res_locations2a, res_locations2b);
    vabsdiff(res_locations3a, res_locations3b);

    vabsdiff(res_sign1a, res_sign1b);
    vabsdiff(res_sign2a, res_sign2b);
    vabsdiff(res_sign3a, res_sign3b);

    vabsdiff(res_sign1a, sign_1);
    vabsdiff(res_sign2a, sign_2);
    vabsdiff(res_sign3a, sign_3);
    #endif


    #ifndef CATCH_CONFIG_MAIN
    return 0;
    #endif

}