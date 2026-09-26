#include "tests.hpp"
#include "../include/barry/models/defm.hpp"

// Regression tests for defm_motif_parser:
//
//  - Motif variables were searched in the whole formula, including the
//    covariate, so a covariate like 'Day1' silently added y1 to the motif
//    (or threw a misleading error).
//  - The syntax check accepted '}x Covar' but the covariate extraction
//    required whitespace before 'x', silently dropping the interaction.
//  - The LHS of a two-group transition accepted rows equal to m_order, and
//    explicit rows were not range-checked when m_order == 1.

BARRY_TEST_CASE("DEFM motif formula parser edge cases", "[DEFM formula parser]") {

    std::vector< size_t > loc;
    std::vector< bool > sgn;
    std::string covar;

    // Order 1, 4 columns: location = col * 2 + row

    // Covariate names containing y<digit> are not motif terms
    defm::defm_motif_parser("{y0} > {y0} x Day1", loc, sgn, 1, 4, covar);
    REQUIRE(loc == std::vector< size_t >({0u, 1u}));
    REQUIRE(sgn == std::vector< bool >({true, true}));
    REQUIRE(covar == "Day1");

    defm::defm_motif_parser("{y0} > {y1} x Day1", loc, sgn, 1, 4, covar);
    REQUIRE(loc == std::vector< size_t >({0u, 3u}));
    REQUIRE(covar == "Day1");

    defm::defm_motif_parser("{y0} x Day5", loc, sgn, 1, 4, covar);
    REQUIRE(loc == std::vector< size_t >({1u}));
    REQUIRE(covar == "Day5");

    // Same for the label
    defm::defm_motif_parser("{y0, 0y2} x Day1 (y3 label)", loc, sgn, 1, 4, covar);
    REQUIRE(loc == std::vector< size_t >({1u, 5u}));
    REQUIRE(sgn == std::vector< bool >({true, false}));
    REQUIRE(covar == "Day1");

    // Order 2, multi-group mode (location = col * 3 + row)
    defm::defm_motif_parser("{y0} > {y1} > {y2} x Hwy3", loc, sgn, 2, 4, covar);
    REQUIRE(loc == std::vector< size_t >({0u, 4u, 8u}));
    REQUIRE(covar == "Hwy3");

    // No whitespace before 'x' keeps the covariate
    defm::defm_motif_parser("{y0} > {y1}x Female", loc, sgn, 1, 4, covar);
    REQUIRE(loc == std::vector< size_t >({0u, 3u}));
    REQUIRE(covar == "Female");

    defm::defm_motif_parser("{y0}x Female", loc, sgn, 1, 4, covar);
    REQUIRE(loc == std::vector< size_t >({1u}));
    REQUIRE(covar == "Female");

    // Trailing whitespace is not part of the covariate name
    defm::defm_motif_parser("{y0} x Female  ", loc, sgn, 1, 4, covar);
    REQUIRE(covar == "Female");

    // No covariate resets the name
    defm::defm_motif_parser("{y0} > {y1}", loc, sgn, 1, 4, covar);
    REQUIRE(covar == "");

    // 'x' glued to the covariate name is ambiguous: syntax error
    REQUIRE_THROWS_AS(
        defm::defm_motif_parser("{y0} xFemale", loc, sgn, 1, 4, covar),
        std::logic_error
    );

    // LHS cannot include the current time (order 1)
    REQUIRE_THROWS_AS(
        defm::defm_motif_parser("{y0_1} > {y1}", loc, sgn, 1, 4, covar),
        std::logic_error
    );

    // LHS cannot include the current time (order 2, two-group mode)
    REQUIRE_THROWS_AS(
        defm::defm_motif_parser("{y0_2} > {y1}", loc, sgn, 2, 4, covar),
        std::logic_error
    );

    // Out-of-range rows with order 1 throw (instead of indexing out of bounds)
    REQUIRE_THROWS_AS(
        defm::defm_motif_parser("{y0_5} > {y1}", loc, sgn, 1, 4, covar),
        std::logic_error
    );

    REQUIRE_THROWS_AS(
        defm::defm_motif_parser("{y0} > {y1_5}", loc, sgn, 1, 4, covar),
        std::logic_error
    );

}
