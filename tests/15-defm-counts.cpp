#include "tests.hpp"

BARRY_TEST_CASE("DEFM counts work", "[DEFM counts]") {

    using namespace barry::counters::defm;
  
    /** Array to check:
     * - Smoke
     * - Drink
     * - Marijuana
     * - Three times
     * 
     * 1 1 0
     * 0 1 1
     * 0 0 1 
     * 
     * Includes fixed effect of age (double) 31.5, 
     * and another datum that is negative
     */
    std::vector< double > covars = {
        31.5, -.26,
        31.5, -.26,
        31.5, -.26
        };
    DEFMArray A1(3, 3);
    A1.set_data(new DEFMData(&covars, 0, 3, 2), true);
    A1(0, 0) = 1;
    A1(0, 1) = 1;
    A1(1, 1) = 1;
    A1(1, 2) = 1;
    A1(2, 2) = 1;

    // Adding counters
    DEFMStatsCounter<> counter(&A1);
    counter_fixed_effect(counter.get_counters(), 0, 1.0);
    counter_fixed_effect(counter.get_counters(), 0, 0.5);
    counter_fixed_effect(counter.get_counters(), 1, 1.0);
    counter_fixed_effect(counter.get_counters(), 1, 2.0);
    counter_ones(counter.get_counters());
    counter_ones(counter.get_counters(), 0);
    counter_transition(counter.get_counters(), {0, 3, 4});
    counter_transition(counter.get_counters(), {0, 1, 4});
    counter_transition(counter.get_counters(), {0, 3, 4}, 0);
    

    std::vector< double > ans_observed = counter.count_all();
    std::vector< double > ans_expected = {
        31.5,
        std::sqrt(31.5),
        -.26,
        std::pow(-.26, 2.0),
        5.0,
        5.0 * 31.5,
        1.0,
        0.0,
        31.5
    };
    
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE_THAT(ans_expected, Catch::Approx(ans_observed).epsilon(.001));
    #else
    return 0;
    #endif


}

