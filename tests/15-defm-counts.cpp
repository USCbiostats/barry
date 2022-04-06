#include "tests.hpp"
#include "../include/barry/models/defm.hpp"
// #include "../defm.hpp"

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
        31.5, 31.5, 31.5,
        -.26, -.26, -.26
        };

    DEFMArray A1(3, 3);
    A1.set_data(new DEFMData(&A1, &covars[0u], 0, 2, 3), true);
    A1(0, 0) = 1;
    A1(0, 1) = 1;
    A1(1, 1) = 1;
    A1(1, 2) = 1;
    A1(2, 2) = 1;

    // Adding counters
    DEFMStatsCounter counter(&A1);
    counter_fixed_effect(counter.get_counters(), 0, 1.0);
    counter_fixed_effect(counter.get_counters(), 0, 0.5);
    counter_fixed_effect(counter.get_counters(), 1, 1.0);
    counter_fixed_effect(counter.get_counters(), 1, 2.0);
    counter_ones(counter.get_counters());
    counter_ones(counter.get_counters(), 0);
    counter_transition(counter.get_counters(), {0, 3, 4}, {});
    counter_transition(counter.get_counters(), {0, 1, 4}, {});
    counter_transition(counter.get_counters(), {0, 3, 4}, {}, 0);
    

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
    #endif

    // Creating models
    std::vector< int > id = {0,0,0,0,1,1,1,1};
    std::vector< int > Y  = {
        0, 0, 1, 1, 0, 0, 1, 1,
        0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 0, 1, 0
    };

    std::vector< double > X = {
        .1, .5, .1, 1, 1, 2, 2, 2,  // Covariate 1
        .1, 1.5, 3, 1, 1, 4, 4, 5   // Covariate 2
    };

    // Creating the model, need to pass the data
    DEFM model(&id[0u], &Y[0u], &X[0u], 8, 3, 2, 2);
    model.get_model().store_psets();
    
    // Generating the model specification
    defmcounters::counter_ones(model.get_model().get_counters());
    defmcounters::counter_ones(model.get_model().get_counters(), 0);
    defmcounters::counter_fixed_effect(model.get_model().get_counters(), 0, 1.0);

    model.init();

    std::vector< double > par0 = {.5, -.1, .1};

    #ifndef CATCH_CONFIG_MAIN

    #define GET_Y(a,b,c,d) \
        a [b * 12 + c * 6 + d]
    
    std::vector< int > out_sim(Y.size(), -1);
    model.simulate(par0, &out_sim[0u]);

    size_t ncell = 0u;
    for (size_t i = 0u; i < 2u; ++i)
    {
        for (size_t t = 0u; t < 4u; ++t)
        {
            printf("Begin % 2li: ", i);
            for (size_t j = 0u; j < 3u; ++j)
                printf("% 2i ", out_sim[ncell++]);
            printf("end\n");
        }
    }

    #undef GET_Y

    #endif

    #ifndef CATCH_CONFIG_MAIN
    auto res = model.get_model().likelihood_total(par0, true);
    model.get_model().print();
    model.get_model().print_stats(0u);
    model.get_model().print_stats(1u);
    model.get_model().print_stats(2u);
    return 0;
    #endif


}

