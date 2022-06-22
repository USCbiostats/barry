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
    counter_transition(counter.get_counters(), {0, 3, 4}, {}, 2, 3);
    counter_transition(counter.get_counters(), {0, 1, 4}, {}, 2, 3);
    counter_transition(counter.get_counters(), {0, 3, 4}, {}, 2, 3, 0);

    // With formula
    counter_transition_formula(counter.get_counters(), "{y0_0, 0y2_0} > {0y0, y2}", 2, 3);
    counter_transition(counter.get_counters(), {0, 6, 2, 8}, {true, false, false, true}, 2, 3);
    

    std::vector< double > ans_observed = counter.count_all();
    std::vector< double > ans_expected = {
        31.5,
        std::sqrt(31.5),
        -.26,
        std::pow(-.26, 2.0),
        1.0,
        31.5,
        1.0,
        0.0,
        31.5, 
        1.0, 
        1.0
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


    std::vector< int > id2(5, 1);
    std::vector< int > Y2 = {
        1, 1, 1, 0, 0,
        0, 1, 1, 0, 0,
        0, 1, 0, 0, 1,
        0, 0, 1, 0, 0,
        0, 1, 0, 0, 1,
        0, 1, 0, 1, 1,
        0, 1, 1, 0, 0,
        0, 0, 0, 0, 1,
        0, 1, 1, 1, 0,
        0, 1, 0, 0, 1
    };

    std::vector< double > X2 = {.4, .1, 1, .1, 1};
    DEFM model2(&id2[0u], &Y2[0u], &X2[0u], 5, 10, 1, 1);
    for (size_t t = 0u; t < 9; ++t)
        defmcounters::counter_transition(
            model2.get_model().get_counters(),
            {0 + 2 * t, 1 + 2 * t, 2 + 2 * t, 3 + 2 * t},
            {true, false, false, true}, 1, 10
            );

    defmcounters::counter_transition(
        model2.get_model().get_counters(),
        {18,19,0,1}, {false, true, true, false}, 1, 10
    );
    
    defmcounters::counter_ones(model2.get_model().get_counters());

    model2.init();

    std::vector< int > out_sim2(Y2.size(), -1);
    std::vector< double > params2(model2.get_model().nterms(), 0.0);
    model2.simulate(params2, &out_sim2[0u]);

    #ifndef CATCH_CONFIG_MAIN
    auto res = model.get_model().likelihood_total(par0, true);
    model.get_model().print();
    model.get_model().print_stats(0u);
    model.get_model().print_stats(1u);
    model.get_model().print_stats(2u);
    (void) model.get_model().get_stats_target();
    return 0;
    #endif


}

