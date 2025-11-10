#include "tests.hpp"
#include "../include/barry/models/defm.hpp"
// #include "../defm.hpp"

BARRY_TEST_CASE("DEFM counts work", "[DEFM counts]") {

    using namespace defm;
  
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
    A1.set_data(new DEFMData(&A1, &covars[0u], 0, 2, 3, true), true);
    A1(0, 0) = 1;
    A1(0, 1) = 1;
    A1(1, 1) = 1;
    A1(1, 2) = 1;
    A1(2, 2) = 1;

    // Adding counters
    DEFMStatsCounter counter(&A1);
    counter_ones(counter.get_counters());
    counter_ones(counter.get_counters(), 0);
    counter_generic(counter.get_counters(), {0, 3, 4}, {}, 2, 3);
    counter_generic(counter.get_counters(), {0, 1, 4}, {}, 2, 3);
    counter_generic(counter.get_counters(), {0, 3, 4}, {}, 2, 3, 0);
    // counter_logit_intercept(counter.get_counters(), 3, {}, 0);

    // With formula
    counter_formula(
        counter.get_counters(), "{y0_0, 0y2_0} > {0y0, y2}", 2, 3
    );
    counter_generic(counter.get_counters(), {0, 6, 2, 8}, {true, false, false, true}, 2, 3);
    

    std::vector< double > ans_observed = counter.count_all();
    std::vector< double > ans_expected = {
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

    model.set_names({"A", "B", "C"}, {"X1", "X2"});
    
    // Generating the model specification
    counter_ones(model.get_model().get_counters());
    counter_ones(model.get_model().get_counters(), 0);

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
            printf("Begin % 2i: ", static_cast<int>(i));
            for (size_t j = 0u; j < 3u; ++j)
                printf("% 2i ", static_cast<int>(out_sim[ncell++]));
            printf("end\n");
        }
    }

    #undef GET_Y

    model.print();

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
        counter_generic(
            model2.get_model().get_counters(),
            {0 + 2 * t, 1 + 2 * t, 2 + 2 * t, 3 + 2 * t},
            {true, false, false, true}, 1, 10
            );

    counter_generic(
        model2.get_model().get_counters(),
        {18,19,0,1}, {false, true, true, false}, 1, 10
    );
    
    counter_ones(model2.get_model().get_counters());

    model2.init();

    std::vector< int > out_sim2(Y2.size(), -1);
    std::vector< double > params2(model2.get_model().nterms(), 0.0);
    model2.simulate(params2, &out_sim2[0u]);

    // Checking formulas -------------------------------------------------------

    // Creating the model, need to pass the data
    DEFM model3(&id[0u], &Y[0u], &X[0u], 8, 3, 2, 2);
    model3.get_model().store_psets();

    model3.set_names({"A", "B", "C"}, {"X1", "X2"});
    
    // Generating the model specification
    counter_ones(model3.get_model().get_counters());
    counter_formula(
        model3.get_counters(), "{y0}", 2, 3,
        &model3.get_X_names(), &model3.get_Y_names()
        );
    counter_formula(
        model3.get_counters(), "{y0} x X2", 2, 3,
        &model3.get_X_names(), &model3.get_Y_names()
        );

    counter_formula(
        model3.get_counters(), "{0y0_0} > {1y0, 1y2} x X2(Space 1)", 2, 3,
        &model3.get_X_names(), &model3.get_Y_names()
        );

    counter_formula(
        model3.get_counters(), "{0y0_0} > {1y0, 1y2} x X1(excess)", 2, 3,
        &model3.get_X_names(), &model3.get_Y_names()
        );


    model3.init();

    DEFM model3b(&id[0u], &Y[0u], &X[0u], 8, 3, 2, 2);
    model3b.get_model().store_psets();

    model3b.set_names({"A", "B", "C"}, {"X1", "X2"});
    
    // Generating the model specification
    counter_ones(model3b.get_model().get_counters());
    counter_formula(
        model3b.get_counters(), "{y0}", 2, 3,
        &model3b.get_X_names(), &model3b.get_Y_names()
        );
    counter_formula(
        model3b.get_counters(), "{y0} x X2", 2, 3,
        &model3b.get_X_names(), &model3b.get_Y_names()
        );

    counter_formula(
        model3b.get_counters(), "{0y0_0} > {1y0, 1y2} x X2", 2, 3,
        &model3b.get_X_names(), &model3b.get_Y_names()
        );

    counter_formula(
        model3b.get_counters(), "{0y0_0} > {1y0, 1y2} x X1", 2, 3,
        &model3b.get_X_names(), &model3b.get_Y_names()
        );

    model3b.init();

    #ifndef CATCH_CONFIG_MAIN
    auto res = model.get_model().likelihood_total(par0, true);
    model.get_model().print();
    model.get_model().print_stats(0u);
    model.get_model().print_stats(1u);
    model.get_model().print_stats(2u);
    (void) model.get_model().get_stats_target();
    model.print();
    model3.print();
    model3b.print();
    return 0;
    #else

    auto terms3 = model3.get_counters()->get_names();
    auto terms3b = model3b.get_counters()->get_names();

    REQUIRE_THAT(terms3, Catch::Equals(terms3b));

    #endif


}

