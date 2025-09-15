#include "tests.hpp"
#include "../include/barry/models/defm.hpp"

BARRY_TEST_CASE("DEFM likelihood computation", "[DEFM likelihood]") {

    using namespace defm;

    // Creating test data similar to 15-defm-counts.cpp
    std::vector< int > id = {0, 0, 0, 0, 1, 1, 1, 1};
    std::vector< int > Y  = {
        0, 0, 1, 1, 0, 0, 1, 1,
        0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 1, 0, 0, 1, 0
    };

    std::vector< double > X = {
        .1, .5, .1, 1, .1, .5, 2, 2,  // Covariate 1
        .1, 1.5, 3, 1, .1, 1.5, 4, 5   // Covariate 2
    };

    // Creating the model
    DEFM model0a(&id[0u], &Y[0u], &X[0u], 8, 3, 2, 1);
    DEFM model0b(&id[0u], &Y[0u], &X[0u], 8, 3, 2, 1);
    
    // Adding counters similar to the ones in 15-defm-counts.cpp
    auto model_builder = [](auto & model_) -> void {

        model_.get_model().store_psets();
        model_.set_names({"A", "B", "C"}, {"X1", "X2"});

        counter_ones(model_.get_model().get_counters());
        counter_ones(model_.get_model().get_counters(), 0);
        counter_fixed_effect(model_.get_model().get_counters(), 0, 1.0);
        counter_fixed_effect(model_.get_model().get_counters(), 1, 1.0);
        counter_transition(model_.get_model().get_counters(), {0, 1, 2}, {}, 2, 3);
        return;
    };

    model_builder(model0a);
    model_builder(model0b);

    model0a.init(true);
    model0b.init(false);
    
    std::cout << "------------ Model 0a -----------" << std::endl;
    model0a.print();
    std::cout << "------------ Model 0b -----------" << std::endl;
    model0b.print();

    // Test parameters
    std::vector< double > par0 = {.5, -.1, .1, .2, -.05};
    std::vector< double > par1 = {-0.2, 0.3, -0.1, 0.15, 0.08};

    // Computing likelihoods using total likelihood method
    std::vector< double > logs_total0a(2), logs_total0b(2);

    // Total likelihood for all arrays
    logs_total0a[0u] = model0a.get_model().likelihood_total(par0, true);
    logs_total0a[1u] = model0a.get_model().likelihood_total(par1, true);
    logs_total0b[0u] = model0b.get_model().likelihood_total(par0, true);
    logs_total0b[1u] = model0b.get_model().likelihood_total(par1, true);

    // Printing results (should match)
    std::cout << "Total likelihoods (model 0a): ";
    for (size_t i = 0u; i < logs_total0a.size(); ++i) {
        printf("%.5f", logs_total0a[i]);
        if (i < logs_total0a.size() - 1u) printf(", ");
    }
    std::cout << std::endl;
    std::cout << "Total likelihoods (model 0b): ";
    for (size_t i = 0u; i < logs_total0b.size(); ++i) {
        printf("%.5f", logs_total0b[i]);
        if (i < logs_total0b.size() - 1u) printf(", ");
    }
    std::cout << std::endl;

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(
        model0b.get_model().support_size() <
        model0a.get_model().support_size()
    );

    // Likelihood should be finite (not NaN or infinite)
    for (size_t i = 0u; i < logs_total0a.size(); ++i) {
        REQUIRE(std::isfinite(logs_total0a[i]));
    }
    
    // Different parameters should generally give different likelihoods
    REQUIRE(logs_total0a[0u] != logs_total0a[1u]);

    REQUIRE_THAT(logs_total0a, Catch::Approx(logs_total0b).epsilon(1e-5));
    #endif

    // Test with a more complex model using formulas
    DEFM model1a(&id[0u], &Y[0u], &X[0u], 8, 3, 2, 2);
    DEFM model1b(&id[0u], &Y[0u], &X[0u], 8, 3, 2, 2);

    auto model_builder1 = [](auto & model_) -> void {
        model_.get_model().store_psets();
        model_.set_names({"A", "B", "C"}, {"X1", "X2"});

        // Using formula-based counters
        counter_ones(model_.get_model().get_counters());
        counter_transition_formula(
            model_.get_counters(), "{y0}", 2, 3, -1,
            "", &model_.get_X_names(), &model_.get_Y_names()
        );
        counter_transition_formula(
            model_.get_counters(), "{y0} x X2", 2, 3, -1,
            "", &model_.get_X_names(), &model_.get_Y_names()
        );
        counter_transition_formula(
            model_.get_counters(), "{0y0_0} > {1y0, 1y2} x X1", 2, 3, -1,
            "", &model_.get_X_names(), &model_.get_Y_names()
        );
    };

    model_builder1(model1a);
    model_builder1(model1b);
    model1a.init(true);
    model1b.init(false);

    std::cout << "------------ Model 1a -----------" << std::endl;
    model1a.print();
    std::cout << "------------ Model 1b -----------" << std::endl;
    model1b.print();

    // Checking that the support sizes in model 1b is smaller than in model 1a

    // Test parameters for formula model
    std::vector< double > par2(model1a.get_model().nterms(), 0.1);
    std::vector< double > par3(model1a.get_model().nterms(), -0.05);

    // Computing likelihoods using total likelihood method
    std::vector< double > logs_total1a(2), logs_total1b(2);

    // Total likelihood for all arrays
    logs_total1a[0u] = model1a.get_model().likelihood_total(par2, true);
    logs_total1a[1u] = model1a.get_model().likelihood_total(par3, true);
    logs_total1b[0u] = model1b.get_model().likelihood_total(par2, true);
    logs_total1b[1u] = model1b.get_model().likelihood_total(par3, true);

    // Printing results (should match)
    std::cout << "Total likelihoods (model 1a): ";
    for (size_t i = 0u; i < logs_total1a.size(); ++i) {
        printf("%.5f", logs_total1a[i]);
        if (i < logs_total1a.size() - 1u) printf(", ");
    }
    std::cout << std::endl;
    std::cout << "Total likelihoods (model 1b): ";
    for (size_t i = 0u; i < logs_total1b.size(); ++i) {
        printf("%.5f", logs_total1b[i]);
        if (i < logs_total1b.size() - 1u) printf(", ");
    }
    std::cout << std::endl;

    #ifdef CATCH_CONFIG_MAIN
    // REQUIRE(
    //     model1b.get_model().support_size() <
    //     model1a.get_model().support_size()
    // );

    // Likelihood should be finite (not NaN or infinite)
    for (size_t i = 0u; i < logs_total1a.size(); ++i) {
        REQUIRE(std::isfinite(logs_total1a[i]));
    }
    
    // Different parameters should generally give different likelihoods
    REQUIRE(logs_total1a[0u] != logs_total1a[1u]);

    REQUIRE_THAT(logs_total1a, Catch::Approx(logs_total1b).epsilon(1e-5));
    #endif

    #ifndef CATCH_CONFIG_MAIN
    return 0;
    #endif

}