#include "tests.hpp"
#include "../include/barry/models/geese.hpp"

BARRY_TEST_CASE("Phylo model likelihood", "[phylo likelihood]")
{

    using namespace geese;

    // More interesting experiment
    std::vector< std::vector<size_t> > ann = {
        {9, 9, 9},
        {9, 9, 9},
        {9, 0, 0},
        {1, 1, 1},
        {0, 9, 9}
    };
    // ann.resize(3u, {9, 9});

    // std::vector< size_t > geneid = {0, 1, 2, 3, 4};
    // std::vector< int >  parent = {4, 3, 3, 4, -1};
    std::vector< size_t > geneid = {0, 1, 2, 3, 4};
    std::vector< int >  parent = {-1, 0, 0, 1, 1};

    std::vector< bool > duplication(geneid.size(), true);

    Geese model(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model.get_counters(), {0, 1});
    counter_maxfuns(model.get_counters(), 2, 2);
    counter_prop_genes_changing(model.get_counters(), 1);
    counter_prop_genes_changing(model.get_counters(), 0);
    counter_k_genes_changing(model.get_counters(), 1, 0);
    counter_k_genes_changing(model.get_counters(), 1, 1);

    model.init();
    model.set_seed(100);

    // Tryingout likelihood
    std::vector<double> params = {1, 1, -1, -5, 1, 1, 1, 1, -1, -2.5};

    double ans0a = model.likelihood(params); // 0.057089548117690694
    double ans1a = model.likelihood_exhaust(params); // 0.057089548117690736
    
    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::abs(ans0a - ans1a) < .0000001);
    #else
    std::cout << "a. Likelihood: " << ans0a << std::endl;
    std::cout << "a. Likelihood(exhaust): " << ans1a << std::endl;
    #endif

    // Updating an annotation --------------------------------------------------
    ann[3u][1u] = 9u;
    ann[4u][1u] = 0u;
    Geese model2(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model2.get_counters(), {0, 1});
    counter_maxfuns(model2.get_counters(), 2, 2);
    counter_prop_genes_changing(model2.get_counters(), 1);
    counter_prop_genes_changing(model2.get_counters(), 0);
    counter_k_genes_changing(model2.get_counters(), 1, 0);
    counter_k_genes_changing(model2.get_counters(), 1, 1);

    model2.init();
    model2.set_seed(100);

    // Modifying the annotations
    model.update_annotations(3u, {1u, 9u, 1u});
    model.update_annotations(4u, {0u, 0u, 9u});

    double ans0b = model.likelihood(params);  // 0.056110910038269665
    double ans1b = model2.likelihood(params); // 0.056110910038269665

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::abs(ans0b - ans1b) < .0000001);
    #else
    std::cout << "b. Likelihood: " << ans0b << std::endl;
    std::cout << "b. Likelihood(exhaust): " << ans1b << std::endl;
    #endif

    // Going back
    model.update_annotations(3u, {1u, 1u, 1u});
    model.update_annotations(4u, {0u, 9u, 9u});

    double ans0c = model.likelihood(params); // 0.057089548117690694

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::abs(ans0a - ans0c) < .0000001);
    #endif

    std::cout << "Test for updating annotations" << std::endl;
    std::cout << "Likelihood(baseline): " << ans0a << std::endl;
    std::cout << "Likelihood(changed) : " << ans0b << std::endl;

    // Case for the jmor tree --------------------------------------------------
    // More interesting experiment
    std::vector< std::vector<size_t> > ann_jmor = {
        {0}, // a (leaf)
        {1}, // b (leaf)
        {1}, // c (leaf)
        {9}, // d (interior)
        {9}  // e (root)
    };
    // ann.resize(3u, {9, 9});

    std::vector< size_t > geneid_jmor = {0, 1, 2, 3, 4};
    std::vector< int >  parent_jmor = {3, 3, 4, 4, -1};

    std::vector< bool > duplication_jmor(geneid_jmor.size(), true);

    Geese model_jmor(ann_jmor, geneid_jmor, parent_jmor, duplication_jmor);

    counter_gains(model_jmor.get_counters(), {0});
    counter_loss(model_jmor.get_counters(), {0});

    model_jmor.init();

    std::vector<double> params_jmor = {
        -2, // Theta gains
        -4, // Theta losses
        -6  // Theta root
        };

    double par_gain = params_jmor[0];
    double par_loss = params_jmor[1];
    double par_root = 1.0/(1.0 + std::exp(-params_jmor[2]));

    // Normalizing constant (one for each parent state: 0, 1)
    double eta_0 = std::exp(2.0 * par_gain) + 2.0 * std::exp(par_gain) + 1.0;
    double eta_1 = 1.0 + 2.0 * std::exp(par_loss) + std::exp(2.0 * par_loss);

    // The manual calculation of the jmor likelihood
    double ans0_jmor = 
        par_root * (
            std::exp(par_loss) / eta_1 * std::exp(par_gain) / eta_0 +
            std::exp(par_loss) / std::pow(eta_1, 2.0)
        ) + 
        (1.0 - par_root) * (
            std::pow((std::exp(par_gain)/eta_0), 2.0) +
            std::exp(2.0 * par_gain) / eta_0 * std::exp(par_loss) / eta_1 
        );

    double ans1_jmor = model_jmor.likelihood(params_jmor);

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(std::abs(ans0_jmor - ans1_jmor) < .0000001);
    #else
    std::cout << "Etas (eta_0, eta_1)      : " << eta_0 << ", " << eta_1 << std::endl;
    std::cout << "Likelihood(jmor)         : " << ans1_jmor << std::endl;
    std::cout << "Likelihood(jmor, manual) : " << ans0_jmor << std::endl;
    #endif




}