#include "tests.hpp"
#include "../include/barry/models/geese.hpp"

BARRY_TEST_CASE("Phylo model likelihood-benchmark", "[phylo likelihood bench]")
{

    using namespace barry::counters::phylo;

    // More interesting experiment
    std::vector< std::vector<size_t> > ann = {
        {9, 9, 9, 1},
        {9, 9, 9, 1},
        {9, 0, 0, 1},
        {1, 1, 1, 1},
        {0, 9, 9, 1}
    };

    // std::vector< size_t > geneid = {0, 1, 2, 3, 4};
    // std::vector< int >  parent = {4, 3, 3, 4, -1};
    std::vector< size_t > geneid = {0, 1, 2, 3, 4};
    std::vector< int >  parent = {-1, 0, 0, 1, 1};

    std::vector< bool > duplication(geneid.size(), true);

    Geese model(ann, geneid, parent, duplication);

    // Adding terms
    counter_gains(model.get_counters(), {0, 1});
    counter_maxfuns(model.get_counters(), 0, 2);
    counter_prop_genes_changing(model.get_counters(), 1);
    counter_prop_genes_changing(model.get_counters(), 0);
    counter_k_genes_changing(model.get_counters(), 1, 0);
    counter_k_genes_changing(model.get_counters(), 1, 1);
    counter_less_than_p_prop_genes_changing(model.get_counters(), .5, 1);
    counter_less_than_p_prop_genes_changing(model.get_counters(), .5, 0);

    model.init();
    model.set_seed(100);

    for (auto S: *model.get_model()->get_stats_support())
    {
        for (auto s: S)
        {
            if (s < 0.0)
                throw std::domain_error("There shouldn't be negatives!");
        }
    }

    // Tryingout likelihood
    std::vector<double> params = {1, 1, -1, -5, 1, 1, 1, 1, -1, .4, .3,-2.5, 1.5};

    size_t n = params.size() - 1u;
    std::uniform_real_distribution<> r(-1,1);
    for (auto i = 0u; i < 1000; ++i)
    {
        size_t idx = std::floor( std::fabs( n * r( * model.get_rengine())));
        params[idx] += r(* model.get_rengine());
        (void) model.likelihood(params); 

    }

    // Checking math


}