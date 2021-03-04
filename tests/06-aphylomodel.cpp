#include "../include/barry/models/aphylomodel.hpp"

TEST_CASE("Phylo model likelihood", "[phylo likelihood]") {

    // More interesting experiment
    std::vector< std::vector<uint> > ann = {
      {9, 9},
      {9, 9},
      {9, 0},
      {1, 1},
      {0, 9}
    };
    // ann.resize(3u, {9, 9});

    // std::vector< uint > geneid = {0, 1, 2, 3, 4};
    // std::vector< int >  parent = {4, 3, 3, 4, -1};
    
    std::vector< uint > geneid = {0, 1, 2, 3, 4};
    std::vector< int >  parent = {-1, 0, 0, 1, 1};

    std::vector< bool > duplication(geneid.size(), true);

    APhyloModel model(ann, geneid, parent, duplication);

    // Adding terms
    phylocounters::counter_gains(&model.counters, {0, 1});
    phylocounters::counter_maxfuns(&model.counters, 2, 2);

    model.init();
    model.set_seed(100);

    // Tryingout likelihood
    std::vector<double> params = {1, 1, -1, -10, -10};


    double ans0 = model.likelihood(params);
    double ans1 = model.likelihood_exact(params);
    
    REQUIRE(std::abs(ans0 - ans1) < .0000001);
}