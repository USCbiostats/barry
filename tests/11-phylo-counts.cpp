// #include <vector>
// #include "../include/barry.hpp"
// #include "catch.hpp"

TEST_CASE("Phylo counts work", "[phylo counts]") {

    using namespace barry::counters::phylo;
  
    /** PhyloArrays to check:
     * 
     * // Subfun for both Duplication and speciation
     * 1 | 0 1
     * 1 | 1 0
     * 0 | 1 1
     * 
     */
    
    PhyloArray A1_dpl(3, 2);
    A1_dpl.set_data(new NodeData({1.0, 1.5}, {true, true}, true), true);
    A1_dpl(1, 0) = 1;
    A1_dpl(0, 1) = 1;
    A1_dpl(2, 0) = 1;
    A1_dpl(2, 1) = 1;

    PhyloArray A1_spe(A1_dpl, true);
    A1_spe.D()->duplication = false;
    
    // Adding counters
    PhyloStatsCounter counter(&A1_dpl);
    counter_gains(counter.get_counters(), {0u, 1u, 2u}, true);
    counter_gains(counter.get_counters(), {0u, 1u, 2u}, false);
    counter_loss(counter.get_counters(), {0u, 1u, 2u}, true);
    counter_loss(counter.get_counters(), {0u, 1u, 2u}, false);
    counter_overall_changes(counter.get_counters(), true);
    counter_overall_changes(counter.get_counters(), false);

    std::vector< double > ans1_dpl_obs      = counter.count_all();
    std::vector< double > ans1_dpl_expected = {
        0, 0, 2,//counter_gains (dpl)
        0, 0, 0,//counter_gains (spe)
        1, 1, 0,//counter_loss (dpl)
        0, 0, 0,//counter_loss (spe)
        4, //counter_overall_changes (dpl)
        0 // counter_overall_changes (spe)
    };

    REQUIRE_THAT(ans1_dpl_expected, Catch::Approx(ans1_dpl_obs).epsilon(.001));

    counter.reset_array(&A1_spe);
    std::vector< double > ans1_spe_obs      = counter.count_all();
    std::vector< double > ans1_spe_expected = {
        0, 0, 0,//counter_gains (spe)
        0, 0, 2,//counter_gains (dpl)
        0, 0, 0,//counter_loss (spe)
        1, 1, 0,//counter_loss (dpl)
        0, // counter_overall_changes (spe)
        4 //counter_overall_changes (dpl)
    };

    REQUIRE_THAT(ans1_spe_expected, Catch::Approx(ans1_spe_obs).epsilon(.001));

    /** PhyloArrays to check 2:
     * 
     * // Subfun for both Duplication and speciation
     * 1 | 0 0 1
     * 1 | 0 1 1
     * 0 | 1 0 0
     * 
     */
    
    PhyloArray A2_dpl(3, 3);
    A2_dpl.set_data(new NodeData({1.0, 1.5}, {true, true}, true), true);
    A2_dpl(0, 2) = 1;
    A2_dpl(1, 1) = 1;
    A2_dpl(1, 2) = 1;
    A2_dpl(2, 0) = 1;

    PhyloArray A2_spe(A2_dpl, true);
    A2_spe.D()->duplication = false;
    
    // Adding counters
    PhyloStatsCounter counter2(&A2_dpl);
    counter_gains(counter2.get_counters(), {0u, 1u, 2u}, true);
    counter_gains(counter2.get_counters(), {0u, 1u, 2u}, false);
    counter_loss(counter2.get_counters(), {0u, 1u, 2u}, true);
    counter_loss(counter2.get_counters(), {0u, 1u, 2u}, false);
    counter_overall_changes(counter2.get_counters(), true);
    counter_overall_changes(counter2.get_counters(), false);
    counter_genes_changing(counter2.get_counters(), true);
    counter_genes_changing(counter2.get_counters(), false);
    counter_k_genes_changing(counter2.get_counters(), 1, true);
    counter_k_genes_changing(counter2.get_counters(), 2, true);
    counter_k_genes_changing(counter2.get_counters(), 3, true);
    counter_k_genes_changing(counter2.get_counters(), 1, false);
    counter_k_genes_changing(counter2.get_counters(), 2, false);
    counter_k_genes_changing(counter2.get_counters(), 3, false);


    std::vector< double > ans2_dpl_obs      = counter2.count_all();
    std::vector< double > ans2_dpl_expected = {
        0, 0, 1,//counter_gains (dpl)
        0, 0, 0,//counter_gains (spe)
        2, 1, 0,//counter_loss (dpl)
        0, 0, 0,//counter_loss (spe)
        4, //counter_overall_changes (dpl)
        0, // counter_overall_changes (spe)
        2, // counter_genes_changing (dpl)
        0, // counter_genes_changing (spe)
        0, 1, 0, // 1, 2, 3 genes changing (dpl)
        0, 0, 0  // 1, 2, 3 genes changing (spe)
    };

    REQUIRE_THAT(ans2_dpl_expected, Catch::Approx(ans2_dpl_obs).epsilon(.001));

    counter2.reset_array(&A2_spe);
    std::vector< double > ans2_spe_obs      = counter2.count_all();
    std::vector< double > ans2_spe_expected = {
        0, 0, 0,//counter_gains (spe)
        0, 0, 1,//counter_gains (dpl)
        0, 0, 0,//counter_loss (spe)
        2, 1, 0,//counter_loss (dpl)
        0, // counter_overall_changes (spe)
        4, // counter_overall_changes (spe)
        0, // counter_genes_changing (dpl)
        2, // counter_genes_changing (spe)
        0, 0, 0, // 1, 2, 3 genes changing (dpl)
        0, 1, 0  // 1, 2, 3 genes changing (spe)
    };

    REQUIRE_THAT(ans2_spe_expected, Catch::Approx(ans2_spe_obs).epsilon(.001));
    // { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 4.0, 0.0, 2.0 }
    // { 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 4.0, 0.0, 2.0 }
}

