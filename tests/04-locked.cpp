#include "tests.hpp"

// Defining rule to lock the first cell
RULE_FUNCTION(myrule) {

    return i != 0 || j != 0u;

}

// Defining a change statistics ------------------------------------------------

// Ones
COUNTER_FUNCTION(n_ones) {
    return 1.0;
}

// Zeros
COUNTER_FUNCTION(n_zeros_init) {
    return Array.nrow() * Array.ncol();
}

COUNTER_FUNCTION(n_zeros) {
    return -1.0;
}

// Mutuals
COUNTER_FUNCTION(n_mutual) {
    return (double) (
        (i != j) &
        (!Array.is_empty(j, i))
    );
}

BARRY_TEST_CASE("Constrained support", "[const-support]") {

    // Creating the following array:
    // [  0,]  1  . 
    // [  1,]  .  1 
    barry::BArray<> dat(2,2);
    dat(0,0) = true;
    dat(1,1) = true;

    dat.print();

    barry::Support<> sup(dat);
    
    // Rules
    barry::Rule<> rule1(myrule<>, false);
    sup.add_rule(rule1);

    // Counters
    barry::Counter<> count_n_ones(n_ones<>, nullptr, nullptr, true);
    barry::Counter<> count_n_zeros(n_zeros<>, n_zeros_init<>, nullptr, true);
    barry::Counter<> count_mutuals(n_mutual<>, nullptr, nullptr, true);
    
    sup.add_counter(count_n_ones);
    sup.add_counter(count_n_zeros);
    sup.add_counter(count_mutuals);

    sup.calc();
    auto counts = sup.get_counts();
    sup.print();

    // Expected counts
    barry::MapVec_type<double> expected(0u);
    expected.insert({{1,3,0}, 1});
    expected.insert({{4,0,1}, 1});
    expected.insert({{2,2,0}, 3});
    expected.insert({{3,1,1}, 1});
    expected.insert({{3,1,0}, 2});

    barry::MapVec_type<double> observed(0u);
    size_t n_stats  = 3u;
    size_t n_unique = counts.size() / (n_stats + 1u);
    for (size_t i = 0u; i < n_unique; ++i)
    {
        std::vector< double > tmp(n_stats);
        for (size_t j = 0u; j < n_stats; ++j)
            tmp[j] = counts[i * (n_stats + 1u) + j + 1u];

        observed.insert({tmp, counts[i * (n_stats + 1)]});
    }

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(expected.size() == observed.size());
    REQUIRE(expected[{1,3,0}] == observed[{1,3,0}]);
    REQUIRE(expected[{4,0,1}] == observed[{4,0,1}]);
    REQUIRE(expected[{2,2,0}] == observed[{2,2,0}]);
    REQUIRE(expected[{3,1,1}] == observed[{3,1,1}]);
    REQUIRE(expected[{3,1,0}] == observed[{3,1,0}]);
    #else
    return 0;
    #endif 
    
}