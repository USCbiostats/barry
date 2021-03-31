#include <iostream>
#include <string>
#include "../include/barry/barry.hpp"
#include "catch.hpp"

// Defining rule to lock the first cell
RULE_FUNCTION(myrule) {
    if (i == 0 & j == 0u)
        return true;
    return false;
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

TEST_CASE("Constrained support", "[const-support]") {

    // Creating the following array:
    // [  0,]  1  . 
    // [  1,]  .  1 
    barry::BArray<> dat(2,2);
    dat(0,0) = true;
    dat(1,1) = true;

    dat.print();

    barry::Support<> sup(dat);
    
    // Rules
    barry::Rule<> rule1(myrule<>);
    sup.add_rule(rule1);

    // Counters
    barry::Counter<> count_n_ones(n_ones<>);
    barry::Counter<> count_n_zeros(n_zeros<>, n_zeros_init<>);
    barry::Counter<> count_mutuals(n_mutual<>);
    
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
    for (auto iter = counts.begin(); iter != counts.end(); ++iter)
        observed.insert({iter->first, iter->second});

    REQUIRE(expected.size() == observed.size());
    REQUIRE(expected[{1,3,0}] == observed[{1,3,0}]);
    REQUIRE(expected[{4,0,1}] == observed[{4,0,1}]);
    REQUIRE(expected[{2,2,0}] == observed[{2,2,0}]);
    REQUIRE(expected[{3,1,1}] == observed[{3,1,1}]);
    REQUIRE(expected[{3,1,0}] == observed[{3,1,0}]);
    
}