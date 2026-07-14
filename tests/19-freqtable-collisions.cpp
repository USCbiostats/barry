#include "tests.hpp"

BARRY_TEST_CASE("FreqTable keeps distinct stats separate on hash collisions", "[freqtable]") {

    // Regular usage: equal vectors merge, distinct vectors get their own row
    barry::FreqTable<double> tab;
    std::vector< double > a = {1.0, 2.0, 3.0};
    std::vector< double > b = {4.0, 5.0, 6.0};

    tab.add(a, nullptr);
    tab.add(a, nullptr);
    tab.add(b, nullptr);

    REQUIRE(tab.size() == 2u);

    auto counts = tab.as_vector();
    REQUIRE(counts[0u].first == a);
    REQUIRE(counts[0u].second == 2u);
    REQUIRE(counts[1u].first == b);
    REQUIRE(counts[1u].second == 1u);

    // Forced collisions: the precomputed-hash argument lets us feed
    // different stat vectors under the same 64-bit hash, mimicking a
    // genuine make_hash collision. They must remain separate rows.
    barry::FreqTable<double> tab2;
    std::vector< double > c = {7.0, 8.0, 9.0};
    size_t h = 42u;

    tab2.add(a, &h);
    tab2.add(b, &h);
    tab2.add(a, &h);
    tab2.add(c, &h);
    tab2.add(b, &h);

    REQUIRE(tab2.size() == 3u);

    auto counts2 = tab2.as_vector();
    REQUIRE(counts2[0u].first == a);
    REQUIRE(counts2[0u].second == 2u);
    REQUIRE(counts2[1u].first == b);
    REQUIRE(counts2[1u].second == 2u);
    REQUIRE(counts2[2u].first == c);
    REQUIRE(counts2[2u].second == 1u);

    // Clearing must also reset the collision chains
    tab2.clear();
    REQUIRE(tab2.size() == 0u);

    tab2.add(a, &h);
    tab2.add(b, &h);
    REQUIRE(tab2.size() == 2u);

}
