// #include "../include/barry/barry.hpp"
// #include "catch.hpp"

TEST_CASE("Sampling networks (with Model)", "[sampling w model]") {

    using namespace barry::counters::network;

    // Builing model    
    NetModel<Network> m;
    m.store_psets();

    counter_edges<>(m.get_counters());
    counter_mutual<>(m.get_counters());

    // Adding rules
    rules_zerodiag<>(m.get_rules());

    // Adding an array
    Network net(4, 4);
    net.set_data(new NetworkData({0,0,1,0}), true);
    m.add_array(net);

    // Sampling
    std::vector<double> par = {-2, 2*0};
    m.set_seed(133);

    unsigned int nsamp = 1e4;
    std::vector< Network > nets;
    nets.reserve(nsamp);
    Network total(4u, 4u);
    for (auto i = 0u; i < nsamp; ++i) {
        nets.push_back(m.sample(0u, par));
        total += nets.at(nets.size() - 1u);
    }

    // Observing the average
    std::cout << "The average per cell is:" << std::endl;
    (total/=(double) nsamp).print();

    std::cout << "The value of cell (1, 0): " << (double) total(1,0) << std::endl;
    REQUIRE(std::abs(
        static_cast<double>(total(1,0)) - 0.1192029
    ) < .1);

}