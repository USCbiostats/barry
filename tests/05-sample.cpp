// #include "../include/barry/barry.hpp"
// #include "catch.hpp"

TEST_CASE("Sampling networks (with Model)", "[sampling w model]") {


    // Builing model    
    netcounters::NetModel m;
    m.store_psets();

    netcounters::counter_edges(&m.counters);
    netcounters::counter_mutual(&m.counters);

    // Adding rules
    netcounters::rules_zerodiag(&m.rules);

    // Adding an array
    netcounters::Network net(4, 4);
    net.set_data(new netcounters::NetworkData({0,0,1,0}), true);
    m.add_array(net);

    // Sampling
    std::vector<double> par = {-2, 2*0};
    m.set_seed(133);

    unsigned int nsamp = 1e4;
    std::vector< netcounters::Network > nets;
    nets.reserve(nsamp);
    netcounters::Network total(4u, 4u);
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