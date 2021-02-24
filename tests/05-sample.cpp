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

    unsigned int nsamp = 5e4;
    std::vector< netcounters::Network > nets;
    nets.reserve(nsamp);
    netcounters::Network total(4u, 4u);
    for (auto i = 0u; i < nsamp; ++i) {
        nets.push_back(m.sample(0u, par));
        total += nets.at(nets.size() - 1u);
    }

    // Observing the average
    // barry::BArray<double,bool> total(4u, 4u);
    (total/=(double) nsamp).print();
    

}