#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include "aphylomodel.hpp" 

int main() {

    // Defining input data
    /* Setting the following array:
     *   o1 o2   o3 o4
     * a  1  9    0  1
     * b  1  1    0  0
     * c  0  9    9  9
     */

    std::vector< std::vector< uint > > annotations = {
        {0, 0, 1, 9, 9, 9},
        {1, 1, 9, 0, 9, 9},
        {1, 9, 0, 1, 9, 9}
    };

    std::vector< uint > geneid = {0, 1, 2, 3, 4, 5};
    std::vector< uint > parent = {4, 4, 5, 5, 6, 6};

    
    APhyloModel dat(
        annotations, geneid, parent
    );

    // Specifying the terms
    counter_gains(&dat.counters, {0, 1, 2});
    counter_loss(&dat.counters, {0, 1, 2});
    counter_cogain(&dat.counters, 0, 1);
    counter_cogain(&dat.counters, 0, 2);
    counter_cogain(&dat.counters, 1, 2);
    
    // Model parameters
    std::vector< double > par = {
        // Main parameters
        .1, .1, .1, .1, .1, .1, .1, .1, .1,
        // Root probabilities
        .1, .1, .1
        };

    par = {
        27.5028843378887, 118.49940341424, 27.4595445668862, 39.6985832725241,
        58.9044118882134, -0.599039831786341, -43.2677113061893, -30.6160416909189,
        89.0051524467181, -13.9983478226013, -82.7948139476884, -139.317765694915
    };

    dat.init();

    // Stargint to measure time
    auto start = std::chrono::system_clock::now();
    double ans = dat.likelihood(par);
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double,std::milli> diff = end-start;

    printf(      "LogLike        : %.8f\n", std::log(ans));
    std::cout << "Total time (ms): " << diff.count() << std::endl;

    auto simres = dat.simulate(par);

    return 0;
}