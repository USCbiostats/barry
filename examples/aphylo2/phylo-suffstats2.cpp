#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include "aphylomodel.hpp" 

template <class T>
inline void print(std::vector< T > & x) {
  std::cout << "[ " ;
  for (auto i = x.begin(); i != x.end(); ++i)
    std::cout << *i << " ";

  std::cout << "]\n";

  return;
}

int main() {

    // Defining input data
    /* Setting the following array:
     *   o1 o2   o3 o4
     * a  1  9    0  1
     * b  1  1    0  0
     * c  0  9    9  9
     */

    std::vector< std::vector< uint > > annotations = {
        {1, 1, 0},
        {9, 1, 0},
        {0, 9, 1},
        {1, 0, 9},
        {9, 9, 9},
        {9, 9, 9},
        {9, 9, 9}
        // {0, 0, 1, 9, 9, 9},
        // {1, 1, 9, 0, 9, 9},
        // {1, 9, 0, 1, 9, 9}
    };

    std::vector< uint > geneid = {0, 1, 2, 3, 4, 5};
    std::vector< uint > parent = {4, 4, 5, 5, 6, 6};

    
    Geese dat(
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
    print(simres[0u]);
    print(simres[1u]);
    print(simres[2u]);
    print(simres[3u]);
    print(simres[4u]);
    print(simres[5u]);
    print(simres[6u]);

    // More interesting experiment
    std::vector< std::vector<uint> > ann2;
    ann2.resize(7u, {9, 9, 9});
    std::vector< uint > geneid2 = {0u, 1u, 2u, 3u, 4u, 5u};
    std::vector< uint > parent2 = {6u, 4u, 4u, 5u, 5u, 6u};
    Geese model2(ann2, geneid2, parent2);

    // Adding terms
    ::counter_cogain(&model2.counters, 0, 1);
    ::counter_cogain(&model2.counters, 0, 2);
    ::counter_cogain(&model2.counters, 1, 2);
    ::counter_subfun(&model2.counters, 1, 2);
    ::counter_subfun(&model2.counters, 1, 0);
    ::counter_subfun(&model2.counters, 0, 2);
    ::counter_maxfuns(&model2.counters, 2, 2);

    model2.init();
    model2.set_seed(1121);

    std::cout << "In this model, root has no functions, ";
    std::cout << "Pairs (0, 1) and (0, 2) are gained, but ";
    std::cout << "Pair (1, 2) is never observed" << std::endl;
    auto simres2 = model2.simulate({100, 100, -100, -100, -100, 100, 100, -100, -100, -100});
    print(simres2[0u]);
    print(simres2[1u]);
    print(simres2[2u]);
    print(simres2[3u]);
    print(simres2[4u]);
    print(simres2[5u]);
    print(simres2[6u]);

    return 0;
}