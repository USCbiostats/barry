#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include "leaf.hpp"

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

    // Stargint to measure time
    auto start = std::chrono::system_clock::now();

    // Specifying the terms
    PhyloCounters counters;
    counter_gains(&counters, {0, 1, 2});
    counter_loss(&counters, {0, 1, 2});
    counter_cogain(&counters, 0, 1);
    counter_cogain(&counters, 0, 2);
    counter_cogain(&counters, 1, 2);

    Leafs dat(
        annotations, geneid, parent, counters
    );
    
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;

    std::cout << "Total time: " << diff.count() << std::endl;

    dat.print();    
    

    return 0;
}