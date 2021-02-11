#include <iostream>
#include <ostream>
#include "../include/barry.hpp"

typedef std::vector< unsigned int > vuint;

int main() {
  // Creating network of size six with five ties
  netcounters::Network net(
      6, 6,
      {0, 0, 4, 4, 2, 0, 1},
      {1, 2, 0, 2, 4, 0, 1}
  );
  
  // How does this looks like?
  std::cout << "Current view" << std::endl;
  net.print();
  
  // Adding extra ties
  net += {1, 0};
  net(2, 0) = true;
  
  // And removing a couple
  net(0, 0) = false;
  net -= {1, 1};

  std::cout << "New view" << std::endl;  
  net.print();
  
  // Initializing the data. The program deals with freing the memory
  net.set_data(new netcounters::NetworkData, true);

  // Creating counter object for the network and adding stats to count
  netcounters::NetStatsCounter counter(&net);
  netcounters::counter_edges(counter.counters);
  netcounters::counter_ttriads(counter.counters);
  netcounters::counter_isolates(counter.counters);
  netcounters::counter_ctriads(counter.counters);
  netcounters::counter_mutual(counter.counters);
  
  // Counting and printing the results
  std::vector< double > counts = counter.count_all();
  
  std::cout <<
    "Edges             : " << counts[0] << std::endl <<
    "Transitive triads : " << counts[1] << std::endl <<
    "Isolates          : " << counts[2] << std::endl <<
    "C triads          : " << counts[3] << std::endl <<
    "Mutuals           : " << counts[4] << std::endl;
  
  return 0;
}
 
