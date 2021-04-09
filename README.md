[![C/C++ CI](https://github.com/USCbiostats/barry/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/USCbiostats/barry/actions/workflows/c-cpp.yml)
[![Doxygen docs](https://github.com/USCbiostats/barry/actions/workflows/doxy-action.yml/badge.svg)](https://USCbiostats.github.io/barry)
[![Integrative Methods of Analysis for Genetic Epidemiology](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)

<h1>Barry: your to-go motif accountant<img src="logo.svg" style="max-width:200px;width:50%;" align="right"></h1>

[This repository](https://github.com/USCbiostats/barry) contains a C++ template
library that essentially counts sufficient statistics on binary arrays. The idea
of the library is that this can be used together to build exponential family
models as those in Exponential Random Graph Models (ERGMs), but as a
generalization that also deals with non square arrays.

# Examples

## Counting statistics in a graph

In the following code we create an array of size 5x5 of class `Network`
(available in the namespace [netcounters](https://uscbiostats.github.io/barry/namespacebarry_1_1counters_1_1network.html)), add/remove ties, print the
graph, and count common statistics used in ERGMs:

```cpp
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
```

Compiling this program using g++

```bash
g++ -std=c++11 -Wall -pedantic 08-counts.cpp -o counts && ./counts
```

Yields the following output:

```bash
Current view
[  0,]  1  1  1  .  .  . 
[  1,]  .  1  .  .  .  . 
[  2,]  .  .  .  .  1  . 
[  3,]  .  .  .  .  .  . 
[  4,]  1  .  1  .  .  . 
[  5,]  .  .  .  .  .  . 
New view
[  0,]  .  1  1  .  .  . 
[  1,]  1  .  .  .  .  . 
[  2,]  1  .  .  .  1  . 
[  3,]  .  .  .  .  .  . 
[  4,]  1  .  1  .  .  . 
[  5,]  .  .  .  .  .  . 
Edges             : 7
Transitive triads : 3
Isolates          : 2
C triads          : 1
Mutuals           : 3
```

## Code of Conduct

Please note that the `barry` project is released with a
[Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

