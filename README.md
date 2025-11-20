[![C/C++ CI](https://github.com/USCbiostats/barry/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/USCbiostats/barry/actions/workflows/c-cpp.yml)
[![Doxygen docs](https://github.com/USCbiostats/barry/actions/workflows/doxy-action.yml/badge.svg)](https://USCbiostats.github.io/barry)
[![codecov](https://codecov.io/gh/USCbiostats/barry/branch/master/graph/badge.svg?token=qGBTD4GJDL)](https://codecov.io/gh/USCbiostats/barry)
[![Integrative Methods of Analysis for Genetic Epidemiology](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)

<h1>Barry: your to-go motif accountant<img src="design/logo.svg" style="max-width:200px;width:50%;" align="right"></h1>

[This repository](https://github.com/USCbiostats/barry) contains a C++ template library that essentially counts sufficient statistics on binary arrays. Its primary goal is to provide a general framework for building discrete exponential-family models. A particular example is Exponential Random Graph Models (ERGMs), but we can use `barry` to deal with non-square arrays.

Among the key features included in `barry`, we have:

* Sparse arrays.
* User-defined count statistics.
* User-defined constrain of the support set.
* Powerset generation of binary arrays.
* Discrete Exponential Family Models module (DEFMs).
* Pooled DEFMs.

To use barry, you can either download the entire repository or, since it is header-only, the single header version [`barry.hpp`](barry.hpp). 

This library was created and maintained by [Dr. George G. Vega Yon](https://ggvy.cl) as part of
his doctoral dissertation ["Essays on Bioinformatics and Social Network Analysis: Statistical and Computational Methods for Complex Systems."](https://digitallibrary.usc.edu/asset-management/2A3BF1WAN5IH)

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
  net.print("Current view");
  
  // Adding extra ties
  net += {1, 0};
  net(2, 0) = true;
  
  // And removing a couple
  net(0, 0) = false;
  net -= {1, 1};

  net.print("New view");
  
  // Initializing the data. The program deals with freing the memory
  net.set_data(new netcounters::NetworkData, true);

  // Creating counter object for the network and adding stats to count
  netcounters::NetStatsCounter counter(&net);
  netcounters::counter_edges(counter.get_counters());
  netcounters::counter_ttriads(counter.get_counters());
  netcounters::counter_isolates(counter.get_counters());
  netcounters::counter_ctriads(counter.get_counters());
  netcounters::counter_mutual(counter.get_counters());
  
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

# Building and Testing

The project includes a Makefile with various build targets. To see all available options, run:

```bash
make help
```

This will display all available targets including:
- `tests` - Build and run tests with OpenMP
- `test-portable` - Build and run tests without OpenMP (cross-platform)
- `coverage` - Run tests with coverage analysis
- `docs` - Generate documentation

The default target is `help`, so running `make` without arguments will show the available options.

# Features

## Efficient memory usage

One of the key features of `barry` is that it will handle memory efficiently. In the case of pooled-data models, the module for statistical models avoids double-counting support when possible by keeping track of what datasets (networks, for instance) share the same.

<div align="center">
<img src="design/ergm-computing.svg">
</div>


# Documentation

More information can be found in the Doxygen website [here](https://uscbiostats.github.io/barry) and in the PDF version
of the documentation [here](https://github.com/USCbiostats/barry/blob/gh-pages/latex/refman.pdf).

# Code of Conduct

Please note that the `barry` project is released with a
[Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

