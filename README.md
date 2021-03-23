[![Build Status](https://travis-ci.com/USCbiostats/barry.svg?branch=master)](https://travis-ci.com/USCbiostats/barry) [![Integrative Methods of Analysis for Genetic Epidemiology](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)

<h1>Barry: your to-go motif accountant<img src="logo.svg" style="max-width:300px;width:50%;"></h1>


This repository contains a C++ template library that essentially counts sufficient statistics
on binary arrays. The idea of the library is that this can be used together to build exponential
family models as those in Exponential Random Graph Models (ERGMs), but as a generalization that
also deals with non square arrays.

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


<!--
# Design considerations (old)

## Data structures

For start, the main class object should hold the following:

* **The array structure** Right now, we are thinking on dealing with an `std::unordered_map` type
  of structure since search/addition/removel operations have constant average time.
  
* **Pointer to undefined structure** Besides of the graph itself, the data may be acompained
  by other datasets, for example, in the case of genetic annotation we may have the current
  state of some genes, i.e., a binary vector.

* **Structural constrains** Having arrays specifying the blocked spaces of the array. This,
  in principle, could affect all operations related to modifying the array, e.g. if the
  pair `(i, j)` is blocked, then no addition/deletion can be done on that respect.

  The structural constrains may be better reflected as a counterpart: free blocks. This
  way, any algorithm that needs to iterate through cells that can be changed can use its
  counter part. 

  Enumeration of both sets would have in total $n \times m$ unordered pairs. One problem of this
  is the fact that this type of data structure is unefficient as it can grow too fast. Yet,
  the whole idea of this C++ library is to be able to fully enumerate support of arrays, so
  problems that need to deal with larger datasets may not be suitable for this approach.
  
## Algorithms to implement

* **Counters** Users should be able to define counters using change statistics. From the
  ERGM literature, we know that change statistics can be a very efficient way of counting
  when we have a Markov process. In our case, since we will be doing exhaustive ennumeration,
  a good an efficient way of counting statistics is counting as we add/remove zeros.
  
  Counters should have the following implementation:
  
  ```cpp
  double counter_[name](const Array & x, unsigned int row, unsigned int col) {
    // Visit neighbors
    Array::local_iterator iter(x, row, col);
    double counts = 0;
    for (auto i = iter.begin(); i != iter.end(); ++i) {
      ...[do your thing]...
      counts += ...
    }
    
    return counts;
  }
  ```
  
* **Constrained Exhaust enumeration** Exhaust enumeration can be done using recursive algorithms
  activating and deactivating cells in the array. One nice feature would be to allow users
  to specify constrained, essentially blocked, cells in the array. This would go together with 
  the counters.
  
  The constrains can just be `std::unordered_map` objects in which the coordinates of the 
  cells that need to be blocked are specified. Some standard constrains can be:
  
  - Blocks (ranges).
  - Symmetry (in the case of square, undirected graphs).
  
  Furthermore, we should, at least in principle, allow the user to speficy default values
  for the blocked cells (0/1).

* **Array changes** Addition and deletion of 0/1 states. If we use `std::unordered_map`
  this should be straight forward. Perhaps just making an alias or binary operator, e.g.

  ```cpp
  X += (i, j);
  X.add(i, j);
  X.rm(i, j);
  ```

* **Markov Chains** This would be nice to have, but not necesary for now. The idea
  is to have various types of algorithms to implement transitions to new states, 
  for example

  - Randomly adding/removing new pairs.
  - Unconstrained endpoints-switching/swap.
  - Constrained endpoints-swap.

  The constrained component can be related to the graph constrains specified part 
  of the countner.

* **Estimate support size** Do this conditioning on the constrains. This should be
  rather straight forward. The support of the set should be defined by
  
  $$
  2^{(n\times m - |blocked|)}
  $$

  Where $blocked$ is the set of blocked cells.

To-do: https://github.com/codecov/example-cpp11/blob/master/.travis.yml

-->