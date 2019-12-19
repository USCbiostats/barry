This repository contains a C++ template library that essentially counts sufficient statistics
on binary arrays. The idea of the library is that this can be used together to build exponential
family models as those in Exponential Random Graph Models (ERGMs), but as a generalization that
also deals with non square arrays.

# Design considerations

## Data structures

For start, the main class object should hold the following:

* **The array structure** Right now, we are thinking on dealing with an `std::unordered_map` type
  of structure since search/addition/removel operations have constant average time.
  
* **Pointer to undefined structure** Besides of the graph itself, the data may be acompained
  by other datasets, for example, in the case of genetic annotation we may have the current
  state of some genes, i.e., a binary vector.
  
## Algorithms to implement

* **Counters** Users should be able to define counters using change statistics. From the
  ERGM literature, we know that change statistics can be a very efficient way of counting
  when we have a Markov process. In our case, since we will be doing exhaustive ennumeration,
  a good an efficient way of counting statistics is counting as we add/remove zeros.
  
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
