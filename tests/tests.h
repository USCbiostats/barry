#include <vector>
#include <iostream>
#include <type_traits>

#ifndef H_PRUNER_TESTS
#define H_PRUNER_TESTS

template <class T>
inline void print(std::vector< T > & x) {
  std::cout << "[ " ;
  for (auto i = x.begin(); i != x.end(); ++i)
    std::cout << *i << " ";

  std::cout << "]\n";

  return;
}
#endif
