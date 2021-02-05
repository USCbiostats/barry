#include <vector>

#ifndef BARRY_VECTORS_HPP
#define BARRY_VECTORS_HPP 

/**@brief Compares if -a- and -b- are equal
 * @param a,b Two vectors of the same length
 * @return `true` if all elements are equal.
 */
template <typename T>
inline bool vec_equal(const std::vector< T > & a, const std::vector< T > & b) {
  
  if (a.size() != b.size())
    throw std::length_error("-a- and -b- should have the same length.");
  
  unsigned int i = 0;
  while (a.at(i) == b.at(i++)) {
    if (i == a.size())
      return true;
  }
  
  return false;
}

#endif