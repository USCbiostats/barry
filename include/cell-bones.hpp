// #include <vector>
#include "typedefs.hpp"

#ifndef CELL_BONES_HPP 
#define CELL_BONES_HPP 1

/**@brief Entries in BArray.
 * For now, it only has two members:
 * - value: the content
 * - visited: boolean (just a convenient)
 */
template <class Cell_Type > class Cell {
public:
  Cell_Type value;
  bool visited;
  Cell();
  Cell(Cell_Type value_, bool visited_ = false) :
    value(value_), visited(visited_) {};
  ~Cell() {};
  
  // Copy by-reference constructor
  Cell(Cell& arg) : value(arg.value), visited(arg.visited) {};
  
  // This is an explicit declaration since in other cases it seems
  // to try to use the move operator, which I do not intent to use.
  Cell(const Cell& arg) : value(arg.value), visited(arg.visited) {};
  
  // Copy by assignment
  Cell& operator=(Cell& other) {
    this->value   = other.value;
    this->visited = other.visited;
    return *this;
  };
  
  // Move constructor
  Cell(Cell&& arg):
    value(std::move(arg.value)),
    visited(std::move(arg.visited)) {} ;
  
  // Move assign operator
  Cell& operator=(Cell&& other) {
    this->value   = std::move(other.value);
    this->visited = std::move(other.visited);
    return *this;
  };
  
  void add(Cell_Type x);
  
  // Casting operator (implicit and explicit)
  // int x = Cell<int>(1); // returns 1
  operator Cell_Type() const {return this->value;};
  
};

/***
 * Specializations
 */

template <> inline void Cell<double>::add(double x) {
  value += x;
  return;
}

template <> inline void Cell<unsigned int>::add(unsigned int x) {
  value += x;
  return;
}

template <> inline void Cell<int>::add(int x) {
  value += x;
  return;
}

template <> inline void Cell<bool>::add(bool x) {
  value = true;
  return;
}

template<> inline Cell< double >::Cell() : value(1.0), visited(false) {}
template<> inline Cell< bool >::Cell() : value(true), visited(false) {}
#endif
