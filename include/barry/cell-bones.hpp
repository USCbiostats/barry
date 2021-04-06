// #include <vector>
#include "typedefs.hpp"

#ifndef BARRY_CELL_BONES_HPP 
#define BARRY_CELL_BONES_HPP 1

/**
 * @brief Entries in BArray.
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
    
    // This is an explicit declaration since in other cases it seems
    // to try to use the move operator, which I do not intent to use.
    Cell(const Cell<Cell_Type>& arg) : value(arg.value), visited(arg.visited) {};
    
    // Copy by assignment
    Cell<Cell_Type>& operator=(Cell<Cell_Type>& other);
    
    // Move constructor
    Cell(Cell<Cell_Type>&& arg) noexcept:
        value(std::move(arg.value)),
        visited(std::move(arg.visited)) {} ;
    
    // Move assign operator
    Cell<Cell_Type>& operator=(Cell<Cell_Type>&& other) noexcept;
    
    void add(Cell_Type x);
    
    // Casting operator (implicit and explicit)
    // int x = Cell<int>(1); // returns 1
    operator Cell_Type() const {return this->value;};
  
};

#endif
