#include "cell-bones.hpp"

#ifndef BARRY_CELL_MEAT_HPP
#define BARRY_CELL_MEAT_HPP 1

template <typename Cell_Type>
Cell<Cell_Type>& Cell<Cell_Type>::operator=(Cell<Cell_Type>& other) {
    this->value   = other.value;
    this->visited = other.visited;
    return *this;
}

template <typename Cell_Type>
Cell<Cell_Type>& Cell<Cell_Type>::operator=(Cell<Cell_Type>&& other) noexcept {
    this->value   = std::move(other.value);
    this->visited = std::move(other.visited);
    return *this;
}

template<typename Cell_Type>
bool Cell<Cell_Type>::operator==(const Cell<Cell_Type>& rhs ) const {

    if (this == *rhs)
        return true;
    
    return this->value == rhs.value;

}

template<typename Cell_Type>
bool Cell<Cell_Type>::operator!=(const Cell<Cell_Type>& rhs ) const {

    return !this->operator==(rhs);
    
}


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
template<> inline Cell< uint >::Cell() : value(1u), visited(false) {}
template<> inline Cell< int >::Cell() : value(1), visited(false) {}
template<> inline Cell< bool >::Cell() : value(true), visited(false) {}

#endif