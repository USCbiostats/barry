#ifndef H_BARRY_TESTS
#define H_BARRY_TESTS

#include "../include/barry/barry.hpp"

#include <vector>
#include <iostream>
#include <type_traits>
#include <exception>

template <class T>
inline void print(std::vector< T > & x)
{

    std::cout << "[ " ;
    for (auto i = x.begin(); i != x.end(); ++i)
        std::cout << *i << " ";

    std::cout << "]\n";

    return;

}

template <class T>
inline std::vector< T > vabsdiff(std::vector< T > & a, std::vector< T > & b)
{
  
    std::vector< T > ans(a.size());
    if (a.size() != b.size())
        throw std::length_error("-a- and -b- should be of th same length.");
    
    for (auto i = 0u; i < a.size(); ++i)
        ans[i] = std::fabs(a[i] - b[i]);
      
    return ans;

}

#ifndef CATCH_CONFIG_MAIN

    
    #define BARRY_TEST_CASE(desc, tag) \
        int main() 
    
#else

    #define BARRY_TEST_CASE(desc, tag) \
        TEST_CASE(desc, tag) 

#endif


#endif
