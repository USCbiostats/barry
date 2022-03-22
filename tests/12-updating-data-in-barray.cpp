/**
 * @file 12-updating-data-in-barray.cpp
 * @author George G Vega Yon (g.vegayon@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-06-29
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef BARRY_TEST_MODE_ON
    #include "../include/barry/barry.hpp"
    int main() {
#else

#endif

    typedef barry::BArray<bool, std::vector<int> > array;

    // Creating a barray with data
    array A(4,4);
    A.set_data(new std::vector<int>({1, 2, 3}), true);

    // Creating a copy
    array B(A, true);

    for (auto & b : *B.D_ptr())
        std::cout << b << std::endl;

    std::vector< std::vector<int> > states(
        {{2,3}, {2,8},{6,7}});

    std::vector<int> * tmp = B.D_ptr();

    for (int i = 0; i < static_cast<int>(states.size()); ++i)
    {
        *tmp = states[i];
        for (auto & b : *B.D_ptr())
            std::cout << b << std::endl;

    }



#ifndef BARRY_TEST_MODE_ON

        return 0;

    }

#else

#endif