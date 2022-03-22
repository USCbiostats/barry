// #include <vector>
// #include <stdexcept>

#include "rules-bones.hpp"

#ifndef BARRY_RULES_MEAT_HPP
#define BARRY_RULES_MEAT_HPP 1

template <typename Array_Type, typename Data_Type>
inline Rules<Array_Type,Data_Type>::Rules(
    const Rules<Array_Type,Data_Type> & rules_
) {

    // Copy all rules, if a rule is tagged as 
    // to be deleted, then copy the value
    for (auto i = 0u; i != rules_.size(); ++i)
        this->add_rule(rules_.data[i]);

    return;

}

template <typename Array_Type, typename Data_Type>
Rules<Array_Type,Data_Type> Rules<Array_Type,Data_Type>::operator=(
    const Rules<Array_Type,Data_Type> & rules_
) {

    if (this != &rules_) {

        // Copy all rules, if a rule is tagged as 
        // to be deleted, then copy the value
        for (auto i = 0u; i != rules_.size(); ++i)
            this->add_rule(rules_.data[i]);

    }

    return *this;

}

template<typename Array_Type, typename Data_Type>
inline bool Rule<Array_Type,Data_Type>::operator()(const Array_Type & a, uint i, uint j) {
    return fun(a, i, j, dat);
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::add_rule(
        Rule<Array_Type, Data_Type> rule
) {
    
    data.push_back(rule);
    
    return;
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::add_rule(
        Rule_fun_type<Array_Type,Data_Type> rule_,
        Data_Type                           data_
) {
       
    data.push_back(Rule<Array_Type,Data_Type>(
        rule_,
        data_
    ));
    
    return;
    
}

template <typename Array_Type, typename Data_Type>
inline bool Rules<Array_Type,Data_Type>::operator()(
    const Array_Type & a, uint i, uint j
) {
    
    if (data.size()==0u)
        return true;
    
    for (auto & f: data)
        if (!f.operator()(a, i, j))
            return false;
    
    return true;
    
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::get_seq(
    const Array_Type & a,
    std::vector< size_t > * free,
    std::vector< size_t > * locked
) {

    
    uint N = a.nrow();
    uint K = a.ncol();
    
    // Reserving some space
    (void) free->empty();
    (void) free->reserve(2u * N * K);
    
    for (uint i = 0u; i < N; ++i)
    {

        for (uint j = 0u; j < K; ++j)
        {

            // Locked cells are skipped
            if (!this->operator()(a, i, j))
            {

                if (locked != nullptr)
                {

                    locked->push_back(i);
                    locked->push_back(j);

                }

                continue;

            }

            free->push_back(i);
            free->push_back(j);
                
        }

    }
    
    free->shrink_to_fit();

    return;

}

#endif
