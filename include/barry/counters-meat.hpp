#include "counters-bones.hpp"

#ifndef BARRY_COUNTERS_MEAT_HPP
#define BARRY_COUNTERS_MEAT_HPP 1

template <typename Array_Type, typename Data_Type>
inline Counter<Array_Type,Data_Type>::Counter(
    const Counter<Array_Type,Data_Type> & counter_
) : count_fun(counter_.count_fun), init_fun(counter_.init_fun) {

    if (counter_.delete_data) 
    {

        this->data = new Data_Type(*counter_.data);
        delete_data = true;
        
    } 
    else
    {

        this->data = counter_.data;
        delete_data = false;

    }

    this->name = counter_.name;
    this->desc = counter_.desc;

    return;

}

template <typename Array_Type, typename Data_Type>
Counter<Array_Type,Data_Type> Counter<Array_Type,Data_Type>::operator=(
    const Counter<Array_Type,Data_Type> & counter_
)
{

    if (this != &counter_) {

        count_fun = counter_.count_fun;
        init_fun  = counter_.init_fun;

        if (counter_.delete_data)
        {
            this->data = new Data_Type(*counter_.data);
            delete_data = true;
        } 
        else
        {

            this->data = counter_.data;
            delete_data = false;

        }

        this->name = counter_.name;
        this->desc = counter_.desc;

    }

    return *this;

}

template <typename Array_Type, typename Data_Type>
inline Counters<Array_Type,Data_Type>::Counters(
    const Counters<Array_Type,Data_Type> & counter_
)
{

    // Checking which need to be deleted
    std::vector< bool > tbd(counter_.size(), false);
    for (auto& i : counter_.to_be_deleted)
        tbd[i] = true;

    // Copy all counters, if a counter is tagged as 
    // to be deleted, then copy the value
    for (auto i = 0u; i != counter_.size(); ++i)
    {

        if (tbd[i])
            this->add_counter(*counter_.data[i]);
        else
            this->add_counter(counter_.data[i]);

    }

    return;

}

template <typename Array_Type, typename Data_Type>
Counters<Array_Type,Data_Type> Counters<Array_Type,Data_Type>::operator=(
    const Counters<Array_Type,Data_Type> & counter_
) {

    if (this != &counter_) {

        // Checking which need to be deleted
        std::vector< bool > tbd(counter_.size(), false);
        for (auto i : counter_.to_be_deleted)
            tbd[i] = true;

        // Copy all counters, if a counter is tagged as 
        // to be deleted, then copy the value
        for (uint i = 0u; i != counter_.size(); ++i)
        {

            if (tbd[i])
                this->add_counter(*counter_.data[i]);
            else
                this->add_counter(counter_.data[i]);
            
        }

    }

    return *this;

}

template <typename Array_Type, typename Data_Type>
inline double Counter<Array_Type, Data_Type>::count(
    Array_Type & Array, uint i, uint j)
{

    if (count_fun == nullptr)
        return 0.0;

    return count_fun(Array, i, j, data);

}

template <typename Array_Type, typename Data_Type>
inline double Counter<Array_Type, Data_Type>::init(
    Array_Type & Array, uint i, uint j
)
{

    if (init_fun == nullptr)
        return 0.0;

    return init_fun(Array, i, j, data);

}

template <typename Array_Type, typename Data_Type>
inline Counter<Array_Type,Data_Type> &
Counters<Array_Type,Data_Type>::operator[](uint idx) {

    return *data[idx];

}

template <typename Array_Type, typename Data_Type>
inline void Counters<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> & counter
)
{
    
    to_be_deleted.push_back(data.size());
    data.push_back(new Counter<Array_Type, Data_Type>(counter));
    
    return;
}

template <typename Array_Type, typename Data_Type>
inline void Counters<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> * counter
)
{
    
    data.push_back(counter);
    return;
    
}

template <typename Array_Type, typename Data_Type>
inline void Counters<Array_Type,Data_Type>::add_counter(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_,
    Data_Type *                            data_,
    bool                                   delete_data_,
    std::string                            name_,
    std::string                            desc_
)
{
  
    /* We still need to delete the counter since we are using the 'new' operator.
      * Yet, the actual data may not need to be deleted.
      */
    to_be_deleted.push_back(data.size());
    
    data.push_back(new Counter<Array_Type,Data_Type>(
        count_fun_,
        init_fun_,
        data_,
        delete_data_,
        name_,
        desc_
    ));
  
    return;
    
}

template <typename Array_Type, typename Data_Type>
inline void Counters<Array_Type,Data_Type>::clear()
{
    
    for (auto& i : to_be_deleted)
        delete data[i];
    
    to_be_deleted.clear();
    
    return;
    
}

#endif 