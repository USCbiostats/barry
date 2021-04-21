#include "counters-bones.hpp"

#ifndef BARRY_COUNTERS_MEAT_HPP
#define BARRY_COUNTERS_MEAT_HPP 1

template <typename Array_Type, typename Data_Type>
inline Counters<Array_Type,Data_Type>::Counters() {
    this->data = new std::vector<Counter<Array_Type,Data_Type>*>(0u);
    this->delete_data = true;
}

template <typename Array_Type, typename Data_Type>
inline Counter<Array_Type,Data_Type>::Counter(
    const Counter<Array_Type,Data_Type> & counter_
) : count_fun(counter_.count_fun), init_fun(counter_.init_fun) {

    if (counter_.delete_data) 
    {

        this->data = new Data_Type(*counter_.data);
        delete_data = true;
        
    } else {

        this->data = counter_.data;
        delete_data = false;

    }

    this->name = counter_.name;
    this->desc = counter_.desc;

    return;

}


template <typename Array_Type, typename Data_Type>
inline Counter<Array_Type,Data_Type>::Counter(
    Counter<Array_Type,Data_Type> && counter_
    ) noexcept :
    count_fun(std::move(counter_.count_fun)),
    init_fun(std::move(counter_.init_fun)),
    data(std::move(counter_.data)),
    delete_data(counter_.delete_data),
    name(std::move(counter_.name)),
    desc(std::move(counter_.desc))
{

    counter_.delete_data = false;
    counter_.data = nullptr;

} ///< Move constructor

template <typename Array_Type, typename Data_Type>
inline Counter<Array_Type,Data_Type> Counter<Array_Type,Data_Type>::operator=(
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
inline Counter<Array_Type,Data_Type> & Counter<Array_Type,Data_Type>::operator=(
    Counter<Array_Type,Data_Type> && counter_
) noexcept {

    if (this != &counter_) {

        // Data
        if (delete_data)
            delete data;
        
        this->data           = counter_.data;
        this->delete_data    = std::move(counter_.delete_data);
        counter_.data        = nullptr;
        counter_.delete_data = false;

        // Functions
        this->count_fun = std::move(counter_.count_fun);
        this->init_fun = std::move(counter_.init_fun);

        // Descriptions
        this->name = std::move(counter_.name);
        this->desc = std::move(counter_.desc);

    }

    return *this;

} ///< Move assignment

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

    return *(data->operator[](idx));

}

template <typename Array_Type, typename Data_Type>
inline Counters<Array_Type,Data_Type>::Counters(
    const Counters<Array_Type,Data_Type> & counter_
) :
    data(new std::vector< Counter<Array_Type,Data_Type>* >(0u)),
    to_be_deleted(0u),
    delete_data(true)
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
            this->add_counter(*counter_.data->operator[](i));
        else
            this->add_counter(counter_.data->operator[](i));

    }

    return;

}

template <typename Array_Type, typename Data_Type>
inline Counters<Array_Type,Data_Type>::Counters(
    Counters<Array_Type,Data_Type> && counters_
    ) noexcept :
    data(std::move(counters_.data)),
    to_be_deleted(std::move(counters_.to_be_deleted)),
    delete_data(counters_.delete_data)
{

    // Taking care of memory
    counters_.data = nullptr;
    counters_.delete_data = false;

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

        // Removing the data currently stored in the 
        this->clear();
        this->to_be_deleted.resize(0u);

        data = new std::vector< Counter<Array_Type,Data_Type>* >(0u);
        delete_data = true;

        // Copy all counters, if a counter is tagged as 
        // to be deleted, then copy the value
        for (uint i = 0u; i != counter_.size(); ++i)
        {

            if (tbd[i])
                this->add_counter(*counter_.data->operator[](i));
            else
                this->add_counter(counter_.data->operator[](i));
            
        }

    }

    return *this;

}

template <typename Array_Type, typename Data_Type>
inline Counters<Array_Type,Data_Type> & Counters<Array_Type,Data_Type>::operator=(
    Counters<Array_Type,Data_Type> && counters_
    ) noexcept 
{

    if (this != &counters_)
    {        

        // Removing the data currently stored in the 
        this->clear();

        data          = std::move(counters_.data);
        delete_data   = counters_.delete_data;
        to_be_deleted = std::move(counters_.to_be_deleted);

        counters_.data = nullptr;
        counters_.delete_data = false;

    }

    return *this;

}

template <typename Array_Type, typename Data_Type>
inline void Counters<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> & counter
)
{
    
    to_be_deleted.push_back(data->size());
    data->push_back(new Counter<Array_Type, Data_Type>(counter));
    
    return;
}

template <typename Array_Type, typename Data_Type>
inline void Counters<Array_Type,Data_Type>::add_counter(
    Counter<Array_Type, Data_Type> * counter
)
{
    
    data->push_back(counter);
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
    to_be_deleted.push_back(data->size());
    
    data->push_back(new Counter<Array_Type,Data_Type>(
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
        delete data->operator[](i);
    
    to_be_deleted.clear();

    if (delete_data)
        delete data;

    data = nullptr;
    
    return;
    
}

#endif 