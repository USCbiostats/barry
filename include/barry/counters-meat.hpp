#include "counters-bones.hpp"

#ifndef BARRY_COUNTERS_MEAT_HPP
#define BARRY_COUNTERS_MEAT_HPP 1

#define COUNTER_TYPE() Counter<Array_Type,Data_Type>

#define COUNTER_TEMPLATE_ARGS() <typename Array_Type, typename Data_Type>

#define COUNTER_TEMPLATE(a,b) \
    template COUNTER_TEMPLATE_ARGS() inline a COUNTER_TYPE()::b

COUNTER_TEMPLATE(,Counter)(
    const Counter<Array_Type,Data_Type> & counter_
) : count_fun(counter_.count_fun), init_fun(counter_.init_fun) {

    if (counter_.delete_data) 
    {

        this->data = new Data_Type(*counter_.data);
        this->delete_data = true;
        
    }
    else
    {

        this->data = counter_.data;
        this->delete_data = false;

    }

    this->name = counter_.name;
    this->desc = counter_.desc;

    return;

}


COUNTER_TEMPLATE(,Counter)(
    Counter<Array_Type,Data_Type> && counter_
    ) noexcept :
    count_fun(std::move(counter_.count_fun)),
    init_fun(std::move(counter_.init_fun)),
    data(std::move(counter_.data)),
    delete_data(std::move(counter_.delete_data)),
    name(std::move(counter_.name)),
    desc(std::move(counter_.desc))
{

    counter_.data = nullptr;
    counter_.delete_data = false;

} ///< Move constructor

COUNTER_TEMPLATE(COUNTER_TYPE(),operator=)(
    const Counter<Array_Type,Data_Type> & counter_
)
{

    if (this != &counter_) {

        this->count_fun = counter_.count_fun;
        this->init_fun = counter_.init_fun;

        if (counter_.delete_data) 
        {

            this->data = new Data_Type(*counter_.data);
            this->delete_data = true;
            
        } else {

            this->data = counter_.data;
            this->delete_data = false;

        }

        this->name = counter_.name;
        this->desc = counter_.desc;

    }

    return *this;

}

COUNTER_TEMPLATE(COUNTER_TYPE() &,operator=)(
    Counter<Array_Type,Data_Type> && counter_
) noexcept {

    if (this != &counter_)
    {

        // Data
        if (delete_data)
            delete data;
        
        this->data        = std::move(counter_.data);
        this->delete_data = std::move(counter_.delete_data);

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

COUNTER_TEMPLATE(double, count)(Array_Type & Array, uint i, uint j)
{

    if (count_fun == nullptr)
        return 0.0;

    return count_fun(Array, i, j, data);

}

COUNTER_TEMPLATE(double, init)(Array_Type & Array, uint i, uint j)
{

    if (init_fun == nullptr)
        return 0.0;

    return init_fun(Array, i, j, data);

}

COUNTER_TEMPLATE(std::string, get_name)() const {
    return this->name;
}

COUNTER_TEMPLATE(std::string, get_description)() const {
    return this->name;
}

////////////////////////////////////////////////////////////////////////////////
// Counters
////////////////////////////////////////////////////////////////////////////////

#define COUNTERS_TYPE() Counters<Array_Type,Data_Type>

#define COUNTERS_TEMPLATE_ARGS() <typename Array_Type, typename Data_Type>

#define COUNTERS_TEMPLATE(a,b) \
    template COUNTERS_TEMPLATE_ARGS() inline a COUNTERS_TYPE()::b

COUNTERS_TEMPLATE(, Counters)() {
    this->data = new std::vector<Counter<Array_Type,Data_Type>*>(0u);
    this->to_be_deleted = new std::vector< uint >(0u);
    this->delete_data = true;
    this->delete_to_be_deleted = true;
}

COUNTERS_TEMPLATE(COUNTER_TYPE() &, operator[])(uint idx) {

    return *(data->operator[](idx));

}

COUNTERS_TEMPLATE(, Counters)(const Counters<Array_Type,Data_Type> & counter_) :
    data(new std::vector< Counter<Array_Type,Data_Type>* >(0u)),
    to_be_deleted(new std::vector< uint >(0u)),
    delete_data(true),
    delete_to_be_deleted(true)
{

    // Checking which need to be deleted
    std::vector< bool > tbd(counter_.size(), false);
    for (auto& i : *(counter_.to_be_deleted))
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

COUNTERS_TEMPLATE(, Counters)(Counters<Array_Type,Data_Type> && counters_) noexcept :
    data(std::move(counters_.data)),
    to_be_deleted(std::move(counters_.to_be_deleted)),
    delete_data(std::move(counters_.delete_data)),
    delete_to_be_deleted(std::move(counters_.delete_to_be_deleted))
{

    // Taking care of memory
    counters_.data = nullptr;
    counters_.to_be_deleted = nullptr;
    counters_.delete_data = false;
    counters_.delete_to_be_deleted = false;

}

COUNTERS_TEMPLATE(COUNTERS_TYPE(), operator=)(const Counters<Array_Type,Data_Type> & counter_) {

    if (this != &counter_) {

        // Checking which need to be deleted
        std::vector< bool > tbd(counter_.size(), false);
        for (auto i : *(counter_.to_be_deleted))
            tbd[i] = true;

        // Removing the data currently stored in the 
        this->clear();
        

        data = new std::vector< Counter<Array_Type,Data_Type>* >(0u);
        to_be_deleted = new std::vector< uint >(0u);
        delete_data = true;
        delete_to_be_deleted = true;

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

COUNTERS_TEMPLATE(COUNTERS_TYPE() &, operator=)(Counters<Array_Type,Data_Type> && counters_) noexcept 
{

    if (this != &counters_)
    {        

        // Removing the data currently stored in the 
        this->clear();

        data          = std::move(counters_.data);
        to_be_deleted = std::move(counters_.to_be_deleted);
        delete_data   = std::move(counters_.delete_data);
        delete_to_be_deleted   = std::move(counters_.delete_to_be_deleted);

        counters_.data = nullptr;
        counters_.to_be_deleted = nullptr;
        counters_.delete_data = false;
        counters_.delete_to_be_deleted = false;

    }

    return *this;

}

COUNTERS_TEMPLATE(void, add_counter)(Counter<Array_Type, Data_Type> & counter)
{
    
    to_be_deleted->push_back(data->size());
    data->push_back(new Counter<Array_Type, Data_Type>(counter));
    
    return;
}

COUNTERS_TEMPLATE(void, add_counter)(Counter<Array_Type, Data_Type> * counter)
{
    
    data->push_back(counter);
    return;
    
}

COUNTERS_TEMPLATE(void, add_counter)(
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
    to_be_deleted->push_back(data->size());
    
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

COUNTERS_TEMPLATE(void, clear)()
{
        
    if (delete_to_be_deleted)
    {

        #pragma GCC ivdep
        for (auto& i : (*to_be_deleted))
            delete data->operator[](i);

        delete to_be_deleted;

    }

    if (delete_data)
        delete data;

    data          = nullptr;
    to_be_deleted = nullptr;
    
    return;
    
}

COUNTERS_TEMPLATE(std::vector<std::string>, get_names)() const
{

    std::vector< std::string > out(this->size());
    for (unsigned int i = 0u; i < out.size(); ++i)
        out[i] = this->data->at(i)->get_name();

    return out;

}

COUNTERS_TEMPLATE(std::vector<std::string>, get_descriptions)() const
{
    
    std::vector< std::string > out(this->size());
    for (unsigned int i = 0u; i < out.size(); ++i)
        out[i] = this->data->at(i)->get_description();

    return this->name;

}

#undef COUNTER_TYPE
#undef COUNTER_TEMPLATE_ARGS
#undef COUNTER_TEMPLATE
#undef COUNTERS_TYPE
#undef COUNTERS_TEMPLATE_ARGS
#undef COUNTERS_TEMPLATE

#endif 