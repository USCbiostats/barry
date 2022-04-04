#ifndef BARRY_COUNTERS_MEAT_HPP
#define BARRY_COUNTERS_MEAT_HPP 1

#define COUNTER_TYPE() Counter<Array_Type,Data_Type>

#define COUNTER_TEMPLATE_ARGS() <typename Array_Type, typename Data_Type>

#define COUNTER_TEMPLATE(a,b) \
    template COUNTER_TEMPLATE_ARGS() inline a COUNTER_TYPE()::b

COUNTER_TEMPLATE(,Counter)(
    const Counter<Array_Type,Data_Type> & counter_
) : count_fun(counter_.count_fun), init_fun(counter_.init_fun) {

    this->data = counter_.data;
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
    name(std::move(counter_.name)),
    desc(std::move(counter_.desc))
{

} ///< Move constructor

COUNTER_TEMPLATE(COUNTER_TYPE(),operator=)(
    const Counter<Array_Type,Data_Type> & counter_
)
{

    if (this != &counter_) {

        this->count_fun = counter_.count_fun;
        this->init_fun = counter_.init_fun;

        
        this->data = counter_.data;
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

        this->data = std::move(counter_.data);

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

COUNTERS_TEMPLATE(, Counters)() : data(0u) {}

COUNTERS_TEMPLATE(COUNTER_TYPE() &, operator[])(uint idx) {

    return data[idx];

}

COUNTERS_TEMPLATE(, Counters)(const Counters<Array_Type,Data_Type> & counter_) :
    data(counter_.data) {}

COUNTERS_TEMPLATE(, Counters)(Counters<Array_Type,Data_Type> && counters_) noexcept :
    data(std::move(counters_.data)) {}

COUNTERS_TEMPLATE(COUNTERS_TYPE(), operator=)(const Counters<Array_Type,Data_Type> & counter_) {

    if (this != &counter_)
        data = counter_.data;

    return *this;

}

COUNTERS_TEMPLATE(COUNTERS_TYPE() &, operator=)(Counters<Array_Type,Data_Type> && counters_) noexcept 
{

    if (this != &counters_)
        data = std::move(counters_.data);

    return *this;

}

COUNTERS_TEMPLATE(void, add_counter)(Counter<Array_Type, Data_Type> counter)
{
    
    data.push_back(counter);
    
    return;
}

COUNTERS_TEMPLATE(void, add_counter)(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_,
    Data_Type                              data_,
    std::string                            name_,
    std::string                            desc_
)
{
  
    data.push_back(Counter<Array_Type,Data_Type>(
        count_fun_,
        init_fun_,
        data_,
        name_,
        desc_
    ));
  
    return;
    
}

COUNTERS_TEMPLATE(std::vector<std::string>, get_names)() const
{

    std::vector< std::string > out(this->size());
    for (unsigned int i = 0u; i < out.size(); ++i)
        out[i] = this->data.at(i).get_name();

    return out;

}

COUNTERS_TEMPLATE(std::vector<std::string>, get_descriptions)() const
{
    
    std::vector< std::string > out(this->size());
    for (unsigned int i = 0u; i < out.size(); ++i)
        out[i] = data.at(i).get_description();

    return out;

}

#undef COUNTER_TYPE
#undef COUNTER_TEMPLATE_ARGS
#undef COUNTER_TEMPLATE
#undef COUNTERS_TYPE
#undef COUNTERS_TEMPLATE_ARGS
#undef COUNTERS_TEMPLATE

#endif 