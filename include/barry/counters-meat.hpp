#ifndef BARRY_COUNTERS_MEAT_HPP
#define BARRY_COUNTERS_MEAT_HPP 1

#define COUNTER_TYPE() Counter<Array_Type,Data_Type>

#define COUNTER_TEMPLATE_ARGS() <typename Array_Type, typename Data_Type>

#define COUNTER_TEMPLATE(a,b) \
    template COUNTER_TEMPLATE_ARGS() inline a COUNTER_TYPE()::b

COUNTER_TEMPLATE(,Counter)(
    const Counter<Array_Type,Data_Type> & counter_
) : count_fun(counter_.count_fun), init_fun(counter_.init_fun), hasher_fun(counter_.hasher_fun) {

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
    hasher_fun(std::move(counter_.hasher_fun)),
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
        this->hasher_fun = counter_.hasher_fun;

        
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
        this->hasher_fun = std::move(counter_.hasher_fun);

        // Descriptions
        this->name = std::move(counter_.name);
        this->desc = std::move(counter_.desc);

    }

    return *this;

} ///< Move assignment

COUNTER_TEMPLATE(double, count)(Array_Type & Array, size_t i, size_t j)
{

    if (count_fun == nullptr)
        return 0.0;

    return count_fun(Array, i, j, data);

}

COUNTER_TEMPLATE(double, init)(Array_Type & Array, size_t i, size_t j)
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

COUNTER_TEMPLATE(void, set_name)(std::string new_name) {
    name = new_name;
}

COUNTER_TEMPLATE(void, set_description)(std::string new_desc) {
    desc = new_desc;
}

COUNTER_TEMPLATE(void, set_hasher)(Hasher_fun_type<Array_Type,Data_Type> fun) {
    hasher_fun = fun;
}

#define TMP_HASHER_CALL Hasher_fun_type<Array_Type,Data_Type>
COUNTER_TEMPLATE(TMP_HASHER_CALL, get_hasher)() {
    return hasher_fun;
}
#undef TMP_HASHER_CALL

////////////////////////////////////////////////////////////////////////////////
// Counters
////////////////////////////////////////////////////////////////////////////////

#define COUNTERS_TYPE() Counters<Array_Type,Data_Type>

#define COUNTERS_TEMPLATE_ARGS() <typename Array_Type, typename Data_Type>

#define COUNTERS_TEMPLATE(a,b) \
    template COUNTERS_TEMPLATE_ARGS() inline a COUNTERS_TYPE()::b

COUNTERS_TEMPLATE(, Counters)() : data(0u), hasher(nullptr) {}

COUNTERS_TEMPLATE(COUNTER_TYPE() &, operator[])(size_t idx) {

    return data[idx];

}

COUNTERS_TEMPLATE(, Counters)(const Counters<Array_Type,Data_Type> & counter_) :
    data(counter_.data), hasher(counter_.hasher) {}

COUNTERS_TEMPLATE(, Counters)(Counters<Array_Type,Data_Type> && counters_) noexcept :
    data(std::move(counters_.data)), hasher(std::move(counters_.hasher)) {}

COUNTERS_TEMPLATE(COUNTERS_TYPE(), operator=)(const Counters<Array_Type,Data_Type> & counter_) {

    if (this != &counter_)
    {
        data = counter_.data;
        hasher = counter_.hasher;
    }

    return *this;

}

COUNTERS_TEMPLATE(COUNTERS_TYPE() &, operator=)(Counters<Array_Type,Data_Type> && counters_) noexcept 
{

    if (this != &counters_) {
        data = std::move(counters_.data);
        hasher = std::move(counters_.hasher);
    }

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
    Hasher_fun_type<Array_Type,Data_Type>  hasher_fun_,
    Data_Type                              data_,
    std::string                            name_,
    std::string                            desc_
)
{
  
    data.emplace_back(Counter<Array_Type,Data_Type>(
        count_fun_,
        init_fun_,
        hasher_fun_,
        data_,
        name_,
        desc_
    ));
  
    return;
    
}

COUNTERS_TEMPLATE(std::vector<std::string>, get_names)() const
{

    std::vector< std::string > out;
    out.reserve(this->size());
    for (size_t i = 0u; i < this->size(); ++i)
        out.push_back(this->data.at(i).get_name());

    return out;

}

COUNTERS_TEMPLATE(std::vector<std::string>, get_descriptions)() const
{
    
    std::vector< std::string > out;
    out.reserve(this->size());
    for (size_t i = 0u; i < this->size(); ++i)
        out.push_back(data.at(i).get_description());

    return out;

}

COUNTERS_TEMPLATE(std::vector<double>, gen_hash)(
    const Array_Type & array,
    bool add_dims
)
{
    std::vector<double> res;
    
    // Iterating over the counters
    for (auto & c: data)
    {

        // If there's a hasher function, then use it!
        if (c.get_hasher())
        {

            for (auto v: c.get_hasher()(array, &(c.data)))
                res.push_back(v);

        }

    }

    // Do we need to add the dims?
    if (add_dims)
    {
        res.push_back(array.nrow());
        res.push_back(array.ncol());
    }

    // Ading the global hasher, if one exists
    if (hasher)
    {
        for (auto i: hasher(array, nullptr))
            res.push_back(i);
    }

    // We have to return something...
    if (res.size() == 0u)
        res.push_back(0.0);

    return res;

}

COUNTERS_TEMPLATE(void, add_hash)(
    Hasher_fun_type<Array_Type,Data_Type> fun_
) {

    hasher = fun_;

}

#undef COUNTER_TYPE
#undef COUNTER_TEMPLATE_ARGS
#undef COUNTER_TEMPLATE
#undef COUNTERS_TYPE
#undef COUNTERS_TEMPLATE_ARGS
#undef COUNTERS_TEMPLATE

#endif 