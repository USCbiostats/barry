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
inline Data_Type & Rule<Array_Type,Data_Type>::D()
{
    return dat;
}

template<typename Array_Type, typename Data_Type>
inline bool Rule<Array_Type,Data_Type>::operator()(const Array_Type & a, size_t i, size_t j) {
    return fun(a, i, j, dat);
}

template<typename Array_Type, typename Data_Type>
inline std::string & Rule<Array_Type,Data_Type>::get_name()
{
    return name;
}

template<typename Array_Type, typename Data_Type>
inline std::string & Rule<Array_Type,Data_Type>::get_description()
{
    return desc;
}

template<typename Array_Type, typename Data_Type>
inline std::string Rule<Array_Type,Data_Type>::get_name() const
{
    return name;
}

template<typename Array_Type, typename Data_Type>
inline std::string Rule<Array_Type,Data_Type>::get_description() const
{
    return desc;
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
        Data_Type                           data_,
        std::string name_,
        std::string description_
) {
       
    data.push_back(Rule<Array_Type,Data_Type>(
        rule_,
        data_,
        name_,
        description_
    ));
    
    return;
    
}

template <typename Array_Type, typename Data_Type>
inline bool Rules<Array_Type,Data_Type>::operator()(
    const Array_Type & a, size_t i, size_t j
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

    
    size_t N = a.nrow();
    size_t K = a.ncol();
    
    // Reserving some space
    (void) free->empty();
    (void) free->reserve(2u * N * K);
    
    for (size_t i = 0u; i < N; ++i)
    {

        for (size_t j = 0u; j < K; ++j)
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

template<typename Array_Type, typename Data_Type>
inline std::vector<std::string> Rules<Array_Type, Data_Type>::get_names() const
{

    std::vector< std::string > out(this->size());
    for (size_t i = 0u; i < out.size(); ++i)
        out[i] = this->data.at(i).get_name();

    return out;

}

template<typename Array_Type, typename Data_Type>
inline std::vector<std::string> Rules<Array_Type, Data_Type>::get_descriptions() const
{
    
    std::vector< std::string > out(this->size());
    for (size_t i = 0u; i < out.size(); ++i)
        out[i] = data.at(i).get_description();

    return out;

}

#endif
