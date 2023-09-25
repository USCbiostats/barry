#ifndef BARRY_POWERSET_MEAT_HPP
#define BARRY_POWERSET_MEAT_HPP 1

template <typename Array_Type, typename Data_Rule_Type>
inline PowerSet<Array_Type,Data_Rule_Type>::PowerSet(
    const Array_Type & array
) : EmptyArray(array), data(0u),
        rules(new Rules<Array_Type,Data_Rule_Type>()), N(array.nrow()), M(array.ncol()) {

}

template <typename Array_Type, typename Data_Rule_Type>
inline PowerSet<Array_Type,Data_Rule_Type>::~PowerSet() {
    if (!this->rules_deleted)
        delete rules;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::init_support()
{
    
    // Computing the locations
    coordinates_free.clear();
    coordinates_locked.clear();
    rules->get_seq(EmptyArray, &coordinates_free, &coordinates_locked);

    n_free   = coordinates_free.size() / 2u;
    n_locked = coordinates_locked.size() / 2u;
    
    // Computing initial statistics
    if (EmptyArray.nnozero() > 0u)
    {

        if (EmptyArray.is_dense())
        {

            for (size_t i = 0u; i < n_free; ++i) 
                EmptyArray(
                    coordinates_free[i * 2u],
                    coordinates_free[i * 2u + 1u]
                    ) = 0;

        }
        else
        {

            for (size_t i = 0u; i < n_free; ++i) 
                EmptyArray.rm_cell(
                    coordinates_free[i * 2u],
                    coordinates_free[i * 2u + 1u],
                    false,
                    true
                );


        }
            
    }

    // EmptyArray.clear(true);
    // EmptyArray.reserve();
    
    // Resizing support
    data.reserve(pow(2.0, n_free)); 

    // Adding the empty array to the set
    data.push_back(EmptyArray);
    
    return;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type, Data_Rule_Type>::calc_backend_sparse(
    size_t pos
)
{
    
    // Did we reached the end??
    if (pos >= n_free)
        return;
            
    // We will pass it to the next step, if the iteration makes sense.
    calc_backend_sparse(pos + 1u);
        
    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray.insert_cell(
        coordinates_free[pos * 2u],
        coordinates_free[pos * 2u + 1u],
        EmptyArray.default_val().value,
        false, false
        );

    data.push_back(EmptyArray);

    #ifdef BARRY_USER_INTERRUPT
    if (data.size() % 1000u == 0u)
    {
        #BARRY_USER_INTERRUPT
    }
    #endif
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc_backend_sparse(pos + 1u);
    
    // We need to restore the state of the cell
    EmptyArray.rm_cell(
        coordinates_free[pos * 2u],
        coordinates_free[pos * 2u + 1u],
        false, false
        );  
    
    return;
    
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type, Data_Rule_Type>::calc_backend_dense(
    size_t pos
)
{
    
    // Did we reached the end??
    if (pos >= n_free)
        return;
            
    // We will pass it to the next step, if the iteration makes sense.
    calc_backend_dense(pos + 1u);
        
    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray(coordinates_free[pos * 2u], coordinates_free[pos * 2u + 1u]) = 1;

    data.push_back(EmptyArray);
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc_backend_dense(pos + 1u);
    
    // We need to restore the state of the cell
    EmptyArray(coordinates_free[pos * 2u], coordinates_free[pos * 2u + 1u]) = 0;
    
    return;
    
}


/***
  * Function to generate the powerset of the 
  */
template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type, Data_Rule_Type>::calc() {

    // Generating sequence
    this->init_support();

    // Recursive function to count
    if (EmptyArray.is_dense())
        calc_backend_dense(0u);
    else
        calc_backend_sparse(0u);

    return;
    
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::reset(
        size_t N_,
        size_t M_
) {
    
    data.empty();
    N = N_, M = M_;
    
    return;

}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
        Rule<Array_Type, Data_Rule_Type> rule
) {
    
    rules->add_rule(rule);
    return;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
        Rule_fun_type<Array_Type,Data_Rule_Type> rule_fun_,
        Data_Rule_Type data_
) {
    
    rules->add_rule(
        rule_fun_,
        data_
    );
    
    return;
    
}

#endif