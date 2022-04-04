#ifndef BARRY_SUPPORT_MEAT
#define BARRY_SUPPORT_MEAT_HPP 1

#define SUPPORT_TEMPLATE_ARGS() <typename Array_Type, typename \
    Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>

#define SUPPORT_TYPE() Support<Array_Type,Data_Counter_Type,Data_Rule_Type,\
    Data_Rule_Dyn_Type>

#define SUPPORT_TEMPLATE(a,b) template SUPPORT_TEMPLATE_ARGS() \
    inline a SUPPORT_TYPE()::b

SUPPORT_TEMPLATE(void, init_support)(
    std::vector< Array_Type > * array_bank,
    std::vector< double > * stats_bank
) {
    
    // Computing the locations
    coordinates_free.clear();
    coordinates_locked.clear();
    rules->get_seq(EmptyArray, &coordinates_free, &coordinates_locked);

    coordiantes_n_free   = coordinates_free.size() / 2u;
    coordiantes_n_locked = coordinates_locked.size() / 2u;
    n_counters           = counters->size();

    hashes.resize(coordiantes_n_free, 0u);
    hashes_initialized.resize(coordiantes_n_free, false);
    
    // Computing initial statistics
    if (EmptyArray.nnozero() > 0u)
    {

        for (uint i = 0u; i < coordiantes_n_free; ++i)
            EmptyArray.rm_cell(
                coordinates_free[i * 2u],
                coordinates_free[i * 2u + 1u],
                false, true
                );
                
    }

    // Looked coordinates should still be removed if these are
    // equivalent to zero
    for (unsigned int i = 0u; i < coordiantes_n_locked; ++i)
    {

        if (static_cast<int>(EmptyArray(
            coordinates_locked[i * 2u], coordinates_locked[i * 2u + 1u]
            )) == 0)

            EmptyArray.rm_cell(
                coordinates_locked[i * 2u],
                coordinates_locked[i * 2u + 1u],
                false, true
                );

    }

    // Do we have any counter?
    if (n_counters == 0u)
        throw std::logic_error("No counters added: Cannot compute the support without knowning what to count!");

    // Initial count (including constrains)
    if (coordiantes_n_locked)
    {

        StatsCounter<Array_Type,Data_Counter_Type> tmpcount(&EmptyArray);
        tmpcount.set_counters(counters);
        current_stats = tmpcount.count_all();

    }
    else
    {

        current_stats.resize(n_counters, 0.0);

        // Initialize counters
        for (uint n = 0u; n < n_counters; ++n)
        {

            current_stats[n] = counters->operator[](n).init(
                EmptyArray,
                coordinates_free[0u],
                coordinates_free[1u]
                );

        }

    }

    // Resizing support
    data.reserve(
        pow(2.0, static_cast<double>(coordiantes_n_free)),
        counters->size()
        ); 

    // Adding to the overall count
    bool include_it = rules_dyn->operator()(EmptyArray, 0u, 0u);
    if (include_it)
        data.add(current_stats, nullptr);

    change_stats.resize(coordiantes_n_free * n_counters, 0.0);
        
    if (include_it && (array_bank != nullptr)) 
        array_bank->push_back(EmptyArray);
    
    if (include_it && (stats_bank != nullptr))
        std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

    return;
}

SUPPORT_TEMPLATE(void, reset_array)() {
    
    data.clear();
    
}

SUPPORT_TEMPLATE(void, reset_array)(const Array_Type & Array_) {
    
    data.clear();
    EmptyArray = Array_;
    N = Array_.nrow();
    M = Array_.ncol();
    
}

SUPPORT_TEMPLATE(void, calc_backend_sparse)(
        uint pos,
        std::vector< Array_Type > * array_bank,
        std::vector< double > * stats_bank
    ) {
    
    // Did we reached the end??
    if (pos >= coordiantes_n_free)
        return;
            
    // We will pass it to the next step, if the iteration makes sense.
    calc_backend_sparse(pos + 1u, array_bank, stats_bank);
    
    // Once we have returned, everything will be back as it used to be, so we
    // treat the data as if nothing has changed.
    const size_t & coord_i = coordinates_free[pos * 2u];
    const size_t & coord_j = coordinates_free[pos * 2u + 1u];

    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray.insert_cell(
        coord_i,
        coord_j,
        EmptyArray.default_val().value,
        false, false
        );

    // Counting
    // std::vector< double > change_stats(counters.size());
    double tmp_chng;
    unsigned int change_stats_different = hashes_initialized[pos] ? 0u : 1u;
    for (uint n = 0u; n < n_counters; ++n)
    {

        tmp_chng = counters->operator[](n).count(
            EmptyArray,
            coord_i,
            coord_j
            );
        
        if ((tmp_chng < DBL_MIN) & (tmp_chng > -DBL_MIN))
        {

            change_stats[pos * n_counters + n] = 0.0;

        }
        else
        {

            change_stats_different++;
            current_stats[n] += tmp_chng;
            change_stats[pos * n_counters + n] = tmp_chng;

        }

    }
    
    // Adding to the overall count
    BARRY_CHECK_SUPPORT(data, max_num_elements)
    if (rules_dyn->size() > 0u)
    {
        
        if (rules_dyn->operator()(
            EmptyArray,
            coord_i,
            coord_j
            ))
        {

            if (change_stats_different > 0u)
                hashes[pos] = data.add(current_stats, nullptr);
            else
                (void) data.add(current_stats, &hashes[pos]);

            // Need to save?
            if (array_bank != nullptr)
                array_bank->push_back(EmptyArray);
            
            if (stats_bank != nullptr)
                std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

        }
            

    } else {

        if (change_stats_different > 0u)
            hashes[pos] = data.add(current_stats, nullptr);
        else
            (void) data.add(current_stats, &hashes[pos]);

        // Need to save?
        if (array_bank != nullptr)
            array_bank->push_back(EmptyArray);
        
        if (stats_bank != nullptr)
            std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

    }
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc_backend_sparse(pos + 1u, array_bank, stats_bank);
    
    // We need to restore the state of the cell
    EmptyArray.rm_cell(
        coord_i,
        coord_j,
        false, false
        );
    
    if (change_stats_different > 0u)
    {
        #ifdef __OPENMP
        #pragma omp simd
        #else
        #pragma GCC ivdep
        #endif
        for (uint n = 0u; n < n_counters; ++n) 
            current_stats[n] -= change_stats[pos * n_counters + n];
    }
        
    
    return;
    
}

SUPPORT_TEMPLATE(void, calc_backend_dense)(
        uint pos,
        std::vector< Array_Type > * array_bank,
        std::vector< double > * stats_bank
    ) {
    
    // Did we reached the end??
    if (pos >= coordiantes_n_free)
        return;
            
    // We will pass it to the next step, if the iteration makes sense.
    calc_backend_dense(pos + 1u, array_bank, stats_bank);
    
    // Once we have returned, everything will be back as it used to be, so we
    // treat the data as if nothing has changed.
    const size_t & coord_i = coordinates_free[pos * 2u];
    const size_t & coord_j = coordinates_free[pos * 2u + 1u];

    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray.insert_cell(coord_i, coord_j, 1, false, false);

    // Counting
    // std::vector< double > change_stats(counters.size());
    double tmp_chng;
    unsigned int change_stats_different = hashes_initialized[pos] ? 0u : 1u;
    for (uint n = 0u; n < n_counters; ++n)
    {

        tmp_chng = counters->operator[](n).count(
            EmptyArray,
            coord_i,
            coord_j
            );

        if ((tmp_chng < DBL_MIN) & (tmp_chng > -DBL_MIN))
        {

            change_stats[pos * n_counters + n] = 0.0;

        }
        else
        {
            if (std::isnan(tmp_chng))
                throw std::domain_error("Undefined number.");

            change_stats_different++;
            current_stats[n] += tmp_chng;
            change_stats[pos * n_counters + n] = tmp_chng;

        }

    }
    
    // Adding to the overall count
    BARRY_CHECK_SUPPORT(data, max_num_elements)
    if (rules_dyn->size() > 0u)
    {
        
        if (rules_dyn->operator()(EmptyArray, coord_i, coord_j))
        {

            if (change_stats_different > 0u)
                hashes[pos] = data.add(current_stats, nullptr);
            else
                (void) data.add(current_stats, &hashes[pos]);

            // Need to save?
            if (array_bank != nullptr)
                array_bank->push_back(EmptyArray);
            
            if (stats_bank != nullptr)
                std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

        }
            

    }
    else
    {

        if (change_stats_different > 0u)
            hashes[pos] = data.add(current_stats, nullptr);
        else
            (void) data.add(current_stats, &hashes[pos]);

        // Need to save?
        if (array_bank != nullptr)
            array_bank->push_back(EmptyArray);
        
        if (stats_bank != nullptr)
            std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

    }
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc_backend_dense(pos + 1u, array_bank, stats_bank);
    
    // We need to restore the state of the cell
    EmptyArray.rm_cell(coord_i, coord_j, false, false);
    
    if (change_stats_different > 0u)
    {
        #ifdef __OPENMP
        #pragma omp simd
        #else
        #pragma GCC ivdep
        #endif
        for (uint n = 0u; n < n_counters; ++n) 
            current_stats[n] -= change_stats[pos * n_counters + n];
    }
    
    return;
    
}

SUPPORT_TEMPLATE(void, calc)(
        std::vector< Array_Type > * array_bank,
        std::vector< double > * stats_bank,
        unsigned int max_num_elements_
) {

    if (max_num_elements_ != 0u)
        this->max_num_elements = max_num_elements_;

    // Generating sequence
    this->init_support(array_bank, stats_bank);

    // Recursive function to count
    if (EmptyArray.is_dense())
        calc_backend_dense(0u, array_bank, stats_bank);
    else
        calc_backend_sparse(0u, array_bank, stats_bank);

    change_stats.clear();

    if (max_num_elements_ != 0u)
        this->max_num_elements = BARRY_MAX_NUM_ELEMENTS;


    return;
    
}

SUPPORT_TEMPLATE(void, add_counter)(
        Counter<Array_Type,Data_Counter_Type> f_
) {
    
    counters->add_counter(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, set_counters)(
        Counters<Array_Type,Data_Counter_Type> * counters_
) {
    
    // Cleaning up before replacing the memory
    if (delete_counters)
        delete counters;
    delete_counters = false;
    counters = counters_;
    
    return;
    
}

/////////////////////////////

SUPPORT_TEMPLATE(void, add_rule)(
        Rule<Array_Type, Data_Rule_Type> * f_
) {
    
    rules->add_rule(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, add_rule)(
        Rule<Array_Type,Data_Rule_Type> f_
) {
    
    rules->add_rule(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, set_rules)(
        Rules<Array_Type,Data_Rule_Type> * rules_
) {
    
    // Cleaning up before replacing the memory
    if (delete_rules)
        delete rules;
    delete_rules = false;
    rules = rules_;
    
    return;
    
}

SUPPORT_TEMPLATE(void, add_rule_dyn)(
        Rule<Array_Type, Data_Rule_Dyn_Type> * f_
) {
    
    rules_dyn->add_rule(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, add_rule_dyn)(
        Rule<Array_Type,Data_Rule_Dyn_Type> f_
) {
    
    rules_dyn->add_rule(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, set_rules_dyn)(
        Rules<Array_Type,Data_Rule_Dyn_Type> * rules_
) {
    
    // Cleaning up before replacing the memory
    if (delete_rules_dyn)
        delete rules_dyn;
    delete_rules_dyn = false;
    rules_dyn = rules_;
    
    return;
    
}

SUPPORT_TEMPLATE(bool, eval_rules_dyn)(
    const std::vector< double > & counts,
    const uint & i,
    const uint & j
) {

    if (rules_dyn->size() == 0u)
        return true;

    // Swapping pointers for a while
    std::vector< double > tmpstats = current_stats;
    current_stats = counts;

    bool rule_res = rules_dyn->operator()(EmptyArray, i, j);
    current_stats = tmpstats;

    return rule_res;

}

// SUPPORT_TEMPLATE(bool, eval_rules_dyn)(
//     const double * counts,
//     const uint & i,
//     const uint & j
// ) {

//     if (rules_dyn->size() == 0u)
//         return true;

//     // Swapping pointers for a while
//     std::vector< double > tmpstats = current_stats;
//     current_stats = counts;

//     bool rule_res = rules_dyn->operator()(EmptyArray, i, j);
//     current_stats = tmpstats;

//     return rule_res;

// }

//////////////////////////

SUPPORT_TEMPLATE(std::vector< double >, get_counts)() const {
    
    return data.get_data(); 
    
}

// SUPPORT_TEMPLATE(const MapVec_type<> *, get_counts_ptr)() const {
    
//     return data.get_data_ptr();
      
// }

SUPPORT_TEMPLATE(std::vector< double > *, get_current_stats)() {
    return &this->current_stats;
}

SUPPORT_TEMPLATE(void, print)() const {

    // Starting from the name of the stats
    printf_barry("Position of variables:\n");
    for (uint i = 0u; i < n_counters; ++i) {
        printf_barry("[% 2i] %s\n", i, counters->operator[](i).name.c_str());
    }

    data.print();
}

SUPPORT_TEMPLATE(const FreqTable<double> &, get_data)() const {
    return this->data;
}

template SUPPORT_TEMPLATE_ARGS()
inline Counters<Array_Type,Data_Counter_Type> * SUPPORT_TYPE()::get_counters() {
    return this->counters;
}   
    
template SUPPORT_TEMPLATE_ARGS()
inline Rules<Array_Type,Data_Rule_Type> * SUPPORT_TYPE()::get_rules() {
    return this->rules;
}

template SUPPORT_TEMPLATE_ARGS()
inline Rules<Array_Type,Data_Rule_Dyn_Type> * SUPPORT_TYPE()::get_rules_dyn() {
    return this->rules_dyn;
}

#undef SUPPORT_TEMPLATE_ARGS
#undef SUPPORT_TYPE
#undef SUPPORT_TEMPLATE

#endif