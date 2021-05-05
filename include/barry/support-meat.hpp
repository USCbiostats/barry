#include "support-bones.hpp"

#ifndef BARRY_SUPPORT_MEAT
#define BARRY_SUPPORT_MEAT_HPP 1

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type>::init_support(
    std::vector< Array_Type > * array_bank,
    std::vector< std::vector< double > > * stats_bank
) {
    
    // Computing the locations
    coordinates_free.clear();
    coordinates_locked.clear();
    rules->get_seq(EmptyArray, &coordinates_free, &coordinates_locked);
    
    // Computing initial statistics
    if (EmptyArray.nnozero() > 0u) {
        for (uint i = 0u; i < coordinates_free.size(); ++i) 
            EmptyArray.rm_cell(coordinates_free[i].first, coordinates_free[i].second, false, true);
    }

    // Do we have any counter?
    if (counters->size() == 0u)
        throw std::logic_error("No counters added: Cannot compute the support without knowning what to count!");

    // Initial count (including constrains)
    if (coordinates_locked.size()) {

        StatsCounter<Array_Type,Data_Counter_Type> tmpcount(&EmptyArray);
        tmpcount.set_counters(counters);
        current_stats = tmpcount.count_all();

    } else {

        current_stats.resize(counters->size(), 0.0);

        const auto & cordfree = coordinates_free[0u];

        // Initialize counters
        for (uint n = 0u; n < counters->size(); ++n)
        {

            current_stats[n] = counters->operator[](n).init(
                EmptyArray, cordfree.first, cordfree.second
                );

        }

    }

    // Resizing support
    if (coordinates_free.size() > 9u)
        data.reserve(pow(2.0, static_cast<double>(coordinates_free.size() - 4u))); 

    // Adding to the overall count
    bool include_it = rules_dyn->operator()(EmptyArray, 0u, 0u);
    if (include_it)
        data.add(current_stats);

    change_stats.resize(coordinates_free.size(), current_stats);
        
    if (include_it && (array_bank != nullptr)) 
        array_bank->push_back(EmptyArray);
    
    if (include_it && (stats_bank != nullptr))
        stats_bank->push_back(current_stats);

    return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type>::reset_array() {
    
    data.clear();
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type>::reset_array(const Array_Type & Array_) {
    
    data.clear();
    EmptyArray = Array_;
    N = Array_.nrow();
    M = Array_.ncol();
    // init_support();
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type>::calc_backend(
        uint                                   pos,
        std::vector< Array_Type > *            array_bank,
        std::vector< std::vector< double > > * stats_bank
    ) {
    
    // Did we reached the end??
    if (pos >= coordinates_free.size())
        return;
            
    // We will pass it to the next step, if the iteration makes sense.
    calc_backend(pos + 1u, array_bank, stats_bank);
    
    // Once we have returned, everything will be back as it used to be, so we
    // treat the data as if nothing has changed.
    
    const std::pair<uint,uint> & cfree = coordinates_free[pos];

    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray.insert_cell(
        cfree.first, cfree.second,
        EmptyArray.default_val().value,
        false, false
        );

    // Counting
    // std::vector< double > change_stats(counters.size());
    for (uint n = 0u; n < counters->size(); ++n) {
        change_stats[pos][n] = counters->operator[](n).count(
            EmptyArray,
            cfree.first,
            cfree.second
            );
        current_stats[n] += change_stats[pos][n];
    }
    
    // Adding to the overall count
    BARRY_CHECK_SUPPORT(data, max_num_elements)
    if (rules_dyn->size() > 0u)
    {
        if (rules_dyn->operator()(EmptyArray, cfree.first, cfree.second))
        {

            data.add(current_stats);

            // Need to save?
            if (array_bank != nullptr)
                array_bank->push_back(EmptyArray);
            
            if (stats_bank != nullptr)
                stats_bank->push_back(current_stats);

        }
            

    } else {

        data.add(current_stats);
        // Need to save?
        if (array_bank != nullptr)
            array_bank->push_back(EmptyArray);
        
        if (stats_bank != nullptr)
            stats_bank->push_back(current_stats);

    }
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc_backend(pos + 1u, array_bank, stats_bank);
    
    // We need to restore the state of the cell
    EmptyArray.rm_cell(
        cfree.first,
        cfree.second,
        false, false
        );
    
    for (uint n = 0u; n < counters->size(); ++n) 
        current_stats[n] -= change_stats[pos][n];
    
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type>::calc(
        std::vector< Array_Type > *            array_bank,
        std::vector< std::vector< double > > * stats_bank,
        unsigned int max_num_elements_
) {

    if (max_num_elements_ != 0u)
        this->max_num_elements = max_num_elements_;

    // Generating sequence
    this->init_support(array_bank, stats_bank);

    // Recursive function to count
    calc_backend(0u, array_bank, stats_bank);

    change_stats.clear();

    if (max_num_elements_ != 0u)
        this->max_num_elements = BARRY_MAX_NUM_ELEMENTS;


    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::add_counter(
        Counter<Array_Type, Data_Counter_Type> * f_
    ) {
    
    counters->add_counter(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::add_counter(
        Counter<Array_Type,Data_Counter_Type> f_
) {
    
    counters->add_counter(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::set_counters(
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

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::add_rule(
        Rule<Array_Type, Data_Rule_Type> * f_
) {
    
    rules->add_rule(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::add_rule(
        Rule<Array_Type,Data_Rule_Type> f_
) {
    
    rules->add_rule(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::set_rules(
        Rules<Array_Type,Data_Rule_Type> * rules_
) {
    
    // Cleaning up before replacing the memory
    if (delete_rules)
        delete rules;
    delete_rules = false;
    rules = rules_;
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::add_rule_dyn(
        Rule<Array_Type, Data_Rule_Dyn_Type> * f_
) {
    
    rules_dyn->add_rule(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::add_rule_dyn(
        Rule<Array_Type,Data_Rule_Dyn_Type> f_
) {
    
    rules_dyn->add_rule(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::set_rules_dyn(
        Rules<Array_Type,Data_Rule_Dyn_Type> * rules_
) {
    
    // Cleaning up before replacing the memory
    if (delete_rules_dyn)
        delete rules_dyn;
    delete_rules_dyn = false;
    rules_dyn = rules_;
    
    return;
    
}

//////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Counts_type Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::get_counts() const {
    
    return data.as_vector(); 
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const MapVec_type<> * Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::get_counts_ptr() const {
    
    return data.get_data_ptr();
      
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< double > * Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::get_current_stats() {
    return &this->current_stats;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::print() const {

    // Starting from the name of the stats
    std::cout << "Position of variables:" << std::endl;
    for (uint i = 0u; i < counters->size(); ++i) {
        std::cout << "[" << i << "] " << counters->operator[](i).name << std::endl;
    }

    data.print();
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const FreqTable<> & Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::get_data() const {
    return this->data;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Counters<Array_Type,Data_Counter_Type> * Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::get_counters() {
    return this->counters;
}   
    
template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Rules<Array_Type,Data_Rule_Type> * Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::get_rules() {
    return this->rules;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Rules<Array_Type,Data_Rule_Dyn_Type> * Support<Array_Type,Data_Counter_Type, Data_Rule_Type,Data_Rule_Dyn_Type>::get_rules_dyn() {
    return this->rules_dyn;
}

#endif