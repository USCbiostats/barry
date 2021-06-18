#include "statscounter-bones.hpp"

#ifndef BARRY_STATSCOUNTER_MEAT_HPP
#define BARRY_STATSCOUNTER_MEAT_HPP 1

template <typename Array_Type, typename Data_Type>
inline StatsCounter<Array_Type,Data_Type>::~StatsCounter() {
    if (!counter_deleted)
        delete counters;
    return;
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type,Data_Type>::reset_array(
    const Array_Type * Array_
) {
    
    Array = Array_;
    EmptyArray = *Array_;
    
    return;
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type,Data_Type>::add_counter(
        Counter<Array_Type,Data_Type> * f_
    ) {
    
    counters->add_counter(f_);
    return;
    
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type,Data_Type>::add_counter(
        Counter<Array_Type,Data_Type> f_
) {
    
    counters->add_counter(f_);
    
    return;
    
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type,Data_Type>::set_counters(
        Counters<Array_Type,Data_Type> * counters_
) {
    
    // Cleaning up before replacing the memory
    if (!counter_deleted)
        delete counters;
    counter_deleted = true;
    counters = counters_;
    
    return;
    
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type, Data_Type>::count_init(
        uint i,
        uint j
    ) {
    
    // Do we have any counter?
    if (counters->size() == 0u)
        throw std::logic_error("No counters added: Cannot count without knowning what to count!");
    
    // Iterating through the functions, and updating the set of
    // statistics.
    current_stats.resize(counters->size(), 0.0);
    // change_stats.resize(counters->size(), 0.0);
    for (uint n = 0u; n < counters->size(); ++n) 
        current_stats[n] = counters->operator[](n).init(EmptyArray, i, j);
    
    return;
}

template <typename Array_Type, typename Data_Type>
inline void StatsCounter<Array_Type, Data_Type>::count_current(
        uint i,
        uint j
    ) {
    
    // Iterating through the functions, and updating the set of
    // statistics.
    for (uint n = 0u; n < counters->size(); ++n) {
        // change_stats[n]   = counters->operator[](n).count(EmptyArray, i, j);
        // current_stats[n] += change_stats[n];
        current_stats[n] += counters->operator[](n).count(EmptyArray, i, j);
    }

    return;
    
}

template <typename Array_Type, typename Data_Type>
inline std::vector< double > StatsCounter<Array_Type, Data_Type>::count_all() {
    
    // Initializing the counter on the empty array
    count_init(0u, 0u);
    
    // Setting it to zero.
    EmptyArray.clear(false);
    
    // Start iterating through the data
    for (uint i = 0; i < Array->nrow(); ++i) {
        
        const auto & row = Array->row(i, false);

        // Any element?
        if (row.size() == 0u)
            continue;
        
        // If there's one, then update the statistic, by iterating
        for (auto& col: row) {

            // We only insert if it is different from zero
            if (static_cast<int>(col.second.value) == 0)
                continue;
            
            // Adding a cell
            EmptyArray.insert_cell(i, col.first, col.second, false, false);

            // Computing the change statistics
            count_current(i, col.first);
          
        } 
        
    }
    
    // Adding to the sufficient statistics
    return current_stats;
    
}

template <typename Array_Type, typename Data_Type>
inline Counters<Array_Type,Data_Type> * StatsCounter<Array_Type, Data_Type>::get_counters() {
    return this->counters;
}


#endif