#include "statscounter-bones.hpp"

#ifndef BARRY_STATSCOUNTER_MEAT_HPP
#define BARRY_STATSCOUNTER_MEAT_HPP 1

#define STATSCOUNTER_TYPE() StatsCounter<Array_Type,Data_Type>

#define STATSCOUNTER_TEMPLATE_ARGS() <typename Array_Type, typename Data_Type>

#define STATSCOUNTER_TEMPLATE(a,b) \
    template STATSCOUNTER_TEMPLATE_ARGS() inline a STATSCOUNTER_TYPE()::b

STATSCOUNTER_TEMPLATE(,~StatsCounter)()
{
    if (!counter_deleted)
        delete counters;
    return;
}

STATSCOUNTER_TEMPLATE(void, reset_array)(const Array_Type * Array_)
{
    
    Array = Array_;
    EmptyArray = *Array_;
    
    return;
}

STATSCOUNTER_TEMPLATE(void, add_counter)(Counter<Array_Type,Data_Type> * f_)
{
    
    counters->add_counter(f_);
    return;
    
}

STATSCOUNTER_TEMPLATE(void, add_counter)(Counter<Array_Type,Data_Type> f_)
{
    
    counters->add_counter(f_);
    
    return;
    
}

STATSCOUNTER_TEMPLATE(void, set_counters)(Counters<Array_Type,Data_Type> * counters_)
{
    
    // Cleaning up before replacing the memory
    if (!counter_deleted)
        delete counters;
    counter_deleted = true;
    counters = counters_;
    
    return;
    
}

STATSCOUNTER_TEMPLATE(void, count_init)(uint i,uint j)
{
    
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

STATSCOUNTER_TEMPLATE(void, count_current)(uint i, uint j)
{
    
    // Iterating through the functions, and updating the set of
    // statistics.
    for (uint n = 0u; n < counters->size(); ++n) {
        // change_stats[n]   = counters->operator[](n).count(EmptyArray, i, j);
        // current_stats[n] += change_stats[n];
        current_stats[n] += counters->operator[](n).count(EmptyArray, i, j);
    }

    return;
    
}

STATSCOUNTER_TEMPLATE(std::vector< double >, count_all)()
{
    
    // Initializing the counter on the empty array
    count_init(0u, 0u);
    
    // Setting it to zero.
    EmptyArray.clear(false);

    #ifdef BARRY_DEBUG_LEVEL
        #if BARRY_DEBUG_LEVEL > 0
            BARRY_DEBUG_MSG("Initializing -count_all- debug. get_names():")
            BARRY_DEBUG_VEC_PRINT<std::string>(this->get_names());
        #endif
    #endif
    
    // Start iterating through the data
    for (uint i = 0; i < Array->nrow(); ++i)
    {
        
        const auto & row = Array->row(i, false);

        // Any element?
        if (row.size() == 0u)
            continue;
        
        // If there's one, then update the statistic, by iterating
        for (auto& col: row)
        {

            // We only insert if it is different from zero
            if (static_cast<int>(col.second.value) == 0)
                continue;
            
            // Adding a cell
            EmptyArray.insert_cell(i, col.first, col.second, false, false);

            #ifdef BARRY_DEBUG_LEVEL
                #if (BARRY_DEBUG_LEVEL >= 1)
                    BARRY_DEBUG_MSG("================================================================================")
                    BARRY_DEBUG_MSG("Debugging Stats counter: current_stats (before)")
                    std::string tmpmgs = "Inserting cell (" +
                        std::to_string(i) + ", " + std::to_string(col.first) + ")";
                    BARRY_DEBUG_MSG(tmpmgs.c_str());
                    BARRY_DEBUG_VEC_PRINT(current_stats);
                    #if (BARRY_DEBUG_LEVEL >= 2)
                        BARRY_DEBUG_MSG("Debugging Stats counter: EmptyArray")
                        EmptyArray.print();
                    #endif
                #endif
            #endif 

            // Computing the change statistics
            count_current(i, col.first);
            #ifdef BARRY_DEBUG_LEVEL
                #if (BARRY_DEBUG_LEVEL >= 1)
                    BARRY_DEBUG_MSG("Debugging Stats counter: current_stats (after)")
                    BARRY_DEBUG_VEC_PRINT(current_stats);
                #endif
            #endif
          
        } 
        
    }
    
    // Adding to the sufficient statistics
    return current_stats;
    
}

template STATSCOUNTER_TEMPLATE_ARGS()
inline Counters<Array_Type,Data_Type> * STATSCOUNTER_TYPE()::get_counters() {
    return this->counters;
}

STATSCOUNTER_TEMPLATE(std::vector< std::string >, get_names)() const
{
    return this->counters->get_names();
}

STATSCOUNTER_TEMPLATE(std::vector< std::string >, get_descriptions)() const
{
    return this->counters->get_descriptions();
}

#undef STATSCOUNTER_TYPE
#undef STATSCOUNTER_TEMPLATE_ARGS
#undef STATSCOUNTER_TEMPLATE

#endif