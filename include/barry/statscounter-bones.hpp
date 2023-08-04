#ifndef BARRY_STATSCOUNTER_BONES_HPP 
#define BARRY_STATSCOUNTER_BONES_HPP 1

class NetworkDense;
class NetCounterData;

/**
 * @brief Count stats for a single Array.
 * 
 * Users can a list of functions that can be used with this. The baseline set of
 * arguments is a pointer to a binary array and a dataset to add the counts to.
 */ 
template <typename Array_Type, typename Data_Type>
class StatsCounter {

private:

    // Should receive an array
    const Array_Type *               Array;
    Array_Type                       EmptyArray;
    std::vector< double >            current_stats;
      
    // We will save the data here
    Counters<Array_Type,Data_Type> * counters;
    bool                             counter_deleted  = false;

    std::vector< double > count_all_dense();
    std::vector< double > count_all_sparse();

public:
        
    /**
     * @brief Creator of a `StatsCounter`
     * 
     * @param Array_ A const pointer to a `BArray`.
     */
    StatsCounter(const Array_Type * Array_) :
        Array(Array_), EmptyArray(*Array_),
        counters(new Counters<Array_Type,Data_Type>()) {
        
        // We are removing the entries without freeing the memory. This should
        // make the insertion faster.
        EmptyArray.clear(false);
        
        return;
    }

    /**
     * @brief Copy constructor
     * 
     * @param counter 
     */
    StatsCounter(const StatsCounter<Array_Type,Data_Type> & counter);
    
    /**
     * @brief Can be created without setting the array.
     * 
     */
    StatsCounter() : Array(nullptr), EmptyArray(0u,0u),
        counters(new Counters<Array_Type,Data_Type>()) {};
    ~StatsCounter();
    
    /**
     * @brief Changes the reference array for the counting.
     * 
     * @param Array_ A pointer to an array of class `Array_Type`.
     */
    void reset_array(const Array_Type * Array_);
    
    void add_counter(Counter<Array_Type,Data_Type> f_);
    void set_counters(Counters<Array_Type,Data_Type> * counters_);
    
    /**
     * @brief Counter functions
     * This function recurses through the entries of `Array` and at each step of
     * adding a new cell it uses the functions to list the statistics.
     */
    void count_init(size_t i, size_t j);
    void count_current(size_t i, size_t j);
    std::vector< double > count_all();

    Counters<Array_Type,Data_Type> * get_counters();
    std::vector< std::string > get_names() const;
    std::vector< std::string > get_descriptions() const;

    size_t size() const {return counters->size();};
    
};

#endif
