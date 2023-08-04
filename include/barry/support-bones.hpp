#ifndef BARRY_SUPPORT_BONES_HPP 
#define BARRY_SUPPORT_BONES_HPP 1

template <typename Cell_Type, typename Data_Type>
class BArray;

template <typename Tdat>
class FreqTable;

template <typename Array_Type, typename Data_Counter_Type>
class Counters;

template <typename Array_Type, typename Data_Rule_Type>
class Rules;

template<typename Array_Type, typename Data_Type>
class Rule;

/**
 * @brief Compute the support of sufficient statistics
 * 
 * Given an array and a set of counters, this object iterates throughout the
 * support set of the Array while at the same time computing the support of
 * the sufficient statitics.
 * 
 * The members `rule` and `rule_dyn` allow constraining the support. The first
 * will establish which cells of the array will be used to iterate, for example,
 * in the case of social networks, self-loops are not allowed, so the entire
 * diagonal would be fixed to zero, reducing the size of the support.
 * 
 * In the case of `rule_dyn`, the function will stablish dynamically whether
 * the current state will be included in the counts or not. For example, this
 * set of rules can be used to constrain the support to networks that have a
 * prescribed degree sequence. 
 */ 
template <
    typename Array_Type         = BArray<bool, bool>,
    typename Data_Counter_Type  = bool,
    typename Data_Rule_Type     = bool,
    typename Data_Rule_Dyn_Type = bool 
    >
class Support {
    
private:
    void calc_backend_sparse(
        size_t pos = 0u,
        std::vector< Array_Type > * array_bank = nullptr,
        std::vector< double > * stats_bank = nullptr
    );

    void calc_backend_dense(
        size_t pos = 0u,
        std::vector< Array_Type > * array_bank = nullptr,
        std::vector< double > * stats_bank = nullptr
    );

    /**
     * @brief Reference array to generate the support.
     */
    Array_Type                               EmptyArray; ///< Temp array used to iterate through the support.
    FreqTable<>                              data;       ///< Table with the support.
    Counters<Array_Type,Data_Counter_Type> * counters;   ///< Vector of couter functions.
    Rules<Array_Type,Data_Rule_Type> *       rules;      ///< Vector of static rules (cells to iterate).
    Rules<Array_Type,Data_Rule_Dyn_Type> *   rules_dyn;  ///< Vector of dynamic rules (to include/exclude a realizaton).
    
public:
    
    size_t N, M;
    bool delete_counters  = true;
    bool delete_rules     = true;
    bool delete_rules_dyn = true;
    size_t max_num_elements = BARRY_MAX_NUM_ELEMENTS;
    
    // Temp variables to reduce memory allocation
    std::vector< double >                current_stats;
    std::vector< size_t >                coordinates_free;
    std::vector< size_t >                coordinates_locked;
    size_t coordiantes_n_free;
    size_t coordiantes_n_locked;
    std::vector< double > change_stats;
    std::vector< size_t > hashes;
    std::vector< bool   > hashes_initialized;
    size_t n_counters;
    
    /**@brief Constructor passing a reference Array.
      */
    Support(const Array_Type & Array_) :
        EmptyArray(Array_),
        counters(new Counters<Array_Type,Data_Counter_Type>()),
        rules(new Rules<Array_Type,Data_Rule_Type>()),
        rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
        N(Array_.nrow()), M(Array_.ncol()), current_stats() {};
    
    /**@brief Constructor specifying the dimensions of the array (empty).
      */
    Support(size_t N_, size_t M_) :
        EmptyArray(N_, M_),
        counters(new Counters<Array_Type,Data_Counter_Type>()),
        rules(new Rules<Array_Type,Data_Rule_Type>()),
        rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
        N(N_), M(M_), current_stats() {};
    
    Support() :
        EmptyArray(0u, 0u),
        counters(new Counters<Array_Type,Data_Counter_Type>()),
        rules(new Rules<Array_Type,Data_Rule_Type>()),
        rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
        N(0u), M(0u), current_stats() {};
    
    ~Support() {
        
        if (delete_counters)
            delete counters;
        if (delete_rules)
            delete rules;
        if (delete_rules_dyn)
            delete rules_dyn;

    };
    
    void init_support(
        std::vector< Array_Type > * array_bank = nullptr,
        std::vector< double > * stats_bank = nullptr
    );
    
    /**
     * @name Resets the support calculator
     * 
     * If needed, the counters of a support object can be reused.
     * 
     * @param Array_ New array over which the support will be computed.
     */
    ///@{
    void reset_array();
    void reset_array(const Array_Type & Array_);
    ///@}
    
    /**
     * @name Manage counters 
     * 
     * @param f_ A counter to be added.
     * @param counters_ A vector of counters to be added.
     */
    ///@{
    void add_counter(Counter<Array_Type,Data_Counter_Type> f_);
    void set_counters(Counters<Array_Type,Data_Counter_Type> * counters_);
    ///@}
    
    /**
     * @name Manage rules 
     * 
     * @param f_ A rule to be added.
     * @param counters_ A vector of rules to be added.
     */
    void add_rule(Rule<Array_Type, Data_Rule_Type> * f_);
    void add_rule(Rule<Array_Type,Data_Rule_Type> f_);
    void set_rules(Rules<Array_Type,Data_Rule_Type> * rules_);
    void add_rule_dyn(Rule<Array_Type, Data_Rule_Dyn_Type> * f_);
    void add_rule_dyn(Rule<Array_Type,Data_Rule_Dyn_Type> f_);
    void set_rules_dyn(Rules<Array_Type,Data_Rule_Dyn_Type> * rules_);
    bool eval_rules_dyn(const std::vector<double> & counts, const size_t & i, const size_t & j);
    // bool eval_rules_dyn(const double * counts, const size_t & i, const size_t & j);
    ///@}

    /**
     * @brief Computes the entire support
     * 
     * Not to be used by the user. Sets the starting point in the array
     * (column-major).
     *  
     * @param array_bank If specified, the counter will add to the vector each 
     * possible state of the array, as it counts.
     * 
     * @param stats_bank If specified, the counter will add to the vector each
     * possible set of statistics, as it counts.
     * 
     */
    void calc(
        std::vector< Array_Type > * array_bank = nullptr,
        std::vector< double > * stats_bank = nullptr,
        size_t max_num_elements_ = 0u
    );
    
    std::vector< double > get_counts() const;
    std::vector< double > * get_current_stats(); ///< List current statistics.
    void print() const;
    
    const FreqTable< double > &              get_data() const;
    Counters<Array_Type,Data_Counter_Type> * get_counters();   ///< Vector of couter functions.
    Rules<Array_Type,Data_Rule_Type> *       get_rules();      ///< Vector of static rules (cells to iterate).
    Rules<Array_Type,Data_Rule_Dyn_Type> *   get_rules_dyn();  ///< Vector of dynamic rules (to include/exclude a realizaton).
    
};


#endif
