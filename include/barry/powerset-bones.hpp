#ifndef BARRY_POWERSET_BONES_HPP 
#define BARRY_POWERSET_BONES_HPP 1

/**
 * @brief Powerset of a binary array
 * 
 * @tparam Array_Type 
 * @tparam Data_Rule_Type 
 */
template <typename Array_Type = BArray<>, typename Data_Rule_Type = bool> 
class PowerSet {
    
private:
    void calc_backend_sparse(size_t pos = 0u);  
    void calc_backend_dense(size_t pos = 0u);  

public:
    Array_Type                         EmptyArray;
    std::vector< Array_Type >          data;
    Rules<Array_Type,Data_Rule_Type> * rules;

    size_t N, M;
    bool rules_deleted   = false;

    // Tempvars
    std::vector< size_t >  coordinates_free;
    std::vector< size_t >  coordinates_locked;
    size_t n_free;
    size_t n_locked;
    
    /**
     * @name Construct and destroy a PowerSet object
     * 
     */
    ///@{
    PowerSet() : 
    EmptyArray(), data(0u), rules(new Rules<Array_Type,Data_Rule_Type>()), N(0u), M(0u) {};
    PowerSet(size_t N_, size_t M_) :
        EmptyArray(N_, M_), data(0u),
        rules(new Rules<Array_Type,Data_Rule_Type>()), N(N_), M(M_) {};
    PowerSet(const Array_Type & array);

    ~PowerSet();
    ///@}
    
    void init_support();
    void calc();
    void reset(size_t N_, size_t M_);
    
    /**
     * @name Wrappers for the `Rules` member. 
     * @details These will add rules to the model, which are shared by the
     * support and the actual counter function.
     */
    ///@{
    void add_rule(Rule<Array_Type, Data_Rule_Type> rule);
    void add_rule(
        Rule_fun_type<Array_Type,Data_Rule_Type> count_fun_,
        Data_Rule_Type data_
    );
    ///@}
    

    /** @name Getter functions */
    ///@{
    const std::vector< Array_Type > * get_data_ptr() const {return &data;};
    std::vector< Array_Type > get_data() const {return data;};
    typename std::vector< Array_Type >::iterator begin() {return data.begin();};
    typename std::vector< Array_Type >::iterator end() {return data.end();};
    std::size_t size() const noexcept {return data.size();};
    const Array_Type& operator[](const size_t & i) const {return data.at(i);};
    ///@}
    
};

#endif
