#ifndef BARRY_RULES_BONES_HPP
#define BARRY_RULES_BONES_HPP 1

template <typename Array_Type, typename Data_Type>
bool rule_fun_default(const Array_Type * array, size_t i, size_t j, Data_Type * dat) {
    return false;
}

/**
  * @brief
  * Rule for determining if a cell should be included in a sequence
  * @details
  * Rules can be used together with `Support` and `PowerSet` to determine
  * which cells should be included when enumerating all possible realizations of
  * a binary array.
  * @tparam Array_Type An object of class `BArray`.
  * @tparam Data_Type Any type.
  */
template<typename Array_Type = BArray<>, typename Data_Type = bool>
class Rule {
    
private:
    Rule_fun_type<Array_Type,Data_Type> fun;
    Data_Type dat;
    
    std::string  name = "";
    std::string  desc = "";
    
public:

    /**
     * @name Construct a new Rule object
     * @brief Construct a new Rule object
     * 
     * @param fun_ A function of type `Rule_fun_type`.
     * @param dat_ Data pointer to be passed to `fun_`
     * @param delete_dat_ When `true`, the `Rule` destructor will delete the
     * pointer, if defined.
     */
    ///@{
    Rule() : fun(rule_fun_default<Array_Type,Data_Type>) {};
    Rule(
        Rule_fun_type<Array_Type,Data_Type> fun_,
        Data_Type dat_,
        std::string name_        = "",   
        std::string desc_        = ""
        ) : fun(fun_), dat(dat_), name(name_), desc(desc_) {};
    ///@}

    ~Rule() {};

    Data_Type & D(); ///< Read/Write access to the data.
    
    bool operator()(const Array_Type & a, size_t i, size_t j);

    std::string & get_name();
    std::string & get_description();

    std::string get_name() const;
    std::string get_description() const;
    
};

/**
  * @brief Vector of objects of class Rule
  * 
  * @tparam Array_Type An object of class `BArray`
  * @tparam Data_Type Any type.
  */
template<typename Array_Type, typename Data_Type>
class Rules {

private:
    std::vector< Rule<Array_Type,Data_Type> > data;
    
public:
    Rules() {};

    Rules(const Rules<Array_Type,Data_Type> & rules_);
    Rules<Array_Type,Data_Type> operator=(const Rules<Array_Type,Data_Type> & rules_);

    ~Rules() {};

    size_t size() const noexcept {
        return data.size();
    };
    
    /**
      * @name Rule adding
      * 
      * @param rule 
      */
    ///@{
    void add_rule(Rule<Array_Type, Data_Type> rule);
    void add_rule(
        Rule_fun_type<Array_Type,Data_Type> rule_,
        Data_Type data_,
        std::string name_ = "",
        std::string description_ = ""
    );
    ///@}

    /**
     * @brief Check whether a given cell is free or locked
     * 
     * @param a A `BArray` object
     * @param i row position
     * @param j col position
     * @return true If the cell is locked
     * @return false If the cell is free
     */

    bool operator()(const Array_Type & a, size_t i, size_t j);
    
    /**
     * @brief Computes the sequence of free and locked cells in an BArray
     * 
     * @param a An object of class `BArray`.
     * @param free Pointer to a vector of pairs (i, j) listing the free cells.
     * @param locked (optional) Pointer to a vector of pairs (i, j) listing the
     * locked cells.
     * @return Nothing.
     */
    void get_seq(
        const Array_Type & a,
        std::vector< size_t > * free,
        std::vector< size_t > * locked = nullptr
    );

    std::vector< std::string > get_names() const;
    std::vector< std::string > get_descriptions() const;

    // Iterator
    typename std::vector< Rule<Array_Type,Data_Type> >::iterator begin() {
        return data.begin();
    };
    typename std::vector< Rule<Array_Type,Data_Type> >::iterator end() {
        return data.end();
    };
    
};


#endif
