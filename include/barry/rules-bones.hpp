#include "typedefs.hpp"

#ifndef BARRY_RULES_BONES_HPP
#define BARRY_RULES_BONES_HPP 1

template <typename Array_Type, typename Data_Type>
bool rule_fun_default(const Array_Type * array, uint i, uint j, Data_Type * dat) {
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
    Data_Type * dat = nullptr;
    bool delete_dat = false;
    
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
        Data_Type * dat_ = nullptr,
        bool delete_dat_ = false
        ) : fun(fun_), dat(dat_), delete_dat(delete_dat_) {};
    ///@}

    ~Rule() {
        if (delete_dat)
            delete dat;
        return;
    }

    Data_Type * D(); ///< Read/Write access to the data.
    
    bool operator()(const Array_Type & a, uint i, uint j);
    
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
    std::vector< Rule<Array_Type,Data_Type> * > data = {};
    std::vector< uint > to_be_deleted                = {};
    
public:
    Rules() {};

    Rules(const Rules<Array_Type,Data_Type> & rules_);
    Rules<Array_Type,Data_Type> operator=(const Rules<Array_Type,Data_Type> & rules_);

    ~Rules() {
        this->clear();
        return;
    }

    uint size() const noexcept {
        return data.size();
    };
    
    /**
      * @name Rule adding
      * 
      * @param rule 
      */
    ///@{
    void add_rule(Rule<Array_Type, Data_Type> & rule);
    void add_rule(Rule<Array_Type, Data_Type> * rule);
    void add_rule(
        Rule_fun_type<Array_Type,Data_Type> rule_,
        Data_Type *                         data_        = nullptr,
        bool                                delete_data_ = false
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

    bool operator()(const Array_Type & a, uint i, uint j);
    
    void clear();
    
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
        std::vector< std::pair<uint,uint> > * free,
        std::vector< std::pair<uint,uint> > * locked = nullptr
    );
    
};


#endif
