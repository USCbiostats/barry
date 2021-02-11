// #include <vector>
// #include <unordered_map>
#include "typedefs.hpp"
#include "barray-bones.hpp"
#include "statsdb.hpp"
#include "counters-bones.hpp"
#include "rules-bones.hpp"

#ifndef BARRY_SUPPORT_BONES_HPP 
#define BARRY_SUPPORT_BONES_HPP 1

/**@brief Compute the support of sufficient statistics
 * 
 * Given an array and a set of counters, this object iterates throughout the
 * support set of the Array while at the same time computing the support of
 * the sufficient statitics.
 */ 
template <
  typename Array_Type = BArray<>,
  typename Data_Counter_Type = bool,
  typename Data_Rule_Type = bool
  >
class Support {
  
private:
  void calc_backend(
      uint pos = 0u,
      std::vector< Array_Type > * array_bank = nullptr,
      std::vector< std::vector< double > > * stats_bank = nullptr
  );
  
public:
  
  /**@brief Reference array to generate the support.
   */
  Array_Type                                    EmptyArray;
  FreqTable<>                                   data;
  Counters<Array_Type,Data_Counter_Type> * counters;
  Rules<Array_Type,Data_Rule_Type>            * rules;

  uint N, M;
  bool counter_deleted = false;
  bool rules_deleted   = false;
  
  // Temp variables to reduce memory allocation
  std::vector< double >                current_stats;
  std::vector< std::pair<uint,uint> >  coordinates_free;
  std::vector< std::pair<uint,uint> >  coordinates_locked;
  std::vector< std::vector< double > > change_stats;
  // std::vector< double > change_stats;
  
  /**@brief Constructor passing a reference Array.
   */
  Support(const Array_Type * Array_) :
    EmptyArray(*Array_),
    counters(new Counters<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    N(Array_->nrow()), M(Array_->ncol()) {
    // init_support();
    return;
    
  };
  
  /**@brief Constructor specifying the dimensions of the array (empty).
   */
  Support(uint N_, uint M_) :
    EmptyArray(N_, M_),
    counters(new Counters<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    N(N_), M(M_) {
    // init_support();
    return;
  };
  
  Support() :
    EmptyArray(0u, 0u),
    counters(new Counters<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    N(0u), M(0u) {
    // init_support();
    return;
  };
  
  ~Support() {
    if (!counter_deleted)
      delete counters;
    if (!rules_deleted)
      delete rules;
  };
  
  void init_support(
    std::vector< Array_Type > * array_bank = nullptr,
    std::vector< std::vector< double > > * stats_bank = nullptr
  );
  
  
  /**@brief Resets the support calculator
   * 
   * If needed, the counters of a support object can be reused.
   * 
   * @param Array_ New array over which the support will be computed.
   */
  void reset_array();
  void reset_array(const Array_Type * Array_);
  
  void add_counter(Counter<Array_Type, Data_Counter_Type> * f_);
  void add_counter(Counter<Array_Type,Data_Counter_Type> f_);
  void set_counters(Counters<Array_Type,Data_Counter_Type> * counters_);
  
  void add_rule(Rule<Array_Type, Data_Rule_Type> * f_);
  void add_rule(Rule<Array_Type,Data_Rule_Type> f_);
  void set_rules(Rules<Array_Type,Data_Rule_Type> * rules_);
  
  /**@brief Computes the entire support
   * 
   * Not to be used by the user. Sets the starting point in the array
   * (column-major).
   * 
   * @param diag When `true`, includes the diagonal (i=j) in the counts.
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
      std::vector< std::vector< double > > * stats_bank = nullptr
    );
  
  
  Counts_type           get_counts() const;
  const MapVec_type<> * get_counts_ptr() const;
  void print() const;
  
};


#endif
