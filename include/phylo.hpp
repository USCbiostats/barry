#include "barray.hpp"

#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H

typedef unsigned int uint;

template <typename T>
using Vec = std::vector< T >;

// template <typename Cell_type>
class PhyloNode : public barray::BArray< bool > {
  
private:
  unsigned int nfun;
  unsigned int noff;
  
  /**@brief Branch length of the offsprings
   * 
   */
  Vec< double > blengths;
  
  /**@brief State of the parent node
   * 
   */
  Vec< bool > states;
  
public:
  // PhyloNode() : BArray<bool>(0u, 0u) {} ;
  PhyloNode() : BArray<bool>(0u,0u), nfun(0u), noff(0u), blengths(0u),
  states(0u) {};
  
  /**@brief Constructor assuming parent has all states in 1.
   * 
   */
  PhyloNode(
    const uint & nfun_,
    const uint & noff_
    ) : BArray<bool>(nfun_, noff_), nfun(nfun_), noff(noff_),
    blengths(noff_, 1.0), states(nfun_, true) {};
  
  /**@brief Constructor specifying the state of the parent node
   * 
   */
  PhyloNode(
    const Vec< bool > & states_,
    const uint & noff_
  ) : BArray<bool>(states_.size(), noff_), nfun(states_.size()), noff(noff_),
  blengths(noff_, 1.0), states(states_) {};
  
  ~PhyloNode() {};
  
};

/**@brief Extension of a simple counter.
 * 
 * It allows specifying extra arguments, in particular, the corresponding
 * sets of rows to which this statistic may be relevant. This could be important
 * in the case of, for example, counting correlation type statistics between
 * function 1 and 2, and between function 1 and 3.
 * 
 * 
 */
class PhyloCounter : public barray::Counter< bool > {
  
public:
  Vec< uint > fun_seq;
  
  /**@brief Sinplest case
   * 
   * Only the counter is provided, this is equivalent to having used a Counter
   * object instead
   * 
   */
  PhyloCounter(
    barray::Counter_type< bool > count_fun_
  ) : barray::Counter< bool >(count_fun_), fun_seq(0u) {};
  
  /**@brief Intermediate constructor
   * 
   * In this case, we are passing an init function, a counter and a function
   * sequence that can be used within the counter.
   */
  PhyloCounter(
    barray::Counter_type< bool > count_fun_,
    const Vec< uint > & fun_seq_
  ) : barray::Counter< bool >(count_fun_), fun_seq(fun_seq_) {};
  
  /**@brief More complex constructor
   * 
   * In this case, we are passing an init function, a counter and a function
   * sequence that can be used within the counter.
   */
  PhyloCounter(
    barray::Counter_type< bool > count_fun_,
    barray::Counter_type< bool > init_fun_,
    const Vec< uint > & fun_seq_
  ) : barray::Counter< bool >(count_fun_, init_fun_), fun_seq(fun_seq_) {};
  
  ~PhyloCounter() {};
  
  void set_seq(const Vec< uint > & fun_seq_) {this->fun_seq = fun_seq_;};
  
};

namespace phylo_counters {

  
  // Functional gains
  inline double count_gains(
      PhyloNode * Array,
      uint i,
      uint j,
      PhyloCounter * counter
  ) {
    
    if (i == counter->fun_seq[0u])
      return 1.0;
    else
      return 0.0;
  }
  
  inline double init_count_gains(
      PhyloNode * Array,
      uint i, uint j,
      PhyloCounter * counter
    ) {
    return (double) (Array->N);
  }
  
  Vec<uint> x(1, 0u);
  PhyloCounter gains(count_gains, init_count_gains);
}

#endif
