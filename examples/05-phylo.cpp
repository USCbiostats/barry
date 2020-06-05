#include <Rcpp.h>
#include "../include/barray.hpp"

using namespace Rcpp;

// #ifndef BARRAY_PHYLO_H
// #define BARRAY_PHYLO_H

class NodeData;
typedef barray::BArray<bool, NodeData> PhyloArray;

typedef unsigned int uint;

template <typename T>
using Vec = std::vector< T >;

class NodeData {
public:
  Vec< double > blengths;
  Vec< bool > states;
  NodeData() : blengths(0u), states(0u) {};
  NodeData(
    Vec< double > & blengths_,
    Vec< bool > & states_
  ) : blengths(blengths_), states(states_) {};
  NodeData(
    Vec< double > blengths_,
    Vec< bool > states_
  ) : blengths(blengths_), states(states_) {};
  ~NodeData() {};
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


#define PHYLO_COUNTER(a) inline double (a) (PhyloArray * Array, uint i, \
  uint j, Vec<uint> * data)                                         \
    
namespace phylo_counters {

   
  // Functional gains ----------------------------------------------------------
  PHYLO_COUNTER(count_gains) {
    return (!Array->data->states[j]) && (j == (*data)[0u]) ? 1.0 : 0.0;
  }

  PHYLO_COUNTER(init_count_gains) {
    
    if (data == nullptr)
      throw std::logic_error("count_gains needs to initialize the data.");
    
    return 0.0;
  }
  
  barray::Counter<PhyloArray, Vec<uint>> gains(count_gains, init_count_gains);

  // Functional loss ----------------------------------------------------------
  
  PHYLO_COUNTER(count_loss) {
    return (Array->data->states[j]) && (j == (*data)[0u]) ? -1.0 : 0.0;
  }
  
  PHYLO_COUNTER(init_count_loss) {
    
    if (data == nullptr)
      throw std::logic_error("count_loss needs to initialize the data.");
    
    return Array->data->states[(*data)[0u]]? Array->M : 0.0;
  }
  
  barray::Counter<PhyloArray, Vec<uint>> loss(count_loss, init_count_loss);
  
  // Sub-functionalization ----------------------------------------------------
  // It requires to specify data = {funA, funB}
  PHYLO_COUNTER(count_subfun) {
    
    // Are we looking at either of the relevant functions?
    if (((*data)[0u] != i) && ((*data)[1u] != i))
      return 0.0;
    
    // Are A and B existant? if not, no change
    if (!Array->data->states[(*data)[0u]] | !Array->data->states[(*data)[1u]])
      return 0.0;
    
    // Figuring out which is the first (reference) function
    uint * first;
    uint * second;
    uint * other;
    if (i == (*data)[0u]) {
      
      first  = &((*data)[0u]);
      second = &((*data)[1u]);
      other  = &(*second);
      
    } else {
      
      first  = &((*data)[1u]);
      second = &((*data)[0u]);
      other  = &(*first);
      
    }
    
    double res = 0.0;
    // There are 4 cases: (first x second) x (had the second function)
    if (!Array->is_empty(*other, j)) {
    
      for (uint off = 0u; off < Array->M; ++off) {
        
        // Not on self
        if (off == j)
          continue;
        
        if (!(Array->is_empty(*first, off)) && Array->is_empty(*second, off))
          res -= 1.0;
        
      }
      
    } else {
      
      for (uint off = 0u; off < Array->M; ++off) {
        
        // Not on self
        if (off == j)
          continue;
        
        if (Array->is_empty(*first, off) && !(Array->is_empty(*second, off)))
          res += 1.0;
        
      }
      
    }
    
    
    return res;
  }
  
  PHYLO_COUNTER(init_count_subfun) {
    
    if (data == nullptr)
      throw std::logic_error("subfun needs to initialize the data.");
    
    if (data->size() != 2u)
      throw std::logic_error("subfun should be initialized with a vec of size 2.");
    
    return 0.0;
  }
  
  barray::Counter<PhyloArray, Vec<uint>> subfun(count_subfun, init_count_subfun);

}

// [[Rcpp::export]]
List counter_phylo(const LogicalVector & x, int nfun, int noffspring) {
  
  // Initializing the node
  PhyloArray tree(nfun, noffspring);
  NodeData data({1.0, 1.0}, as<Vec<bool>>(x)); 
  tree.data = &data;
  
  // Setting counters, one per function
  barray::Counter<PhyloArray, Vec<uint>> counter0 = phylo_counters::gains;
  barray::Counter<PhyloArray, Vec<uint>> counter1 = phylo_counters::gains;
  counter0.data = new Vec<uint>({0u});
  counter1.data = new Vec<uint>({1u});
  
  barray::Counter<PhyloArray, Vec<uint>> counter2 = phylo_counters::loss;
  barray::Counter<PhyloArray, Vec<uint>> counter3 = phylo_counters::loss;
  counter2.data = new Vec<uint>({0u});
  counter3.data = new Vec<uint>({1u});
  
  barray::Counter<PhyloArray, Vec<uint>> counter4 = phylo_counters::subfun;
  counter4.data = new Vec<uint>({0u, 1u});
  
  
  barray::Support<PhyloArray, Vec<uint>> support(&tree);
  support.add_counter(counter0);
  support.add_counter(counter1);
  support.add_counter(counter2);
  support.add_counter(counter3);
  support.add_counter(counter4);
  
  // Computing and retrieving
  support.calc(0u, true);
  
  delete counter0.data;
  delete counter1.data;
  delete counter2.data;
  delete counter3.data;
  delete counter4.data;
  counter0.data = nullptr;
  counter1.data = nullptr;
  counter2.data = nullptr;
  counter3.data = nullptr;
  counter4.data = nullptr;
  
  // Generating the entries
  barray::Counts_type ans = support.support.get_entries();
  
  List res(ans.size());
  for (unsigned int i = 0u; i < res.size(); ++i) {
    res[i] = List::create(
      _["x"] = ans.at(i).first,
      _["count"] = ans.at(i).second
    );
  }
  
  return res;
  
}


/***R

# Case 1 (0,0)
ans <- counter_phylo(c(FALSE, FALSE), nfun = 2, noffspring = 2)
do.call(rbind, lapply(ans, unlist))

# Case 2 (1, 0)
ans <- counter_phylo(c(TRUE, FALSE), nfun = 2, noffspring = 2)
do.call(rbind, lapply(ans, unlist))

# Case 3 (0, 1)
ans <- counter_phylo(c(FALSE, TRUE), nfun = 2, noffspring = 2)
do.call(rbind, lapply(ans, unlist))

# Case 3 (1, 1)
ans <- counter_phylo(c(TRUE, TRUE), nfun = 2, noffspring = 2)
(ans <- do.call(rbind, lapply(ans, unlist)))
sum(ans[,"count"])



*/

// 
// #endif
