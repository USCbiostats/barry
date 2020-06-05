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
  uint j, uint * data)                                         \
    
namespace phylo_counters {

   
  // Functional gains ----------------------------------------------------------
  PHYLO_COUNTER(count_gains) {
    return (!Array->data->states[i]) && (j == *data) ? 1.0 : 0.0;
  }

  PHYLO_COUNTER(init_count_gains) {
    
    if (data == nullptr)
      throw std::logic_error("count_gains needs to initialize the data.");
    
    return 0.0;
  }
  
  barray::Counter<PhyloArray, uint> gains(count_gains, init_count_gains);

  // Functional loss ----------------------------------------------------------
  
  PHYLO_COUNTER(count_loss) {
    return (Array->data->states[i]) && (j == *data) ? -1.0 : 0.0;
  }
  
  PHYLO_COUNTER(init_count_loss) {
    
    if (data == nullptr)
      throw std::logic_error("count_loss needs to initialize the data.");
    
    return Array->data->states[*data]? Array->M : 0.0;
  }
  
  barray::Counter<PhyloArray, uint> loss(count_loss, init_count_loss);

}

// [[Rcpp::export]]
List counter_phylo(const LogicalVector & x) {
  
  // Initializing the node
  PhyloArray tree(2, 2);
  NodeData data({1.0, 1.0}, as<Vec<bool>>(x)); 
  tree.data = &data;
  
  // Setting counter 1
  barray::Counter<PhyloArray, uint> counter0 = phylo_counters::gains;
  counter0.data = new uint(0u);
  barray::Counter<PhyloArray, uint> counter1 = phylo_counters::gains;
  counter1.data = new uint(1u);
  
  barray::Support<PhyloArray, uint> support(&tree);
  support.add_counter(counter0);
  support.add_counter(counter1);
  
  // Computing and retrieving
  support.calc(0u, true);
  
  delete counter0.data;
  delete counter1.data;
  counter0.data = nullptr;
  counter1.data = nullptr;
  
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
ans <- counter_phylo(c(FALSE, FALSE))
do.call(rbind, lapply(ans, unlist))

# Case 2
ans <- counter_phylo(c(TRUE, FALSE))
do.call(rbind, lapply(ans, unlist))


*/

// 
// #endif
