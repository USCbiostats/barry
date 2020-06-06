#include "../counters-bones.hpp"

#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H 1

class NodeData {
public:
  std::vector< double > blengths;
  std::vector< bool > states;
  NodeData() : blengths(0u), states(0u) {};
  NodeData(
    std::vector< double > & blengths_,
    std::vector< bool > & states_
  ) : blengths(blengths_), states(states_) {};
  NodeData(
    std::vector< double > blengths_,
    std::vector< bool > states_
  ) : blengths(blengths_), states(states_) {};
  ~NodeData() {};
};

typedef BArray<bool, NodeData> PhyloArray;


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
  uint j, std::vector<uint> * data)                                         \
    

// Functional gains ----------------------------------------------------------
PHYLO_COUNTER(count_gains) {
  return (!Array->data->states[i]) && (i == (*data)[0u]) ? 1.0 : 0.0;
}
 
PHYLO_COUNTER(init_count_gains) {
  
  if (data == nullptr)
    throw std::logic_error("count_gains needs to initialize the data.");
  
  return 0.0;
}

Counter<PhyloArray, std::vector<uint>> gains(count_gains, init_count_gains);

// Functional loss ----------------------------------------------------------

PHYLO_COUNTER(count_loss) {
  return (Array->data->states[i]) && (i == (*data)[0u]) ? -1.0 : 0.0;
}

PHYLO_COUNTER(init_count_loss) {
  
  if (data == nullptr)
    throw std::logic_error("count_loss needs to initialize the data.");
  
  return Array->data->states[(*data)[0u]]? Array->M : 0.0;
}

Counter<PhyloArray, std::vector<uint>> loss(count_loss, init_count_loss);

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
  uint other = (i == (*data)[0u])? (*data)[1u] : (*data)[0u];
  double res = 0.0;
  // There are 4 cases: (first x second) x (had the second function)
  if (!Array->is_empty(other, j)) {
  
    for (uint off = 0u; off < Array->M; ++off) {
      
      // Not on self
      if (off == j)
        continue;
      
      if (!(Array->is_empty(i, off)) && Array->is_empty(other, off))
        res -= 1.0;
       
    }
    
  } else {
    
    for (uint off = 0u; off < Array->M; ++off) {
      
      // Not on self
      if (off == j)
        continue;
      
      if (Array->is_empty(i, off) && !(Array->is_empty(other, off)))
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

Counter<PhyloArray, std::vector<uint>> subfun(count_subfun, init_count_subfun);

// Co-evolution (joint gain or loss) -----------------------------------------
PHYLO_COUNTER(count_cogain) {
  
  // Is the function in scope relevant?
  if ((i != (*data)[0u]) && (i != (*data)[1u]))
    return 0.0;
  
  // Both parents should either be zero or one
  if (!Array->data->states[0u] && Array->data->states[1u])
    return 0.0;
  
  if (Array->data->states[1u] && !Array->data->states[0u])
    return 0.0;
  
  // If both have it, then nothing is gained at this point
  if (Array->data->states[0u]) 
    return 0.0;
  
  uint other = (i == (*data)[0u])? (*data)[1u] : (*data)[0u];
  if (!Array->is_empty(other, j)) {
    return 1.0;
  }
  
  return 0.0;
  
}

PHYLO_COUNTER(init_count_cogain) {
  
  if (data == nullptr)
    throw std::logic_error("cogain needs to initialize the data.");
  
  if (data->size() != 2u)
    throw std::logic_error("cogain should be initialized with a vec of size 2.");
  
  return 0.0;
}

Counter<PhyloArray, std::vector<uint>> cogain(count_cogain, init_count_cogain);

// Longest branch mutates (either by gain or by loss) ------------------------
PHYLO_COUNTER(count_longest) {
  
  // Only relevant if the 
  double res = 0.0;
  if (Array->data->states[i]) {
    
    for (auto off = data->begin(); off != data->end(); ++off)
      if (*off == j) {
        res -= 1.0;
        break;
      }
    
  } else {
    
    for (auto off = data->begin(); off != data->end(); ++off)
      if (*off == j) {
        res += 1.0;
        break;
      }
    
  }
    
  return res;
  
}

PHYLO_COUNTER(init_count_longest) {
  
  if (Array->data == nullptr)
    throw std::logic_error("longest needs to initialize the data.");
  
  if (data == nullptr)
    throw std::logic_error(
    "data should be initialized in the counter (nullptr right now)."
    );
  
  if (data->size() != 0u)
    throw std::logic_error(
        "data should be initialized as a vector of size 0."
        );

  if (Array->data->blengths.size() != Array->M)
    throw std::logic_error(
    "longest should be initialized with a vec of size Array->M."
    );
  
  // Finding the longest branch (or branches) --
  uint longest_idx = 0u;
  double diff = 0.0;
  data->reserve(Array->M);
  data->push_back(0u);
  for (uint ii = 1u; ii < Array->M; ++ii) {
    
    diff = Array->data->blengths[longest_idx] - Array->data->blengths[ii];
    if (diff > 0.0) {
      continue;
    } else if (diff < 0.0) {
      data->empty();
      data->push_back(ii);
      longest_idx = ii;
    } else if (diff == 0.0) {
      data->push_back(ii);
    }
    
  }
  data->shrink_to_fit();
  
  if (data->size() == 0u)
    throw std::logic_error("The data on the longest branch has size 0.");
  
  // Starting the counter, since all in zero, then this will be equal to
  // the number of functions in 1 x number of longest branches
  double res = 0.0;
  for (uint ii = 0u; ii < Array->N; ++ii)
    if (Array->data->states[ii])
      res += (1.0 * data->size());
  
  
  return res;
}

Counter<PhyloArray, std::vector<uint>> longest(count_longest, init_count_longest);

#endif
