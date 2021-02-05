#include "../counters-bones.hpp"
#include "../support.hpp"
#include "../statscounter.hpp"
#include "../model-bones.hpp"

#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H 1

/**@brief Data definition for the `PhyloArray` class.
 * 
 * This holds basic information about a given node.
 * 
 */
class NodeData {
public:
  
  /**
   * Branch length.
   */
  std::vector< double > blengths;
  
  /**
   * State of the parent node.
   */
  std::vector< bool > states;
  
  /**
   * 
   */
  bool duplication = true;
  
  NodeData() : blengths(0u), states(0u) {};
  
  NodeData(
    const std::vector< double > & blengths_,
    const std::vector< bool > & states_,
    bool duplication_ = true
  ) : blengths(blengths_), states(states_), duplication(duplication_) {};
  
  ~NodeData() {};
  
};

typedef std::vector< uint > PhyloCounterData;
typedef std::vector< std::pair< uint, uint > > PhyloRuleData;

#define PHYLO_C_DATA_IDX(i) (data->operator[](i))
// #define PHYLO_C_DATA_NUM(i) (data->numbers[i])

/**@name Convenient typedefs for Node objects. */
///@{
typedef BArray<bool, NodeData> PhyloArray;
typedef Counter<PhyloArray, PhyloCounterData > PhyloCounter;
typedef CounterVector< PhyloArray, PhyloCounterData> PhyloCounterVector;
typedef Rule<PhyloArray,PhyloRuleData> PhyloRule;
typedef Rules<PhyloArray,PhyloRuleData> PhyloRules;
typedef Support<PhyloArray, PhyloCounterData, PhyloRuleData> PhyloSupport;
typedef StatsCounter<PhyloArray, PhyloCounterData> PhyloStatsCounter;
typedef Model<PhyloArray, PhyloCounterData, PhyloRuleData> PhyloModel;
///@}


/**@brief Extension of a simple counter.
 * 
 * It allows specifying extra arguments, in particular, the corresponding
 * sets of rows to which this statistic may be relevant. This could be important
 * in the case of, for example, counting correlation type statistics between
 * function 1 and 2, and between function 1 and 3.
 * 
 * 
 */
#define PHYLO_COUNTER(a) inline double (a) (const PhyloArray * Array, uint i, \
  uint j, std::vector<uint> * data)

#define PHYLO_COUNTER_LAMBDA(a) Counter_fun_type<PhyloArray, PhyloCounterData> a = \
  [](const PhyloArray * Array, uint i, uint j, PhyloCounterData * data)

/**@name Counters for phylogenetic modeling.
 * @param counters A pointer to a `PhyloCounterVector` object (`CounterVector`<`PhyloArray`, `PhyloCounterData`>).
 */
//@{
// -----------------------------------------------------------------------------
/**@brief Overall functional gains
 * @details Total number of gains (irrespective of the function).
 */
inline void counter_overall_gains(PhyloCounterVector * counters) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    return 1.0;
  };
  
  counters->add_counter(tmp_count);
  return;
  
}

// -----------------------------------------------------------------------------
/**@brief Functional gains for a specific function (`nfun`). */
inline void counter_gains(PhyloCounterVector * counters, std::vector<uint> nfun) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    return (!Array->data->states[i]) && (i == PHYLO_C_DATA_IDX(0u)) ? 1.0 : 0.0;
  };
  
  for (auto i = nfun.begin(); i != nfun.end(); ++i) {
    counters->add_counter(
        tmp_count, nullptr,
        new PhyloCounterData({*i}),
        true
    );
  }
  
  return;
  
}

// -----------------------------------------------------------------------------
/**@brief Overall functional loss */
inline void counter_overall_loss(PhyloCounterVector * counters, uint nfun) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    return -1.0;
  };
  
  PHYLO_COUNTER_LAMBDA(tmp_init) {
    return Array->N * Array->M;
  };
  
  counters->add_counter(
      tmp_count, tmp_init,
      new PhyloCounterData({nfun}),
      true
  );
  
  return;

}

  
// -----------------------------------------------------------------------------
/**@brief Total count of losses for an specific function. */
inline void counter_loss(PhyloCounterVector * counters, std::vector<uint> nfun) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    return (Array->data->states[i]) && (i == PHYLO_C_DATA_IDX(0u)) ? -1.0 : 0.0;
  };
  
  PHYLO_COUNTER_LAMBDA(tmp_init) {
    return Array->data->states[PHYLO_C_DATA_IDX(0u)]? Array->M : 0.0;
  };
  
  for (auto i = nfun.begin(); i != nfun.end(); ++i) {
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({*i}),
        true
    );
  }
  
  return;
  
}


// -----------------------------------------------------------------------------
/**@brief Total count of Sub-functionalization events.
 * @details It requires to specify data = {funA, funB}
 */
inline void counter_subfun(PhyloCounterVector * counters, uint nfunA, uint nfunB) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    // Are we looking at either of the relevant functions?
    if ((PHYLO_C_DATA_IDX(0u) != i) && (PHYLO_C_DATA_IDX(1u) != i))
      return 0.0;
    
    // Are A and B existant? if not, no change
    if (!Array->data->states[PHYLO_C_DATA_IDX(0u)] | !Array->data->states[PHYLO_C_DATA_IDX(1u)])
      return 0.0;
    
    // Figuring out which is the first (reference) function
    uint other = (i == PHYLO_C_DATA_IDX(0u))? PHYLO_C_DATA_IDX(1u) :PHYLO_C_DATA_IDX(0u);
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
  };
  
  counters->add_counter(
      tmp_count, nullptr,
      new PhyloCounterData({nfunA, nfunB}),
      true
  );
  
  return;
  
}

// -----------------------------------------------------------------------------
/**@brief Co-evolution (joint gain or loss)
 * @details Needs to specify pairs of functions (`nfunA`, `nfunB`).
 */
inline void counter_cogain(PhyloCounterVector * counters, uint nfunA, uint nfunB) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    
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

  };
  
  counters->add_counter(
      tmp_count, nullptr,
      new PhyloCounterData({nfunA, nfunB}),
      true
  );
  
  return;
  
}

// -----------------------------------------------------------------------------
/**@brief Longest branch mutates (either by gain or by loss) */
inline void counter_longest(PhyloCounterVector * counters) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
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
  };
  
  PHYLO_COUNTER_LAMBDA(tmp_init) {
    if (Array->data == nullptr)
      throw std::logic_error("longest needs to initialize the data.");
    
    if (data == nullptr)
      throw std::logic_error(
          "data should be initialized in the counter (nullptr right now)."
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
    for (uint ii = 0u; ii < Array->N; ++ii) {
      if (Array->data->states[ii])
        res += (1.0 * data->size());
    }
    
    return res;
  };
  
  counters->add_counter(
      tmp_count, tmp_init,
      new PhyloCounterData({}),
      true
  );
  
  return;
  
}

//------------------------------------------------------------------------------
/**@brief Total number of neofunctionalization events 
 * @details Needs to specify pairs of function.
 */
inline void counter_neofun(PhyloCounterVector * counters, uint nfunA, uint nfunB) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    // Is the function in scope relevant?
    if ((i != (*data)[0u]) && (i != (*data)[1u]))
      return 0.0;
    
    // Checking if the parent has both functions
    if (!Array->data->states[(*data)[0u]] && !Array->data->states[(*data)[1u]]) {
      return 0.0;
    } else if (Array->data->states[(*data)[0u]] && Array->data->states[(*data)[1u]]) {
      return 0.0;
    }
    
    // Figuring out which is the first (reference) function
    uint other = (i == (*data)[0u])? (*data)[1u] : (*data)[0u];
    double res = 0.0;
    
    if (Array->is_empty(other, j, false)) {
      
      for (auto off = 0u; off < Array->M; ++off) {
        if (off == j)
          continue;
        
        if (Array->is_empty(i, off, false) && !Array->is_empty(other, off, false))
          res += 1.0;
        
      }
      
    } else {
      
      for (auto off = 0u; off < Array->M; ++off) {
        if (off == j)
          continue;
        
        if (!Array->is_empty(i, off, false) && Array->is_empty(other, off, false))
          res -= 1.0;
        
      }
      
    }
    
    return res;
  };
  
  counters->add_counter(
      tmp_count, nullptr,
      new PhyloCounterData({nfunA, nfunB}),
      true
  );
  
  return;
  
}
//@}
  
#endif
