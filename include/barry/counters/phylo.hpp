#include "../counters-bones.hpp"
#include "../support-bones.hpp"
#include "../statscounter-bones.hpp"
#include "../model-bones.hpp"
#include "../powerset-bones.hpp"

#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H 1

/**
 * @brief Data definition for the `PhyloArray` class.
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
  
  NodeData() {};
  
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

/**
 * @name Convenient typedefs for Node objects.
 * */
///@{
typedef BArray<uint, NodeData> PhyloArray;
typedef Counter<PhyloArray, PhyloCounterData > PhyloCounter;
typedef Counters< PhyloArray, PhyloCounterData> PhyloCounters;
typedef Rule<PhyloArray,PhyloRuleData> PhyloRule;
typedef Rules<PhyloArray,PhyloRuleData> PhyloRules;
typedef Support<PhyloArray, PhyloCounterData, PhyloRuleData> PhyloSupport;
typedef StatsCounter<PhyloArray, PhyloCounterData> PhyloStatsCounter;
typedef Model<PhyloArray, PhyloCounterData, PhyloRuleData> PhyloModel;
typedef PowerSet<PhyloArray, PhyloRuleData> PhyloPowerSet;
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

/**
 * @name Counters for phylogenetic modeling.
 * @param counters A pointer to a `PhyloCounters` object (`Counters`<`PhyloArray`, `PhyloCounterData`>).
 */
//@{
// -----------------------------------------------------------------------------
/**
 * @brief Overall functional gains
 * @details Total number of gains (irrespective of the function).
 */
inline void counter_overall_gains(PhyloCounters * counters, bool duplication = true) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    if ((data->at(0u) == 1u) & Array->data->duplication)
      return 1.0;
    else if ((data->at(0u) == 0u) & !Array->data->duplication) {
      return 1.0;
    }

    return 0.0;
    
  };
  
  counters->add_counter(
    tmp_count, nullptr,
    new PhyloCounterData({duplication ? 1u : 0u}),
    true
  );

  return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Functional gains for a specific function (`nfun`).
 */
inline void counter_gains(PhyloCounters * counters, std::vector<uint> nfun, bool duplication = true) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {

    if (Array->data->duplication & (data->at(1u) == 0u))
      return 0.0;
    else if (!Array->data->duplication & (data->at(1u) == 1u))
      return 0.0;
    
    return (!Array->data->states[i]) && (i == data->at(0u)) ? 1.0 : 0.0;

  };
  
  for (auto i = nfun.begin(); i != nfun.end(); ++i) {
    counters->add_counter(
        tmp_count, nullptr,
        new PhyloCounterData({*i, duplication ? 1u : 0u}),
        true
    );
  }
  
  return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Overall functional loss
 */
inline void counter_overall_loss(PhyloCounters * counters, bool duplication = true) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    
    if ((data->at(0u) == 1u) & Array->data->duplication)
      return -1.0;
    else if ((data->at(0u) == 0u) & !Array->data->duplication) {
      return -1.0;
    } else {
      return 0.0;
    }
  };
  
  PHYLO_COUNTER_LAMBDA(tmp_init) {

    if ((data->at(0u) == 1u) & Array->data->duplication)
      return static_cast<double>((Array->N * Array->M));
    else if ((data->at(0u) == 0u) & !Array->data->duplication)
      return static_cast<double>((Array->N * Array->M));
    else 
      return 0.0;

  };
  
  counters->add_counter(
      tmp_count, tmp_init,
      new PhyloCounterData({duplication ? 1u : 0u}),
      true
  );
  
  return;

}

// -----------------------------------------------------------------------------
/**
 * @brief Cap the number of functions per gene
 */
inline void counter_maxfuns(
  PhyloCounters * counters,
  uint            lb,
  uint            ub,
  bool duplication = true
  ) {

  PHYLO_COUNTER_LAMBDA(tmp_init) {
    if (data->at(0u) == 0u)
      return static_cast<double>(Array->ncol());
    else {

      double ans = 0.0;
      for (uint j = 0u; j < Array->ncol(); ++j) {

        uint count = 0u;
        for (uint i = 0u; i < Array->nrow(); ++i) {
          if (!Array->is_empty(i, j))
            ++count;
        }

        if (count >= data->at(0u) && count <= data->at(1u))
          ans += 1.0;

      }

      return ans;

    }
  };
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {

    if (Array->data->duplication & (data->at(2u) == 0u))
      return 0.0;
    else if (!Array->data->duplication & (data->at(2u) == 1u))
      return 0.0;
    
    // Does the focal gene has nfun in [lb,ub]?
    if (data->at(1u) == 0u) {

      return
        Array->el_ji.at(j).size() <= data->at(0u) ?
        0.0 : -1.0;

    } else {

      // Right above the lb?
      if (Array->el_ji.at(j).size() == data->at(0u))
        return 1.0;
      else if (Array->el_ji.at(j).size() == (data->at(1u) + 1u))
        return -1.0;
      else
        return 0.0;

    }

  };

  counters->add_counter(
      tmp_count, tmp_init,
      new PhyloCounterData({lb, ub, duplication ? 1u : 0u}),
      true
  );
  
  return;
  
}
  
// -----------------------------------------------------------------------------
/**
 * @brief Total count of losses for an specific function.
 */
inline void counter_loss(PhyloCounters * counters, std::vector<uint> nfun, bool duplication = true) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {

    if ((data->at(1u) == 1u) & !Array->data->duplication)
      return 0.0;
    else if ((data->at(1u) == 0u) & Array->data->duplication)
      return 0.0;
    else
      return (Array->data->states[i]) && (i == data->at(0u)) ? -1.0 : 0.0;

  };
  
  PHYLO_COUNTER_LAMBDA(tmp_init) {

    if ((data->at(1u) == 1u) & !Array->data->duplication)
      return 0.0;
    else if ((data->at(1u) == 0u) & Array->data->duplication)
      return 0.0;
    else
      return Array->data->states[data->at(0u)]? Array->M : 0.0;

  };
  
  for (auto& i : nfun) {
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({i, duplication ? 1u : 0u}),
        true
    );
  }
  
  return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief Total count of Sub-functionalization events.
 * @details It requires to specify data = {funA, funB}
 */
inline void counter_subfun(PhyloCounters * counters, uint nfunA, uint nfunB, bool duplication = true) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {

    // Is this node duplication?
    if ((data->at(2u) == 1u) & !Array->data->duplication)
      return 0.0;
    else if ((data->at(2u) == 0u) & Array->data->duplication)  
      return 0.0;

    // Are we looking at either of the relevant functions?
    if ((data->at(0u) != i) && (data->at(1u) != i))
      return 0.0;
    
    // Are A and B existant? if not, no change
    if (!Array->data->states[data->at(0u)] | !Array->data->states[data->at(1u)])
      return 0.0;
    
    // Figuring out which is the first (reference) function
    uint other = (i == data->at(0u))? data->at(1u) :data->at(0u);
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
      new PhyloCounterData({nfunA, nfunB, duplication ? 1u : 0u}),
      true
  );
  
  return;
  
}

// -----------------------------------------------------------------------------
/**@brief Co-evolution (joint gain or loss)
 * @details Needs to specify pairs of functions (`nfunA`, `nfunB`).
 */
inline void counter_cogain(PhyloCounters * counters, uint nfunA, uint nfunB, bool duplication = true) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {
    
        // Is this node duplication?
    if ((data->at(2u) == 1u) & !Array->data->duplication)
      return 0.0;
    else if ((data->at(2u) == 0u) & Array->data->duplication)  
      return 0.0;

    // Is the function in scope relevant?
    if ((i != data->at(0u)) && (i != data->at(1u)))
      return 0.0;
    
    // None should have it
    if (!Array->data->states[data->at(0u)] && !Array->data->states[data->at(1u)]) {

      uint other = (i == data->at(0u))? data->at(1u) : data->at(0u);

      if (Array->get_cell(other, j) == 1u) {

        return 1.0;

      } else {

        return 0.0;

      }

    } else {

      return 0.0;

    }

  };
  
  counters->add_counter(
      tmp_count, nullptr,
      new PhyloCounterData({nfunA, nfunB, duplication ? 1u : 0u}),
      true
  );
  
  return;
  
}

// -----------------------------------------------------------------------------
/**@brief Longest branch mutates (either by gain or by loss) */
inline void counter_longest(PhyloCounters * counters) {
  
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
inline void counter_neofun(PhyloCounters * counters, uint nfunA, uint nfunB, bool duplication = true) {
  
  PHYLO_COUNTER_LAMBDA(tmp_count) {

    // Is this node duplication?
    if ((data->at(2u) == 1u) & !Array->data->duplication)
      return 0.0;
    else if ((data->at(2u) == 0u) & Array->data->duplication)  
      return 0.0;

    // Is the function in scope relevant?
    if ((i != data->at(0u)) && (i != data->at(1u)))
      return 0.0;
    
    // Checking if the parent has both functions
    if (!Array->data->states[data->at(0u)] && !Array->data->states[data->at(1u)]) {
      return 0.0;
    } else if (Array->data->states[data->at(0u)] && Array->data->states[data->at(1u)]) {
      return 0.0;
    }
    
    // Figuring out which is the first (reference) function
    uint other = (i == data->at(0u))? data->at(1u) : data->at(0u);
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
      new PhyloCounterData({nfunA, nfunB, duplication ? 1u : 0u}),
      true
  );
  
  return;
  
}
//@}
  
#endif
