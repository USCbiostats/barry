#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H 1

// Default value that is used for the counters.
#define DEFAULT_DUPLICATION 2u
#define DUPL_SPEC 0u
#define DUPL_DUPL 1u
#define DUPL_EITH 2u


#define MAKE_DUPL_VARS() bool DPL = Array.D()->duplication; unsigned int DATA_AT = data->at(0u);
#define IS_EITHER() (DATA_AT == DUPL_EITH)
#define IS_DUPLICATION() ((DATA_AT == DUPL_DUPL) && (DPL))
#define IS_SPECIATION() ((DATA_AT == DUPL_SPEC) && (!DPL))

#define IF_MATCHES() MAKE_DUPL_VARS() \
    if (IS_EITHER() || IS_DUPLICATION() || IS_SPECIATION())
#define IF_NOTMATCHES() MAKE_DUPL_VARS() \
    if (!IS_EITHER() && !IS_DUPLICATION() && !IS_SPECIATION())


/**
 * @ingroup counting
 * @details Details about the available counters for `PhyloArray`
 * objects can be found in the \ref counters-phylo section.
 */
///@{

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
    std::vector< double > blengths = {};
    
    /**
     * State of the parent node.
     */
    std::vector< bool > states = {};
    
    /**
     * 
     */
    bool duplication = true;
    
    // NodeData() : blengths(0u), states(0u) {};
    
    NodeData(
        const std::vector< double > & blengths_,
        const std::vector< bool > & states_,
        bool duplication_ = true
    ) : blengths(blengths_), states(states_), duplication(duplication_) {};
    
    // ~NodeData() {};
  
};

// typedef std::vector< uint > PhyloCounterData;
class PhyloCounterData {
private:
    std::vector< uint > data;
    std::vector< double > * counters;

public:
    PhyloCounterData(
        std::vector< uint > data_,
        std::vector< double > * counters_ = nullptr
        ) : data(data_), counters(counters_) {};

    uint at(uint d) {return data.at(d);};
    uint operator()(uint d) {return data.at(d);};
    void reserve(uint x) {return data.reserve(x);};
    void push_back(uint x) {return data.push_back(x);};
    void shrink_to_fit()  {return data.shrink_to_fit();};
    uint size() {return data.size();};

    std::vector< uint >::iterator begin() {return data.begin();};
    std::vector< uint >::iterator end() {return data.end();};

    bool empty() {return data.empty();};
    std::vector< double > * get_counters() {return counters;};

};


typedef std::vector< std::pair< uint, uint > > PhyloRuleData;
class PhyloRuleDynData;

/**
 * @name Convenient typedefs for Node objects.
 * */
///@{
typedef BArray<uint, NodeData> PhyloArray;
typedef Counter<PhyloArray, PhyloCounterData > PhyloCounter;
typedef Counters< PhyloArray, PhyloCounterData> PhyloCounters;

typedef Rule<PhyloArray,PhyloRuleData> PhyloRule;
typedef Rules<PhyloArray,PhyloRuleData> PhyloRules;

typedef Rule<PhyloArray,PhyloRuleDynData> PhyloRuleDyn;
typedef Rules<PhyloArray,PhyloRuleDynData> PhyloRulesDyn;

typedef Support<PhyloArray, PhyloCounterData, PhyloRuleData, PhyloRuleDynData > PhyloSupport;
typedef StatsCounter<PhyloArray, PhyloCounterData> PhyloStatsCounter;
typedef Model<PhyloArray, PhyloCounterData, PhyloRuleData, PhyloRuleDynData > PhyloModel;
typedef PowerSet<PhyloArray, PhyloRuleData> PhyloPowerSet;
///@}


/**
 * @brief Extension of a simple counter.
 * 
 * It allows specifying extra arguments, in particular, the corresponding
 * sets of rows to which this statistic may be relevant. This could be important
 * in the case of, for example, counting correlation type statistics between
 * function 1 and 2, and between function 1 and 3.
 * 
 * 
 */
#define PHYLO_COUNTER_LAMBDA(a) Counter_fun_type<PhyloArray, PhyloCounterData> a = \
    [](const PhyloArray & Array, uint i, uint j, PhyloCounterData * data)

#define PHYLO_RULE_DYN_LAMBDA(a) Rule_fun_type<PhyloArray, PhyloRuleDynData> a = \
    [](const PhyloArray & Array, uint i, uint j, PhyloRuleDynData * data)

#define PHYLO_CHECK_MISSING() if (Array.D() == nullptr) \
    throw std::logic_error("The array data is nullptr."); \
    if (data == nullptr) \
    throw std::logic_error("The counter/rule data is nullptr.")

inline std::string get_last_name(unsigned int d) {return ((d == 1u)? " at duplication" : ((d == 0u)? " at speciation" : ""));}

/**
 * @weakgroup counters-phylo Phylo counters
 * @brief Counters for phylogenetic modeling
 * @param counters A pointer to a `PhyloCounters` object (`Counters`<`PhyloArray`, `PhyloCounterData`>).
 */
///@{
// -----------------------------------------------------------------------------
/**
 * @brief Overall functional gains
 * @details Total number of gains (irrespective of the function).
 */
inline void counter_overall_gains(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {
        IF_MATCHES()
            return 1.0;
        
        return 0.0;
      
    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication}),
        true,
        "Overall gains" + get_last_name(duplication)
    );

    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Functional gains for a specific function (`nfun`).
 */
inline void counter_gains(
    PhyloCounters * counters,
    std::vector<uint> nfun,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_MATCHES()
            return (i == data->at(1u)) ? 1.0 : 0.0;
        
        return 0.0;

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init,
            new PhyloCounterData({duplication, i}),
            true,
            "Gains " + std::to_string(i) + get_last_name(duplication)
        );
    
    return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief k genes gain function nfun
 */
inline void counter_gains_k_offspring(
    PhyloCounters * counters,
    std::vector<uint> nfun,
    uint k = 1u,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

      PHYLO_CHECK_MISSING();
      return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this relevant?
        if (i != data->at(1u))
            return 0.0;

        IF_NOTMATCHES()
            return 0.0;

        // Is there any gain?
        if (Array.D()->states[i])
            return 0.0;

        // Making the counts
        int counts = 0;
        for (uint k = 0u; k < Array.ncol(); ++k)
            if (k != j)
            {
                if (Array(i, k, false) == 1u)
                    ++counts;
            }

        // Three cases: base on the diff
        int diff = static_cast<int>(data->at(2u)) - counts + 1;
        // (a) counts were 1 below k, then +1
        if (diff == 1)
            return -1.0;
            // (b) counts were equal to k, then -1
        else if (diff == 0)
        {
            return 1.0;
        } else 
            // (c) Otherwise, nothing happens
            return 0.0;
      

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init,
            new PhyloCounterData({duplication, i, k}),
            true,
            std::to_string(k) + " genes gain " + std::to_string(i) +
                get_last_name(duplication)
        );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Keeps track of how many genes are changing (either 0, 1, or 2 if dealing
 * with regular trees.)
 */
inline void counter_genes_changing(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        for (uint j0 = 0u; j0 < Array.nrow(); ++j0)
        {

            if (Array.D()->states[j0]) 
                // Yup, we are loosing a function, so break
                return static_cast<double>(Array.ncol());
            
        }

        return 0.0;
      

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;

        // Need to check the other functions
        for (uint k = 0u; k < Array.nrow(); ++k)
        {

            // Nah, this gene was already different.
            if ((k != i) && (Array.D()->states[k] != (Array(k, j, false) == 1u)))
                return 0.0;
            

        }

        // Nope, this gene is now matching its parent, so we need to 
        // take it out from the count of genes that have changed.
        return Array.D()->states[i] ? -1.0 : 1.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication}),
        true,
        "Num. of genes changing" + get_last_name(duplication)
    );

    
    return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief Keeps track of how many genes are changing (either 0, 1, or 2 if dealing
 * with regular trees.)
 */
inline void counter_prop_genes_changing(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.

        for (uint j0 = 0u; j0 < Array.nrow(); ++j0)
        {

            if (Array.D()->states[j0]) 
                // Yup, we are loosing a function, so break
                return 1.0; // static_cast<double>(Array.ncol());
            
        }

        return 0.0;
      

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;

        // Case 1: The parent had the function (then probably need to substract one)
        if (Array.D()->states[i]) {

            // Need to check the other functions
            for (uint k = 0u; k < Array.nrow(); ++k)
            {

                if (k != i)
                {

                    // Nah, this gene was already different.
                    if (Array.D()->states[k] && (Array(k, j, false) == 0u))
                        return 0.0;
                    else if ((!Array.D()->states[k]) && (Array(k, j, false) == 1u))
                        return 0.0;

                }

            }

            // Nope, this gene is now matching its parent, so we need to 
            // take it out from the count of genes that have changed.
            return -1.0/static_cast<double>(Array.ncol());

        }
        else if (!Array.D()->states[i])
        {
            // Case 2: The parent didn't had the function. Probably need to increase
            // by one.


              // Need to check the other functions, where these the same?
              // if these were the same, then we are facing a gene who is changing.
              for (uint k = 0u; k < Array.nrow(); ++k)
              {

                  if (k != i)
                  {
                      // Nah, this gene was already different.
                      if (Array.D()->states[k] && (Array(k, j, false) == 0u))
                          return 0.0;
                      else if ((!Array.D()->states[k]) && (Array(k, j, false) == 1u))
                          return 0.0;
                  }
                  
              }

              // Nope, this gene is now matching its parent, so we need to 
              // take it out from the count of genes that have changed.
              return 1.0/static_cast<double>(Array.ncol());

        } else
            throw std::logic_error(
                "Reach the end of -counter_prop_genes_changing-. This shouldn't happen!"
                );

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication}),
        true,
        "Proportion of genes changing" + get_last_name(duplication)
    );

    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Overall functional loss
 */
inline void counter_overall_loss(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
    )
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_MATCHES()
            return -1.0;
        else 
            return 0.0;
        
    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();

        IF_MATCHES()
            return static_cast<double>((Array.nrow() * Array.ncol()));
        else 
            return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication}),
        true,
        "Overall loses" + get_last_name(duplication)
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
    unsigned int duplication = DEFAULT_DUPLICATION
 )
 {

    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

      PHYLO_CHECK_MISSING();    

      if (data->at(0u) == 0u)
        return static_cast<double>(Array.ncol());
      else
      {

        double ans = 0.0;
        for (uint k = 0u; k < Array.ncol(); ++k)
        {

          // How many functions the k-th offspring has
          uint count = 0u;
          for (uint l = 0u; l < Array.nrow(); ++l)
          {

            if (Array(l, k, false) == 1u)
              ++count;

          }

          if (count >= data->at(0u) && count <= data->at(1u))
            ans += 1.0;

        }

        return ans;

      }

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;
        
        uint counts = 1u;
        for (uint k = 0u; k < Array.nrow(); ++k)
            if (k != j)
                if (Array(k, j, false) == 1u)
                    ++counts;

        // Reached the lower bound
        if (counts == data->at(1u))
            return 1.0;
          // Went outside of the upper bound
        else if (counts == (data->at(2u) + 1u))
            return -1.0;
        else
            return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication, lb, ub}),
        true,
        "Genes with [" + std::to_string(lb) + ", " + std::to_string(ub) +
            "] funs" + get_last_name(duplication)
    );
    
    return;
  
}
  
// -----------------------------------------------------------------------------
/**
 * @brief Total count of losses for an specific function.
 */
inline void counter_loss(
    PhyloCounters * counters, std::vector<uint> nfun,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        if (!Array.D()->states[i])
            return 0.0;
        
        return (i == data->at(1u)) ? -1.0 : 0.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;
        
        if (!Array.D()->states[i])
            return 0.0;
        
        return Array.D()->states[data->at(1u)]? Array.ncol() : 0.0;

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init,
            new PhyloCounterData({duplication, i}),
            true,
            "Loss " + std::to_string(i) + get_last_name(duplication)
        );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Total number of changes. Use this statistic to account for "preservation"
 */
inline void counter_overall_changes(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        if (Array.D()->states[i])
            return -1.0;
        else 
            return 1.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        IF_NOTMATCHES()
            return 0.0;

        PHYLO_CHECK_MISSING();


        // Since we start with all the array at zero,
        // As many chances to change as offspring
        double noff   = static_cast<double> (Array.ncol());
        double counts = 0.0;
        for (uint k = 0u; k < Array.nrow(); ++k)
            if (Array.D()->states[k])
                counts += noff;

        return counts;

        

    };

    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication}),
        true,
        "Overall changes" + get_last_name(duplication)
    );
    
    
    return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief Total count of Sub-functionalization events.
 * @details It requires to specify data = {funA, funB}
 */
inline void counter_subfun(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;
        
        // Are we looking at either of the relevant functions?
        if ((data->at(1u) != i) && (data->at(2u) != i))
            return 0.0;
        
        // Are A and B existant? if not, no change
        if (!Array.D()->states[data->at(1u)] | !Array.D()->states[data->at(2u)])
            return 0.0;
        
        // Figuring out which is the first (reference) function
        uint other = (i == data->at(1u))? data->at(2u) : data->at(1u);
        double res = 0.0;
        // There are 4 cases: (first x second) x (had the second function)
        if (Array(other, j, false) == 1u)
        { 
          
            for (uint off = 0u; off < Array.ncol(); ++off)
            {
                
                // Not on self
                if (off == j)
                    continue;
                
                if ((Array(i, off, false) == 1u) && (Array(other, off, false) == 0u))
                    res -= 1.0;
                
            }
          
        } else {
          
            for (uint off = 0u; off < Array.ncol(); ++off)
            {
              
                // Not on self
                if (off == j)
                    continue;
                
                if ((Array(i, off, false) == 0u) && (Array(other, off, false) == 1u))
                    res += 1.0;
              
            }
          
        }
        
        return res;

    };

    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication, nfunA, nfunB}),
        true,
        "Subfun between " + std::to_string(nfunA) + " and " +
            std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Co-evolution (joint gain or loss)
 * @details Needs to specify pairs of functions (`nfunA`, `nfunB`).
 */
inline void counter_cogain(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        auto d1 = data->at(1u);
        auto d2 = data->at(2u);
      
        // Is the function in scope relevant?
        if ((i != d1) && (i != d2))
            return 0.0;
        
        // None should have it
        if (!Array.D()->states[d1] && !Array.D()->states[d2])
        {

            uint other = (i == d1)? d2 : d1;

            if (Array(other, j, false) == 1u)
                return 1.0;
            else
                return 0.0;

        } else 
            return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication, nfunA, nfunB}),
        true,
        "Co-gains " + std::to_string(nfunA) + " & " + std::to_string(nfunB) +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**@brief Longest branch mutates (either by gain or by loss) */
inline void counter_longest(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
    )
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // Only relevant if the 
        double res = 0.0;
        if (Array.D()->states[i])
        {
            
            for (auto& off : *data)
                if (off == j)
                {

                    res -= 1.0;
                    break;

                }
            
        } else {
            
            for (auto& off : *data)
                if (off == j)
                {

                    res += 1.0;
                    break;

                }
            
        }
        
        return res;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        
        if (Array.D()->blengths.size() != Array.ncol())
            throw std::logic_error(
                "longest should be initialized with a vec of size Array.ncol()."
            );
          
        // Finding the longest branch (or branches) --
        uint longest_idx = 0u;
        double diff      = 0.0;
        data->reserve(Array.ncol()); 
        data->push_back(0u);
        for (uint ii = 1u; ii < Array.ncol(); ++ii)
        {
            
            diff = Array.D()->blengths[longest_idx] - Array.D()->blengths[ii];
            if (diff > 0.0)
                continue;
            else if (diff < 0.0)
            {

                data->empty();
                data->push_back(ii);
                longest_idx = ii;

            }
            else if (diff == 0.0)
                data->push_back(ii);
            
        }

        data->shrink_to_fit();
        
        if (data->size() == 0u)
            throw std::logic_error("The data on the longest branch has size 0.");
        
        // Starting the counter, since all in zero, then this will be equal to
        // the number of functions in 1 x number of longest branches
        double res = 0.0;
        for (uint ii = 0u; ii < Array.nrow(); ++ii)
        {
            
            if (Array.D()->states[ii])
                res += (1.0 * data->size());

        }
        
        return res;
    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication}),
        true,
        "Longest branch mutates" + get_last_name(duplication)
    );
    
    return;
  
}

//------------------------------------------------------------------------------
/**
 * @brief Total number of neofunctionalization events 
 * @details Needs to specify pairs of function.
 */
inline void counter_neofun(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;
        
        // Is the function in scope relevant?
        if ((i != data->at(1u)) && (i != data->at(2u)))
            return 0.0;
        
        // Checking if the parent has both functions
        if (!Array.D()->states[data->at(1u)] && !Array.D()->states[data->at(2u)]) 
            return 0.0;
        else if (Array.D()->states[data->at(1u)] && Array.D()->states[data->at(2u)])
            return 0.0;
        
        // Figuring out which is the first (reference) function
        uint other = (i == data->at(1u))? data->at(2u) : data->at(1u);
        double res = 0.0;
        
        if (Array.is_empty(other, j, false))
        {
            
            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;
                
                if (Array.is_empty(i, off, false) && !Array.is_empty(other, off, false))
                    res += 1.0;
                
            }
          
        } else {
            
            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;
              
                if (!Array.is_empty(i, off, false) && Array.is_empty(other, off, false))
                    res -= 1.0;
              
            }
            
        }
        
        return res;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication, nfunA, nfunB}),
        true,
        "Neofun between " + std::to_string(nfunA) + " and " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
  
}

//------------------------------------------------------------------------------
/**
 * @brief Total number of neofunctionalization events 
 * @details Needs to specify pairs of function.
 */
inline void counter_neofun_a2b(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        if ((data->at(2u) == 1u) & !Array.D()->duplication)
            return 0.0;
        else if ((data->at(2u) == 0u) & Array.D()->duplication)  
            return 0.0;

        const uint & funA = data->at(0u);
        const uint & funB = data->at(1u);
        
        // Checking the parent has funA but not funb
        if ((!Array.D()->states[funA]) | Array.D()->states[funB]) 
            return 0.0;
      
        double res = 0.0;

        if (funA == i)
        {

            if (Array.get_cell(funB, j, true) == 1)
            {

                for (uint k = 0u; k < Array.ncol(); ++k)
                {

                    if (k == j)
                        continue;
                    if ((Array(funA, k, false) == 1u) && (Array(funB, k, false) == 0u))
                        res -= 1.0;

                }

            } else {

                for (uint k = 0u; k < Array.ncol(); ++k) {

                    if (k == j)
                        continue;
                    if ((Array(funA, k, false) == 0u) && (Array(funB, k, false) == 1u))
                        res += 1.0;

                }

            }

        } else {

            if (Array.get_cell(funB, j, true) == 1)
            {

                for (uint k = 0u; k < Array.ncol(); ++k)
                {

                    if (k == j)
                        continue;
                    if ((Array(funA, k, false) == 0u) && (Array(funB, k, false) == 1u))
                        res -= 1.0;

                }

            } else {

                for (uint k = 0u; k < Array.ncol(); ++k)
                {

                    if (k == j)
                        continue;
                    if ((Array(funA, k, false) == 1u) && (Array(funB, k, false) == 0u))
                        res += 1.0;

                }

            }

        }

        return res;
        
    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({nfunA, nfunB, duplication ? 1u : 0u}),
        true,
        "Neofun from " + std::to_string(nfunA) + " to " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
    
}

// -----------------------------------------------------------------------------
/**
 * @brief Function co-opting
 * @details Function co-opting of functions A and B happens when, for example,
 * function B is gained as a new featured leveraging what function A already does;
 * without losing function A. The sufficient statistic is defined as follows:
 * \f[
 * x_{pa}(1 - x_{pb})\sum_{i<j}\left[x_{ia}^p(1 - x_{ib}^p)x_{ja}^px_{jb}^p + x_{ja}^p(1 - x_{jb}^p)x_{ia}^px_{ib}^p\right]
 * \f]
 * This algorithm implements the change statistic.
 */
inline void counter_co_opt(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB, 
    unsigned int duplication = DEFAULT_DUPLICATION
) {
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    { 

        // Checking whether this is for duplication or not
        if ((data->at(2u) == 0u) & Array.D()->duplication)
            return 0.0;
        else if ((data->at(2u) == 1u) & !Array.D()->duplication)
            return 0.0;
        else {

            const unsigned int funA = data->at(0u);
            const unsigned int funB = data->at(1u);

            // If the change is out of scope, then nothing to do
            if ((i != funA) & (i != funB))
                return 0.0;

            // If the parent does not have the initial state, then it makes no sense
            if ((!Array.D()->states[funA]) | Array.D()->states[funB])
                return 0.0;

            // Checking whether function A or function B changed
            if (i == funA) {

                // What was the state of the other function? If B is present, then
                // nothing changes.
                if (Array(funB, j, false) == 1u) 
                    return 0.0;

                // Iterating through the sibs
                double res = 0.0;
                for (auto c = 0u; c < Array.ncol(); ++c)
                    if ((c != j) && (Array(funA, c, false) == 1u) && (Array(funB, c, false) == 1u))
                        res += 1.0;

                return res;

            } else {

                // What was the state of the other function? If A is not present, then
                // nothing changes.
                if (Array(funA, j, false) == 0u) 
                    return 0.0;

                // Iterating through the sibs
                double res = 0.0;
                for (auto c = 0u; c < Array.ncol(); ++c)
                    if ((c != j) && (Array(funA, c, false) == 1u))
                        res += (Array(funB, c, false) == 0u) ? 1.0 : -1.0;

                return res;

            }

        }

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        if (data->size() != 3u)
            throw std::length_error("The counter data should be of length 2.");

        if (data->at(0u) == data->at(1u))
            throw std::logic_error("Functions A and B should be different from each other.");

        if (data->at(0u) >= Array.nrow())
            throw std::length_error("Function A in counter out of range.");

        if (data->at(1u) >= Array.nrow())
            throw std::length_error("Function B in counter out of range.");

        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({nfunA, nfunB, duplication ? 1u : 0u}),
        true,
        "Coopt of " + std::to_string(nfunA) + " by " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Indicator function. Equals to one if \f$k\f$ genes changed and zero
 * otherwise.
 */
inline void counter_k_genes_changing(
    PhyloCounters * counters,
    unsigned int k,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        if (Array.D()->duplication != (data->at(0u) == 1))
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        unsigned int count = 0u;
        for (uint j0 = 0u; j0 < Array.nrow(); ++j0)
            if (Array.D()->states[j0]) 
                count++;

        return (count == data->at(1u)) ? 1.0: 0.0;
      

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        if (Array.D()->duplication & (data->at(0u) == 0u))
            return 0.0;
        else if (!Array.D()->duplication & (data->at(0u) == 1u))
            return 0.0;

        // Setup
        int count = 0u; ///< How many genes diverge the parent
        int k     = static_cast<int>(data->at(1u));
        std::vector< bool > par_state = Array.D()->states;

        // Comparing the focal gene with the parent. Whas this already
        // different?
        bool j_matches = true;
        for (unsigned int f = 0u; f < Array.nrow(); ++f)
        {
            if (f == i)
                continue;

            if (par_state[f] != Array(f, j))
            {
                j_matches = false;
                break;
            }
        }

        // If there was a function other than this not matching, then
        // there's no change.
        if (!j_matches)
            return 0.0;

        // Otherwise, if the parent doesn't have the function, then it means
        // that this gene now diverges (so we increase the counter).
        if (!par_state[i])
            count++;

        // Iterating through genes
        for (unsigned int g = 0u; g < Array.ncol(); ++g) 
        {

            // We already did j
            if (g == j)
                continue;

            // Iterating through the function
            bool changes = false;
            for (unsigned int f = 0u; f < Array.nrow(); ++f)
            {
                
                if ((Array(f, g) == 1u) != par_state[f])
                {
                    changes = true;
                    break;
                }

            }

            // If the branch has changed
            if (!changes)
                ++count;

        }

        // If it matches, then we are now swithing to 1
        if (count == k) 
            return 1.0;
        // Otherwise, if the difference is only one, then it means
        // that we were matching
        else if (std::fabs(count - k) < 1.0)
            return -1.0;
        else
            return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new PhyloCounterData({duplication ? 1u : 0u, k}),
        true,
        std::to_string(k) + " genes changing" + get_last_name(duplication)
    );
  
}   

///@}

/**
 * @weakgroup rules-phylo Phylo rules
 * @brief Rules for phylogenetic modeling
 * @param rules A pointer to a `PhyloRules` object (`Rules`<`PhyloArray`, `PhyloRuleData`>).
 */
///@{

class PhyloRuleDynData {
public:
    const std::vector< double > * counts;
    uint pos;
    uint lb;
    uint ub;
    bool duplication;
    PhyloRuleDynData(
        const std::vector< double > * counts_,
        uint pos_,
        uint lb_,
        uint ub_,
        bool duplication_
        ) :
        counts(counts_), pos(pos_), lb(lb_), ub(ub_), duplication(duplication_) {};
    
    ~PhyloRuleDynData() {};
};

/**
 * @brief Overall functional gains
 * @param support Support of a model.
 * @param pos Position of the focal statistic.
 * @param lb Lower bound
 * @param ub Upper bound
 * @details 
 * @return (void) adds a rule limiting the support of the model.
 */
inline void rule_dyn_limit_changes(
    PhyloSupport * support,
    uint pos,
    uint lb,
    uint ub,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_RULE_DYN_LAMBDA(tmp_rule)
    {

        if (data->duplication != Array.D()->duplication)
            return true;
        
        if (data->counts->operator[](data->pos) < data->lb)
            return false;
        else if (data->counts->operator[](data->pos) > data->ub)
            return false;
        else
            return true;
      
    };
    
    support->get_rules_dyn()->add_rule(
        tmp_rule,
        new PhyloRuleDynData(
            support->get_current_stats(),
            pos, lb, ub, duplication
            ),
        true 
    );

    return;
  
}

///@}

#undef MAKE_DUPL_VARS
#undef IS_EITHER
#undef IS_DUPLICATION
#undef IS_SPECIATION
#undef IF_MATCHES
#undef IF_NOTMATCHES

#undef DEFAULT_DUPLICATION
#undef DUPL_SPEC
#undef DUPL_DUPL
#undef DUPL_EITH


#endif
