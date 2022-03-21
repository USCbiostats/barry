#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H 1

// Default value that is used for the counters.
#define DEFAULT_DUPLICATION 1u
#define DUPL_SPEC 0u
#define DUPL_DUPL 1u
#define DUPL_EITH 2u


#define MAKE_DUPL_VARS() \
    bool DPL = Array.D_ptr()->duplication; \
    unsigned int DATA_AT = data[0u];

#define IS_EITHER()      (DATA_AT == DUPL_EITH)
#define IS_DUPLICATION() ((DATA_AT == DUPL_DUPL) & (DPL))
#define IS_SPECIATION()  ((DATA_AT == DUPL_SPEC) & (!DPL))

#define IF_MATCHES() MAKE_DUPL_VARS() \
    if (IS_EITHER() | IS_DUPLICATION() | IS_SPECIATION())
#define IF_NOTMATCHES() MAKE_DUPL_VARS() \
    if (!IS_EITHER() & !IS_DUPLICATION() & !IS_SPECIATION())


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

    PhyloCounterData() : data(0u) {};

    uint at(uint d) {return data.at(d);};
    uint operator()(uint d) {return data.at(d);};
    uint operator[](uint d) {return data[d];};
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
typedef BArrayDense<uint, NodeData> PhyloArray;
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
    [](const PhyloArray & Array, uint i, uint j, PhyloCounterData & data)

#define PHYLO_RULE_DYN_LAMBDA(a) Rule_fun_type<PhyloArray, PhyloRuleDynData> a = \
    [](const PhyloArray & Array, uint i, uint j, PhyloRuleDynData & data)

#define PHYLO_CHECK_MISSING() if (Array.D_ptr() == nullptr) \
    throw std::logic_error("The array data is nullptr."); \
    
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
        IF_NOTMATCHES()
            return 0.0;
      
        return Array.D_ptr()->states[i] ? 0.0 : 1.0;
      
    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication}),
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

        IF_NOTMATCHES()
            return 0.0;

        double ngains = 0.0;
        auto   k = data[1u];
        auto   s = Array.D_ptr()->states[k];

        if (s)
            return 0.0;

        for (auto o = 0u; o < Array.ncol(); ++o)
        {
            if (!s && (Array(k,o) == 1u))
                ngains += 1.0;
        }
        
        return ngains;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is there any gain?
        if (Array.D_ptr()->states[i])
            return 0.0;

        IF_MATCHES()
            return (i == data[1u]) ? 1.0 : 0.0;
        
        return 0.0;

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init,
            PhyloCounterData({duplication, i}),
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
        if (i != data[1u])
            return 0.0;

        IF_NOTMATCHES()
            return 0.0;

        // Is there any gain?
        if (Array.D_ptr()->states[i])
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
        int diff = static_cast<int>(data[2u]) - counts + 1;
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
            PhyloCounterData({duplication, i, k}),
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
        for (auto s : Array.D_ptr()->states)
        {

            if (s) 
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
            if ((k != i) && (Array.D_ptr()->states[k] != (Array(k, j, false) == 1u)))
                return 0.0;
            

        }

        // Nope, this gene is now matching its parent, so we need to 
        // take it out from the count of genes that have changed.
        return Array.D_ptr()->states[i] ? -1.0 : 1.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication}),
        "Num. of genes changing" + get_last_name(duplication)
    );

    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Keeps track of how many pairs of genes preserve pseudostate.
 */
inline void counter_preserve_pseudogene(
    PhyloCounters * counters,
    unsigned int nfunA,
    unsigned int nfunB,
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
        if (Array.D_ptr()->states[data[1u]] | Array.D_ptr()->states[data[2u]])
            return 0.0;

        double n = static_cast<double>(Array.ncol());
        return n * (n - 1.0) / 2.0;
      

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;

        auto nfunA = data[1u];
        auto nfunB = data[2u];

        if ((i != nfunA) & (i != nfunB))
            return 0.0;

        if (Array.D_ptr()->states[data[1u]] | Array.D_ptr()->states[data[2u]])
            return 0.0;

        unsigned int k = (i == nfunA) ? nfunB : nfunA;

        if (Array(k, j) == 1u)
            return 0.0;

        double res = 0.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {
            if (off == j)
                continue;

            if ((Array(i, off) == 0u) && (Array(k, off) == 0u))
                res -= 1.0;

        }

        return res;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication}),
        "Preserve pseudo gene (" + 
        std::to_string(nfunA) + ", " +
        std::to_string(nfunB) + ")" + get_last_name(duplication)
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
        for (auto s : Array.D_ptr()->states)
        {
            if (s)
                return 1.0;
        }
        
        return 0.0;
      
    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;
        
        // Setup
        bool j_diverges = false;
        const std::vector< bool > & par_state = Array.D_ptr()->states;

        for (unsigned int f = 0u; f < Array.nrow(); ++f)
        {

            // Was the gene annotation different from the parent?
            if (par_state[f] != (Array(f,j) == 1u))
            {
                j_diverges = true;
                break;
            }

        }


        bool j_used_to_diverge = false;
        for (unsigned int f = 0u; f < Array.nrow(); ++f)
        {

            if (f == i)
            {
                if (par_state[f])
                {
                    j_used_to_diverge = true;
                    break;
                }
            }
            else
            {

                if (par_state[f] != (Array(f,j) == 1u))
                {
                    j_used_to_diverge = true;
                    break;
                }

            }

        }

        // Case 1: j hasn't changed
        if ((!j_used_to_diverge & !j_diverges) | (j_used_to_diverge & j_diverges))
            return 0.0;
        // Case 2: j NOW diverges
        else if (j_diverges)
            return 1.0/Array.ncol();
        // Case 3: j USED to diverge, so no more
        else
            return -1.0/Array.ncol();

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication}),
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

        if (!Array.D_ptr()->states[i])
            return 0.0;

        IF_MATCHES()
            return -1.0;
        else 
            return 0.0;
        
    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        IF_NOTMATCHES()
            return 0.0;
        
        double res = 0.0;
        for (auto s : Array.D_ptr()->states)
            if (s)
                res += 1.0;

        return res * static_cast<double>(Array.ncol());

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication}),
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

        IF_NOTMATCHES()
            return 0.0;

        // At first, all are zero, so we need to check if the lower
        // bound is zero
        if (data[1u] == 0)
            return static_cast<double>(Array.ncol());
        
        return 0.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        int count = Array.colsum(j);
        int ub    = data[2u];
        
        // It now matches
        if (count == static_cast<int>(data[1u]))
            return 1.0;

        // Was within, but now outside
        if (count > ub && ((count - ub) == 1))
            return -1.0;

        // Otherwise nothing happens.
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication, lb, ub}),
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
    PhyloCounters * counters,
    std::vector<uint> nfun,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        if (!Array.D_ptr()->states[i])
            return 0.0;
        
        return (i == data[1u]) ? -1.0 : 0.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;
        
        auto f = data[1u];

        if (!Array.D_ptr()->states[f])
            return 0.0;
        
        return static_cast<double>(Array.ncol());

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init,
            PhyloCounterData({duplication, i}),
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

        if (Array.D_ptr()->states[i])
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
            if (Array.D_ptr()->states[k])
                counts += noff;

        return counts;

        

    };

    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication}),
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

        auto funA = data[1u];
        auto funB = data[2u];
        
        // Are we looking at either of the relevant functions?
        if ((funA != i) && (funB != i))
            return 0.0;
        
        // Are A and B existant? if not, no change
        if (!Array.D_ptr()->states[funA] | !Array.D_ptr()->states[funB])
            return 0.0;
        
        // Figuring out which is the first (reference) function
        uint other = (i == funA)? funB : funA;
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
        PhyloCounterData({duplication, nfunA, nfunB}),
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

        auto d1 = data[1u];
        auto d2 = data[2u];
      
        // Is the function in scope relevant?
        if ((i != d1) && (i != d2))
            return 0.0;
        
        // None should have it
        if (!Array.D_ptr()->states[d1] && !Array.D_ptr()->states[d2])
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
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Co-gains " + std::to_string(nfunA) + " & " + std::to_string(nfunB) +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/** @brief Longest branch mutates (either by gain or by loss) */
inline void counter_longest(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
    )
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // Figuring out which match
        std::vector< bool> is_longest(Array.ncol(), false);
        bool j_mutates = false;
        int nmutate = 0;
        int nmutate_longest = 0;

        auto states  = Array.D_ptr()->states;
        
        for (auto off = 0u; off < Array.ncol(); ++off)
        {

            // On the fly, figuring out if it is longest
            for (auto & l : data)
                if (l == off)
                    is_longest[off] = true;

            for (auto f = 0u; f < Array.nrow(); ++f)
            {
                if ((Array(f, off) == 1u) != states[f])
                {
                    
                    // If it happens that j != off and is not longest
                    // then return 0 (a not longest was mutating prev)
                    if (is_longest[off] && (off != j))
                        return 0.0;

                    if (off == j)
                        j_mutates = true;

                    if (is_longest[j])
                        nmutate_longest++;
                    else
                        nmutate++;

                    break;
                }

            }
        }

        // There was already more than one in difference
        // so nothing to change
        if (std::fabs(nmutate - nmutate_longest) > 1)
            return 0.0;

        // Figuring out previously
        bool j_mutates_prev = false;
        for (auto f = 0u; f < Array.nrow(); ++f)
        {
            // Checking the previous function... was it
            // different before?
            if ((f == i) && states[i])
            {
                j_mutates_prev = true;
                break;
            }
            else if ((Array(f, j) == 1u) != states[f])
            {
                j_mutates_prev = true;
                break;
            }

        }

        // Adjusting the previous count
        auto nmutate_prev         = nmutate;
        auto nmutate_longest_prev = nmutate_longest;
        if (j_mutates & !j_mutates_prev)
        {
            if (is_longest[j])
                nmutate_longest_prev--;
            else
                nmutate_prev--;
        }
        else if (!j_mutates & j_mutates)
        {
            if (is_longest[j])
                nmutate_longest_prev++;
            else
                nmutate_prev++;

        }
        
        // Just compute the change statistic directly
        return
            ( ((nmutate == 0) & (nmutate_longest > 0)) ? 1.0 : 0.0 ) +
            ( ((nmutate_prev == 0) & (nmutate_longest_prev > 0)) ? 1.0 : 0.0 );

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        
        if (Array.D_ptr()->blengths.size() != Array.ncol())
            throw std::logic_error(
                "longest should be initialized with a vec of size Array.ncol()."
            );
          
        // Finding the longest branch (or branches) --
        uint longest_idx = 0u;
        double diff      = 0.0;
        data.reserve(Array.ncol()); 
        data.push_back(0u);
        for (uint ii = 1u; ii < Array.ncol(); ++ii)
        {
            
            diff = Array.D_ptr()->blengths[longest_idx] - Array.D_ptr()->blengths[ii];
            if (diff > 0.0)
                continue;
            else if (diff < 0.0)
            {

                data.empty();
                data.push_back(ii);
                longest_idx = ii;

            }
            else if (diff == 0.0)
                data.push_back(ii);
            
        }

        data.shrink_to_fit();
        
        if (data.size() == 0u)
            throw std::logic_error("The data on the longest branch has size 0.");
        
        // Starting the counter, since all in zero, then this will be equal to
        // the number of functions in 1 x number of longest branches
        for (uint ii = 0u; ii < Array.nrow(); ++ii)
        {
            
            if (Array.D_ptr()->states[ii])
                return (1.0 * static_cast<double>(data.size()));

        }
        
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication}),
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
        
        auto funA = data[1u];
        auto funB = data[2u];

        // Is the function in scope relevant?
        if ((i != funA) && (i != funB))
            return 0.0;
        
        // Checking if the parent has both functions
        uint other = (i == funA)? funB : funA;
        bool parent_i     = Array.D_ptr()->states[i];
        bool parent_other = Array.D_ptr()->states[other];
        
        if (!parent_i & !parent_other) 
            return 0.0;
        else if (parent_i & parent_other) 
            return 0.0;
        
        // Figuring out which is the first (reference) function
        double res = 0.0;
        
        if (Array(other, j) == 0u)
        {


            for (auto off = 0u; off < Array.ncol(); ++off)
                if ((Array(i,off) == 0) && (Array(other,off) == 1))
                    res += 1.0;

        }
        else
        {

            for (auto off = 0u; off < Array.ncol(); ++off)
                if ((Array(i,off) == 1) && (Array(other,off) == 0))
                    res -= 1.0;
                
        }
             
        return res;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Neofun between " + std::to_string(nfunA) + " and " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
  
}

//------------------------------------------------------------------------------
/**
 * @brief Total number of neofunctionalization events 
 * sum_u sum_{w < u} [x(u,a)*(1 - x(w,a)) + (1 - x(u,a)) * x(w,a)]
 * change stat: delta{x(u,a): 0->1} = 1 - 2 * x(w,a)
 */
inline void counter_pairwise_neofun_singlefun(
    PhyloCounters * counters,
    uint nfunA,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;
        
        // Is the function in scope relevant?
        if (i != data[1u])
            return 0.0;
        
        // Checking if the parent has the function
        if (Array.D_ptr()->states[i])
            return 0.0;
        
        // Figuring out which is the first (reference) function
        double res = 0.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {

            if (off == j)
                continue;

            if ((Array(i, off) == 0))
                res += 1.0;
            else 
                res -= 1.0;

        }

        return res;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication, nfunA}),
        "Pairwise neofun function " + std::to_string(nfunA) +
        get_last_name(duplication)
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
        IF_NOTMATCHES()
            return 0.0;
        
        const uint & funA = data[1u];
        const uint & funB = data[2u];

        // Checking scope
        if ((i != funA) && (i != funB))
            return 0.0;
        
        // Checking the parent doesn't have funA or has funB
        if (!Array.D_ptr()->states[funA] | Array.D_ptr()->states[funB]) 
            return 0.0;

        double res = 0.0;

        if (i == funA)
        {

            if (Array(funB, j) == 0u)
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {

                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 0u) && (Array(funB, off) == 1u))
                        res += 1.0;

                }

            }
            else
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {
                    
                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 1u) && (Array(funB, off) == 0u))
                        res -= 1.0;

                }

            }

        }
        else
        {

            if (Array(funA, j) == 0u)
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {

                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 1u) && (Array(funB, off) == 0u))
                        res += 1.0;

                }

            }
            else
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {
                    
                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 0u) && (Array(funB, off) == 1u))
                        res -= 1.0;

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
        PhyloCounterData({duplication, nfunA, nfunB}),
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
        IF_NOTMATCHES()
            return 0.0;
        
        const unsigned int funA = data[1u];
        const unsigned int funB = data[2u];

        // If the change is out of scope, then nothing to do
        if ((i != funA) & (i != funB))
            return 0.0;

        // If the parent does not have the initial state, then it makes no sense
        if ((!Array.D_ptr()->states[funA]) | Array.D_ptr()->states[funB])
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

        

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        if (data.size() != 3u)
            throw std::length_error("The counter data should be of length 2.");

        if (data[1u] == data[2u])
            throw std::logic_error("Functions A and B should be different from each other.");

        if (data[1u] >= Array.nrow())
            throw std::length_error("Function A in counter out of range.");

        if (data[2u] >= Array.nrow())
            throw std::length_error("Function B in counter out of range.");

        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication, nfunA, nfunB}),
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

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        for (auto s : Array.D_ptr()->states)
            if (s)
                return Array.ncol() == data[1u] ? 1.0 : 0.0;

        return data[1u] == 0 ? 1.0 : 0.0;
      
    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;
        
        // How many genes diverge the parent
        int              count = 0; 
        bool        j_diverges = false;
        const auto & par_state = Array.D_ptr()->states;

        int k = static_cast<int>(data[1u]);

        for (auto o = 0u; o < Array.ncol(); ++o)
        {

            for (auto f = 0u; f < Array.nrow(); ++f)
            {

                // Was the gene annotation different from the parent?
                if ((Array(f, o) == 1u) != par_state[f])
                {

                    if (o == j)
                        j_diverges = true;

                    count++;
                    break;

                }

            }

        }

        // Counts will only be relevant if (count - k) > 1. Otherwise,
        // having the j gene changed is not relevant
        if (std::abs(count - k) > 1)
            return 0.0;

        // Did it used to diverge?
        bool j_used_to_diverge = false;
        for (auto f = 0u; f < Array.nrow(); ++f)
        {

            if (f == i)
            {
                if (par_state[f]) // Since it is now true, it used to diverge
                {
                    j_used_to_diverge = true;
                    break;
                }
            }
            else
            {

                if (par_state[f] != (Array(f,j) == 1u))
                {
                    j_used_to_diverge = true;
                    break;
                }

            }

        }

        auto count_prev = count;
        // Case 1: j hasn't changed
        if ((!j_used_to_diverge & !j_diverges) | (j_used_to_diverge & j_diverges))
            return 0.0;
        // Case 2: j NOW diverges
        else if (j_diverges)
            count_prev--;
        // Case 3: j USED to diverge
        else
            count_prev++;

        return (count == k ? 1.0 : 0.0) - (count_prev == k ? 1.0 : 0.0);

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication, k}),
        std::to_string(k) + " genes changing" + get_last_name(duplication)
    );
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Indicator function. Equals to one if \f$k\f$ genes changed and zero
 * otherwise.
 */
inline void counter_less_than_p_prop_genes_changing(
    PhyloCounters * counters,
    double p,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        for (auto s : Array.D_ptr()->states)
            if (s)
                return data[1u] == 100 ? 1.0 : 0.0;

        // Only one if it was specified it was zero
        return 1.0;
      
    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;
        
        // Setup
        double count = 0.0; ///< How many genes diverge the parent

        bool j_diverges = false;
        const std::vector< bool > & par_state = Array.D_ptr()->states;

        for (unsigned int o = 0u; o < Array.ncol(); ++o)
        {

            for (unsigned int f = 0u; f < Array.nrow(); ++f)
            {

                // Was the gene annotation different from the parent?
                if ((Array(f, o) == 1u) != par_state[f])
                {

                    if (o == j)
                        j_diverges = true;

                    count += 1.0;
                    break;

                }

            }

        }


        bool j_used_to_diverge = false;
        for (unsigned int f = 0u; f < Array.nrow(); ++f)
        {

            if (f == i)
            {
                if (par_state[f])
                {
                    j_used_to_diverge = true;
                    break;
                }
            }
            else
            {

                if (par_state[f] != (Array(f,j) == 1u))
                {
                    j_used_to_diverge = true;
                    break;
                }

            }

        }

        auto count_prev = count;
        // Case 1: j hasn't changed
        if ((!j_used_to_diverge & !j_diverges) | (j_used_to_diverge & j_diverges))
            return 0.0;
        // Case 2: j NOW diverges
        else if (j_diverges)
            count_prev -= 1.0;
        // Case 3: j USED to diverge
        else
            count_prev += 1.0;

        double ncol = static_cast<double>(Array.ncol());
        double p    = static_cast<double>(data[1u]) / 100.0;

        return ((count/ncol) <= p ? 1.0 : 0.0) - ((count_prev/ncol) <= p ? 1.0 : 0.0);

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication, static_cast<uint>(p * 100)}),
        std::to_string(p) + " prop genes changing" + get_last_name(duplication)
    );
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 */
inline void counter_gains_from_0(
    PhyloCounters * counters,
    std::vector< uint > nfun,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // All must be false
        for (auto s : Array.D_ptr()->states)
        {

            if (s)
                return 0.0;

        }

        // Is this the function?
        if (i != data[1u])
            return 0.0;

        // Now computing the change stats
        double res = static_cast<double>(Array.ncol()) - 1.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {
            if (off  == j)
                continue;

            if (Array(i, off) == 1u)
                res -= 2.0;
        }


        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init,
            PhyloCounterData({duplication, i}),
            "First gain " + std::to_string(i) +
                get_last_name(duplication)
        );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 */
inline void counter_overall_gains_from_0(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // All must be false
        for (auto s : Array.D_ptr()->states)
        {

            if (s)
                return 0.0;

        }

        return 1.0;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication}),
        "Overall first gains" +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 */
inline void counter_pairwise_overall_change(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        unsigned int funpar = Array.D_ptr()->states[i] == 1u;

        // All must be false
        double res = 0.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {
            if (off == j)
                continue;

            if (funpar > Array(i, off))
                res -= 1.0;
            else if (funpar < Array(i, off))
                res += 1.0;
        }
        
        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        double res = 0.0;
        double n   = static_cast<double>(Array.ncol());
        for (auto s : Array.D_ptr()->states)
            if (s)
                res += n * (n - 1.0) / 2.0;

        return res;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication}),
        "Pairs of genes changing" +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 * sum x(a)^3(1-x(b))^3 + x(b)^3(1-x(a))^3 + x(a)^3 * x(b)^3 + (1 - x(a))^3 * (1-x(b))^3
 */
inline void counter_pairwise_preserving(
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

        // Not in the scope
        auto funA = data[1u];
        auto funB = data[2u];
        if ((funA != i) && (funB != i))
            return 0.0;

        unsigned int k = (funA == i) ? funB : funA;

        bool parent_i = Array.D_ptr()->states[i];
        bool parent_k = Array.D_ptr()->states[k];

        // if (!parent_i & !parent_k)
        //     return 0.0;

        double res = 0.0;
        // Case 1: (0,0)
        if (!parent_i & !parent_k)
        {

            if (Array(k, j) == 1u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 0u) && (Array(k, off) == 0u))
                    res -= 1.0;

            }

        }
        else if (parent_i & !parent_k)
        {

            if (Array(k, j) == 1u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 1u) && (Array(k, off) == 0u))
                    res += 1.0;

            }

        }
        else if (!parent_i & parent_k)
        {

            if (Array(k, j) == 0u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 0u) && (Array(k, off) == 1u))
                    res += 1.0;

            }

        }
        else
        {

            if (Array(k, j) == 0u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 1u) && (Array(k, off) == 1u))
                    res += 1.0;

            }
        }

        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {


        IF_NOTMATCHES()
            return 0.0;
        
        PHYLO_CHECK_MISSING();
        
        double n = static_cast< double >(Array.ncol());
        if (!Array.D_ptr()->states[data[1u]] && !Array.D_ptr()->states[data[2u]])
            return n * (n - 1.0) / 2.0;

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Pariwise preserve (" + std::to_string(nfunA) + ", " +
            std::to_string(nfunB) + ")" +get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 * sum x(a)^3(1-x(b))^3 + x(b)^3(1-x(a))^3 + x(a)^3 * x(b)^3 + (1 - x(a))^3 * (1-x(b))^3
 */
inline void counter_pairwise_first_gain(
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

        // Not in the scope
        auto funA = data[1u];
        auto funB = data[2u];
        if ((funA != i) && (funB != i))
            return 0.0;

        unsigned int k = (funA == i) ? funB : funA;

        double res = 0.0;
        if (Array(k, j) == 1)
        {

            for (auto off = 0u; off < Array.ncol(); ++off)
            {
                if (off == j)
                    continue;

                if ((Array(i,off) == 0u) && (Array(k,off) == 0u))
                    res -= 1.0;
            }

        }
        else
        {

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 1u))
                {

                    // j: (0,0)\(1,0) -> (1,0)\(1,0), so less 1
                    if (Array(k, off) == 0u)
                        res -= 1.0;

                }
                else
                {

                    if (Array(k, off) == 1u) 
                    // j: (0,0)\(0,1) -> (1,0)\(0,1), so less 1
                        res -= 1.0;
                    else
                    // j: (0,0)\(0,0) -> (1,0)\(0,0), so plus 1
                        res += 1.0;

                }

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
        PhyloCounterData({duplication, nfunA, nfunB}),
        "First gain (either " + std::to_string(nfunA) + " or " +
            std::to_string(nfunB) + ")" +get_last_name(duplication)
    );
    
    return;
  
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
    uint duplication;

    PhyloRuleDynData(
        const std::vector< double > * counts_,
        uint pos_,
        uint lb_,
        uint ub_,
        uint duplication_
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

        unsigned int rule_type = data.duplication;
        if (rule_type != DUPL_EITH)
        {

            if (Array.D_ptr()->duplication & (rule_type != DUPL_DUPL))
                return true;
            else if (!Array.D_ptr()->duplication & (rule_type != DUPL_SPEC))
                return true;
                
        }

        if (data.counts->operator[](data.pos) < data.lb)
            return false;
        else if (data.counts->operator[](data.pos) > data.ub)
            return false;
        else
            return true;
      
    };
    
    support->get_rules_dyn()->add_rule(
        tmp_rule,
        PhyloRuleDynData(
            support->get_current_stats(),
            pos, lb, ub, duplication
            )
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
