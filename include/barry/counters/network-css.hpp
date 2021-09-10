#ifndef BARRY_CSS_COUNTERS
#define BARRY_CSS_COUNTERS

// n: Net size, 
// s: Start of the i-th network
// e: end of the i-th network
#define CSS_SIZE() \
    uint n = data->indices[0u]; \
    uint s = data->indices[1u]; \
    uint e = data->indices[2u];

// Variables in case that the current cell corresponds to the True
#define CSS_CASE_TRUTH() if ((i < n) && (j < n)) 
#define CSS_TRUE_CELLS() \
    double tji = static_cast<double>(Array(j, i, false)); \
    double pij = static_cast<double>(Array(i + s, j + s, false)); \
    double pji = static_cast<double>(Array(j + s, i + s, false));

// Variables in case that the current cell corresponds to the Perceived
#define CSS_CASE_PERCEIVED() else if (((i >= s) && (i < e)) & ((j >= s) && (j < e)))
#define CSS_PERCEIVED_CELLS() \
    double tji = static_cast<double>(Array(j - s, i - s, false)); \
    double pji = static_cast<double>(Array(j, i, false)); \
    double tij = static_cast<double>(Array(i - s, j - s, false));

// Nothing for else (for now)
#define CSS_CASE_ELSE()

// Checks whether the start and end of the node (perceived network) falls within
// the boundaries of the graph.
#define CSS_CHECK_SIZE_INIT() \
    /* The indices fall within the network */ \
    if ((data->indices.at(0) > Array.ncol()) \
    | (data->indices.at(2) > Array.ncol())) \
        throw std::range_error("The network does not match the prescribed size."); 

#define CSS_CHECK_SIZE() for (uint i = 0u; i < end_.size(); ++i) {\
    if (i == 0u) continue; \
    else if (end_[i] < end_[i-1u]) \
        throw std::logic_error("Endpoints should be specified in order.");}

#define CSS_APPEND(name) std::string name_ = (name);\
    for (uint i = 0u; i < end_.size(); ++i) { \
    std::string tmpname = name_ + " (" + std::to_string(i) + ")";\
    counters->add_counter(tmp_count, tmp_init,\
            new NetCounterData({netsize, i == 0u ? netsize : end_[i-1], end_[i]}, {}),\
            true, tmpname);}


/** @brief Counts errors of commission 
 * @param netsize Size of the reference (true) network 
 * @param end_ Vector indicating one past the ending index of each network. (see details)
 * @details 
 * The `end_` parameter should be of length `N of networks` - 1. It is
 * assumed that the first network ends at `netsize`.
 */
inline void counter_css_partially_false_recip_commi(
        NetCounters * counters,
        uint netsize,
        const std::vector< uint > & end_
    ) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {
        

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            // Checking change stat of the true net
            CSS_TRUE_CELLS()
            return pij * pji * (1.0 - 2.0 * tji);

        } CSS_CASE_PERCEIVED() {

            // Checking change stat of the percieved net
            CSS_PERCEIVED_CELLS()
            return pji * (tij + tji - 2.0 * tij*tji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init) {

        // The reported size doesn't match the true network
        CSS_CHECK_SIZE_INIT()
        
        return 0.0;
        
    };
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Partially false recip (comission)")

    return;
    
}

/** @brief Counts errors of omission */
inline void counter_css_partially_false_recip_omiss(
        NetCounters * counters,
        uint netsize,
        const std::vector< uint > & end_
    ) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {
        
        // Getting the network size
        CSS_SIZE()
        
        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * (pji + pji - 2.0 * pij * pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return tji * tij * (1.0 - 2.0 * pji);

        } CSS_CASE_ELSE()
            return 0.0;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init) {

        // The reported size doesn't match the true network
        CSS_CHECK_SIZE_INIT()
        
        return 0.0;
        
    };
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Partially false recip (omission)")
        
    return;
    
}

/** @brief Counts completely false reciprocity (comission) */
inline void counter_css_completely_false_recip_comiss(
        NetCounters * counters,
        uint netsize,
        const std::vector< uint > & end_
    ) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return -(1.0 - tji) * pij * pji;

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return (1.0 - tij) * (1.0 - tji) * pji;

        } CSS_CASE_ELSE()
            return 0.0;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init) {

        // The reported size doesn't match the true network
        CSS_CHECK_SIZE_INIT()
        
        return 0.0;
        
    };
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Completely false recip (comission)")
        
    return;
    
}

/** @brief Counts completely false reciprocity (omission) */
inline void counter_css_completely_false_recip_omiss(
        NetCounters * counters,
        uint netsize,
        const std::vector< uint > & end_
    ) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * (1.0 - pij) * (1.0 - pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return - tij * tji * (1.0 - pji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init) {

        // The reported size doesn't match the true network
        CSS_CHECK_SIZE_INIT()
        
        return 0.0;
        
    };
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Completely false recip (omission)")
        
    return;
    
}

#undef CSS_CASE_TRUTH
#undef CSS_TRUE_CELLS
#undef CSS_CASE_PERCEIVED
#undef CSS_PERCEIVED_CELLS
#undef CSS_CASE_ELSE
#undef CSS_CHECK_SIZE_INIT
#undef CSS_CHECK_SIZE
#endif
