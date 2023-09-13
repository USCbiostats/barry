#ifndef BARRY_CSS_COUNTERS
#define BARRY_CSS_COUNTERS

// n: Net size, 
// s: Start of the i-th network
// e: end of the i-th network
// ego_id: Ego of the cell (i, j)
#define CSS_SIZE() \
    size_t n      = data.indices[0u]; \
    size_t s      = data.indices[1u]; \
    size_t e      = data.indices[2u]; \
    size_t ctype  = data.indices[3u]; \
    size_t ego_id = data.indices[4u]; \
    if (ctype > 2) \
        throw std::range_error("Counter type should be 0, 1, or 2.");

// Check whether ego_id is involved in the current cell
// ctype: Type of counter
// 0: All cells
// 1: Only if perceiver
// 2: Only if not perceiver
#define CSS_MATCH_TYPE() \
    if (ctype != 0u) { /* all counts */ \
        if (ctype == 1u) { /* Only if perceiver */ \
            if ((i_ != ego_id) && (j_ != ego_id)) return 0.0; \
        } else if (ctype == 2u) { /* Only if not perceiver */ \
            if ((i_ == ego_id) || (j_ == ego_id)) return 0.0; \
        } \
    };

// Variables in case that the current cell corresponds to the True
#define CSS_CASE_TRUTH() if ((i < n) && (j < n)) 

// i_: i-th index of the cell
// j_: j-th index of the cell
// tji: True value of the cell (i, j)
// pij: Perceived value of the cell (i, j)
// pji: Perceived value of the cell (j, i)
#define CSS_TRUE_CELLS() \
    size_t i_ = i; \
    size_t j_ = j; \
    CSS_MATCH_TYPE() \
    double tji = static_cast<double>(Array(j, i, false)); \
    double pij = static_cast<double>(Array(i + s, j + s, false)); \
    double pji = static_cast<double>(Array(j + s, i + s, false));

// Variables in case that the current cell corresponds to the Perceived
#define CSS_CASE_PERCEIVED() else if (((i >= s) && (i < e)) & ((j >= s) && (j < e)))

// i_: i-th index of the cell
// j_: j-th index of the cell
// tji: True value of the cell (i, j)
// pji: Perceived value of the cell (i, j)
// tij: True value of the cell (j, i)
#define CSS_PERCEIVED_CELLS() \
    size_t i_ = i - s; \
    size_t j_ = j - s; \
    CSS_MATCH_TYPE() \
    double tji = static_cast<double>(Array(j - s, i - s, false)); \
    double pji = static_cast<double>(Array(j, i, false)); \
    double tij = static_cast<double>(Array(i - s, j - s, false));



// Nothing for else (for now)
#define CSS_CASE_ELSE()

// Checks whether the start and end of the node (perceived network) falls within
// the boundaries of the graph.
#define CSS_CHECK_SIZE_INIT() \
    /* The indices fall within the network */ \
    if ((data.indices.at(0) > Array.ncol()) \
    | (data.indices.at(2) > Array.ncol())) \
        throw std::range_error("The network does not match the prescribed size."); 

#define CSS_CHECK_SIZE() for (size_t i = 0u; i < end_.size(); ++i) {\
    if (i == 0u) continue; \
    else if (end_[i] < end_[i-1u]) \
        throw std::logic_error("Endpoints should be specified in order.");}

#define CSS_APPEND(name) std::string name_ = (name);\
    for (size_t i = 0u; i < end_.size(); ++i) { \
    std::string tmpname = name_ + " (" + std::to_string(i) + ")" + \
    ((counter_type == 1u) ? " (only perceiver)" : ((counter_type == 2u)? " (only alters)": ""));\
    counters->add_counter(tmp_count, tmp_init, nullptr, \
            NetCounterData({netsize, i == 0u ? netsize : end_[i-1], end_[i], counter_type, i}, {}),\
            tmpname);}

#define CSS_NET_COUNTER_LAMBDA_INIT() NETWORK_COUNTER_LAMBDA(tmp_init) {\
        CSS_CHECK_SIZE_INIT() \
        return 0.0; \
    };


/**
 * @brief Counts errors of commission 
 * @param netsize Size of the reference (true) network 
 * @param end_ Vector indicating one past the ending index of each network. (see details)
 * @param counter_type Size_t indicating the type of counter to use. Possible
 *  values are: 0: Count all, 1: Only count if perceiver is involved, and 
 *  2: Only count if perceiver is not involved.
 * @details 
 * The `end_` parameter should be of length `N of networks` - 1. It is
 * assumed that the first network ends at `netsize`.
 */
template<typename Tnet = Network>
inline void counter_css_partially_false_recip_commi(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {
        

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            // Checking change stat of the true net
            CSS_TRUE_CELLS()
            return pij * pji * (1.0 - 2.0 * tji) - (1.0 - tji)*(
                pij * (1.0 - pji) + (1.0 - pij) * pji
            );

        } CSS_CASE_PERCEIVED() {

            // Checking change stat of the percieved net
            CSS_PERCEIVED_CELLS()
            return pji * (tij * (1.0 - tji) + (1.0 - tij) * tji) +
                (1.0 - tij) * (1.0 - tji) * (1 - 2.0 * pji)
            ;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };

    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Partially false recip (comission)")

    return;
    
}

/** @brief Counts errors of omission */
template<typename Tnet = Network>
inline void counter_css_partially_false_recip_omiss(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {
        
        // Getting the network size
        CSS_SIZE()
        
        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * ((1.0 - pij) * pji + pij * (1.0 - pji)) +
                (1.0 - 2.0 * tji) * (1.0 - pij) * (1.0 - pji)
            ;

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return tji * tij * (1.0 - 2.0 * pji) - 
                (1.0 - pji) * ((1.0 - tij) * tji + tij * (1.0 - tji))
            ;

        } CSS_CASE_ELSE()
            return 0.0;

    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Partially false recip (omission)")
        
    return;
    
}

/** @brief Counts completely false reciprocity (comission) */
template<typename Tnet = Network>
inline void counter_css_completely_false_recip_comiss(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
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
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Completely false recip (comission)")
        
    return;
    
}

/** @brief Counts completely false reciprocity (omission) */
template<typename Tnet = Network>
inline void counter_css_completely_false_recip_omiss(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
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
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Completely false recip (omission)")
        
    return;
    
}

/** @brief Counts mixed reciprocity errors */
template<typename Tnet = Network>
inline void counter_css_mixed_recip(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return (1.0 - tji) * (1.0 - pij) * pji - tji * pij * (1.0 - pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return (1.0 - tij) * tji * (1.0 - pji) - tij * (1.0 - tij) * pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Mixed reciprocity errors")
        
    return;
    
}

/////////////////////////// CENSUS

template<typename Tnet = Network>
inline void counter_css_census01(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return -(1.0 - pij) * (1.0 - pji) * (1.0 - tji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return -(1.0 - tij) * (1.0 - tji) * (1.0 - pji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    // CSS_NET_COUNTER_LAMBDA_INIT()
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {

        CSS_CHECK_SIZE_INIT()
        double n_dbl = static_cast<double>(data.indices[0u]);

        // Discount in case of the type of counter
        size_t ctype = data.indices[3u];

        if (ctype == 1u) /* Only perceiver */
        {

            return (n_dbl - 1.0) * (Array.D().directed ? 2.0 : 1.0);

        } else if (ctype == 2u) /* All but the perceiver */
        {
            // We remove the perceiver from the eq.
            n_dbl -= 1.0;
        }

        // At the beginning is all zero
        return n_dbl * (n_dbl - 1.0)/ (Array.D().directed ? 1.0 : 2.0);

    };
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(01) Accurate null")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census02(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return -(1.0 - tji) * ((1.0 - pij) * pji + pij * (1.0 - pji));

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return (1.0 - tij) * (1.0 - tji) * (1 - 2.0 * pji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(02) Partial false positive (null)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census03(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
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
            return (1.0 - tij) * (1.0 - tji) *pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(03) Complete false positive (null)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census04(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return (1.0 - pij) * (1.0 - pji) * (1.0 - 2.0 * tji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return -(1.0 - pji) * ((1.0 - tij) * tji + tij * (1.0 - tji));

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(04) Partial false negative (assym)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census05(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return pij * (1.0 - tji) * (1.0 - pji) - (1.0 - pij) * tji * pji;

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return tij * (1.0 - tji) * (1.0 - pji) - (1.0 - tij) * tji * pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(05) Accurate assym")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census06(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return (1.0 - pij) * (1.0 - tji) * pji - pij * tji * (1.0 - pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return (1.0 - tij) * tji * (1.0 - pji) - tij * (1.0 - tji) * pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(06) Mixed assym")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census07(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return pij * pji * (1.0 - 2.0 * tji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return pji * (tij * (1.0 - tji) + (1.0 - tij) * tji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(07) Partial false positive (assym)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census08(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
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
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(08) Complete false negative (full)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census09(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * (pij * (1.0 - pji) + (1.0 - pij) * pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return tij * tji * (1.0 - 2.0 * pji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(09) Partial false negative (full)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census10(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * pij * pji;

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return tij * tji * pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(10) Accurate full")
        
    return;

}

#undef CSS_APPEND
#undef CSS_CASE_TRUTH
#undef CSS_TRUE_CELLS
#undef CSS_CASE_PERCEIVED
#undef CSS_PERCEIVED_CELLS
#undef CSS_CASE_ELSE
#undef CSS_CHECK_SIZE_INIT
#undef CSS_CHECK_SIZE
#undef CSS_NET_COUNTER_LAMBDA_INIT
#undef CSS_MATCH_TYPE
#undef CSS_SIZE
#endif
