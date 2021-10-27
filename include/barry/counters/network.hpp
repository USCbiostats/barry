#ifndef BARRAY_NETWORK_H
#define BARRAY_NETWORK_H 1

/**
 * @ingroup counting 
 * @details Details on the available counters for `NetworkData` can be found in
 * the \ref counters-network section.
 * 
 */
///@{

/**
 * @brief Data class for Networks.
 * 
 * This holds information about whether the graph is directed or not, and,
 * if defined, vectors of node (vertex) attributes (`vertex_attr`).
 * 
 */
class NetworkData {
public:
    
    bool directed = true;
    std::vector< std::vector< double > > vertex_attr;
    
    NetworkData() : vertex_attr(0u) {};
    
    /**
     * @brief Constructor using a single attribute
     * @param vertex_attr_ Double vector of length equal to the number of vertices
     * in the data.
     * @param directed_ When `true` the graph as treated as directed.
     */
    NetworkData(
        std::vector< double >  vertex_attr_,
        bool directed_ = true
    ) : directed(directed_), vertex_attr(1u, vertex_attr_) {};
    
    /**
     * @brief Constructor using multiple attributes
     * @param vertex_attr_ Vector of double vectors. The size equals to the number
     * of attributes to be created. Each individual vector should be of length
     * equal to the number of vertices.
     * @param directed_ When `true` the graph as treated as directed.
     */
    NetworkData(
        std::vector< std::vector< double > > vertex_attr_,
        bool directed_ = true
    ) : directed(directed_), vertex_attr(vertex_attr_) {};
    
    
    ~NetworkData() {};
};

/**
  * @brief Data class used to store arbitrary uint or double vectors */
class NetCounterData {
public:
    
    std::vector< uint > indices;
    std::vector< double > numbers;
    
    NetCounterData() : indices(0u), numbers(0u) {};
    NetCounterData(
        const std::vector< uint > indices_,
        const std::vector< double > numbers_
    ): indices(indices_), numbers(numbers_) {};
    
    ~NetCounterData() {};
    
    // const uint get_uint
    
};

#define NET_C_DATA_IDX(i) (data->indices[i])
#define NET_C_DATA_NUM(i) (data->numbers[i])


/**
 * @name Convenient typedefs for network objects.
 */
///@{
typedef BArray<double, NetworkData> Network;
typedef BArrayDense<double, NetworkData> NetworkDense;
template <typename Tnet>
using NetCounter =  Counter<Network, NetCounterData >;
template <typename Tnet>
using NetCounters =  Counters< Network, NetCounterData>;
template <typename Tnet>
using NetSupport =  Support<Network, NetCounterData >;
template <typename Tnet>
using NetStatsCounter =  StatsCounter<Network, NetCounterData>;
template <typename Tnet>
using NetModel =  Model<Network, NetCounterData>;
template <typename Tnet>
using NetRule =  Rule<Network,bool>;
template <typename Tnet>
using NetRules =  Rules<Network,bool>;
///@}

/**@name Macros for defining counters
  */
///@{
/**Function for definition of a network counter function*/
#define NETWORK_COUNTER(a) inline double (a) \
(const Network & Array, uint i, uint j, NetCounterData * data)
/**Lambda function for definition of a network counter function*/
#define NETWORK_COUNTER_LAMBDA(a) Counter_fun_type<Network, NetCounterData> a = \
    [](const Network & Array, uint i, uint j, NetCounterData * data)
///@}


/**@name Macros for defining rules
  */
///@{
/**Function for definition of a network counter function*/
#define NETWORK_RULE(a) inline bool (a) \
(const Network & Array, uint i, uint j, bool * data)
/**Lambda function for definition of a network counter function*/
#define NETWORK_RULE_LAMBDA(a) Rule_fun_type<Network, bool> a = \
[](const Network & Array, uint i, uint j, bool * data)
///@}

/**
  * @weakgroup  counters-network Network counters
  * @brief Counters for network models
  * @param counters A pointer to a `NetCounters` object (`Counters`<`Network`, `NetCounterData`>).
  */
///@{
// -----------------------------------------------------------------------------
/**@brief Number of edges */
inline void counter_edges(NetCounters * counters)
{
    
    NETWORK_COUNTER_LAMBDA(count_edges)
    {
        return 1.0;
    };
    
    counters->add_counter(
        count_edges, nullptr, nullptr, false, 
        "Edge counts", 
        "Number of edges"
        );
    
    return;

}


// -----------------------------------------------------------------------------
/**@brief Number of isolated vertices */
inline void counter_isolates(NetCounters * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        if (i == j)
            return 0.0;
        
        double res = 0.0;
        
        // i is sending its first tie
        if (Array.row(i).size() == 1u && Array.col(i).size() == 0u)
            res -= 1.0;
        
        // j is receiving its first tie, meaning that he
        // has no other tie but i's?
        if (Array.row(j).size() == 0u && Array.col(j).size() == 1u)
            res -= 1.0;
        
        return res;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        return static_cast<double>(Array.nrow());
    };
    
    counters->add_counter(
        tmp_count,
        tmp_init, nullptr, false, "Isolates", "Number of isolate vertices");

    return;
}

// -----------------------------------------------------------------------------
/**@brief Number of mutual ties */
inline void counter_mutual(NetCounters * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        // Is there any tie at ji? If not, then we have a new mutual!
        // but this only makes sence if the jth row and ith column exists
        // if ((Array.nrow() > j) && (Array.ncol() > i)) 
        if (i == j)
            return 0.0;
        
        // printf_barry("Checking if it is empty or not at (%i, %i)... ", i, j);
        if (!Array.is_empty(j, i, false))
        {
            // printf_barry("Yes, mutual.\n");
            return 1.0;
        }
        // printf_barry("No, no mutual.\n");
        
        return 0.0;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {

        if (Array.nrow() != Array.ncol())
            throw std::logic_error("The -mutual- counter only works on square arrays.");
        
        if (Array.D() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!Array.D()->directed)
            throw std::logic_error(
                "The -mutual- counter only works on directed (non-symmetric) arrays."
                );
        
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr, false, "Reciprocity",
        "Number of mutual ties"
        );

    return;

}


// 2-istars --------------------------------------------------------------------
inline void counter_istar2(NetCounters * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        // Need to check the receiving, if he/she is getting a new set of stars
        // when looking at triads
        
        if (Array.col(j).size() == 1u)
            return 0.0;
        
        return static_cast<double>(Array.col(j).size() - 1.0);

    };
    
    counters->add_counter(tmp_count, nullptr, nullptr, false, "Istar 2", "Indegree 2-star");
    
    return ;
}

// 2-ostars --------------------------------------------------------------------
inline void counter_ostar2(NetCounters * counters)
{
   
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        // Need to check the receiving, if he/she is getting a new set of stars
        // when looking at triads
        
        if (Array.row(i).size() == 1u)
            return 0.0;
        
        return static_cast<double>( Array.row(i).size() - 1.0);

    };
    
    counters->add_counter(tmp_count, nullptr, nullptr, false, "Ostar 2", "Outdegree 2-star");

    return ;
    
}


// ttriads ---------------------------------------------------------------------
inline void counter_ttriads(NetCounters * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        // Self ties do not count
        if (i == j)
            return 0.0;
        
        double ans = 0.0;
        
        // Case 1: i-j, i-k, j-k
        if (Array.row(j).size() < Array.row(i).size())
        {
            
            for (auto j_row = Array.row(j).begin(); j_row != Array.row(j).end(); ++j_row) 
                if ((j != j_row->first) && (i != j_row->first) && !Array.is_empty(i, j_row->first, false))
                    ans += 1.0;
                
        } else {
            
            for (auto i_row = Array.row(i).begin(); i_row != Array.row(i).end(); ++i_row) 
                if ((i != i_row->first) && (i_row->first != j) && !Array.is_empty(j, i_row->first, false))
                    ans += 1.0;
                
        }
        
        // Case 2: i-j, i-k, k-j  
        if (Array.row(i).size() > Array.col(j).size())
        {
            
            for (auto j_col = Array.col(j).begin(); j_col != Array.col(j).end(); ++j_col)
                if ((j != j_col->first) && (i != j_col->first) && !Array.is_empty(i, j_col->first, false))
                    ans += 1.0;
                
        } else {
            
            for (auto i_row = Array.row(i).begin(); i_row != Array.row(i).end(); ++i_row) 
                if ((i != i_row->first) && (j != i_row->first) && !Array.is_empty(i_row->first, j, false))
                    ans += 1.0;
                
        }
        
        // Case 3: 
        if (Array.col(i).size() > Array.col(j).size())
        {
            
            for (auto j_col = Array.col(j).begin(); j_col != Array.col(j).end(); ++j_col)
                if ((j != j_col->first) && (i != j_col->first) && !Array.is_empty(j_col->first, i, false))
                    ans += 1.0;
                
        } else {
            
            for (auto i_col = Array.col(i).begin(); i_col != Array.col(i).end(); ++i_col) 
                if ((i != i_col->first) && (j != i_col->first) && !Array.is_empty(i_col->first, j, false))
                    ans += 1.0;
                
        }
        
        // The regular counter double counts
        return ans;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!(Array.D()->directed))
            throw std::invalid_argument("The ttriads counter is only valid for directed networks. This is undirected.");

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr, false, "Balance",
        "Number of directed triangles"
    );
    
    return;

}


// Cycle triads --------------------------------------------------------------
inline void counter_ctriads(NetCounters * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        if (i == j)
            return 0.0;
        
        double ans = 0.0;
        if (Array.col(i).size() < Array.row(j).size())
        {
            
            for (auto i_col = Array.col(i).begin(); i_col != Array.col(i).end(); ++i_col) 
                if ((i != i_col->first) && (j != i_col->first) && !Array.is_empty(j, i_col->first, false))
                    ans += 1.0;
                
        } else {
            
            for (auto j_row = Array.row(j).begin(); j_row != Array.row(j).end(); ++j_row) 
                if ((j != j_row->first) && (i != j_row->first) && !Array.is_empty(j_row->first, i, false))
                    ans += 1.0;
                
        }
        
        return ans;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {

        if (Array.D() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!(Array.D()->directed))
            throw std::invalid_argument(
                "The ctriads counter is only valid for directed networks. This is undirected."
                );

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr, false, "Cyclical triads"
    );

    return;
    
}
    
// Density --------------------------------------------------------------
inline void counter_density(NetCounters * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return
            1.0/(Array.nrow() * (Array.ncol() - 1.0)) / (
                (Array.D()->directed)? 1.0 : 2.0
            );
        
    };
    
    // Preparing the counter data and returning. We make sure that the memory is 
    // released so we set delete_data = true.
    counters->add_counter(
        tmp_count, nullptr, nullptr, false, "Density",
        "Proportion of present ties"
    );

    return ;
    
}

// idegree1.5  -------------------------------------------------------------
inline void counter_idegree15(NetCounters * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        // In case of the first, we need to add
        if (Array.col(j).size() == 1u)
            return 1.0;
        
        return 
            pow(static_cast<double> (Array.col(j).size()), 1.5) -
            pow(static_cast<double> (Array.col(j).size() - 1), 1.5)
            ;
        
    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr, false, "Indegree^(1.5)"
    );

    return;
    
}

// odegree1.5  -------------------------------------------------------------
inline void counter_odegree15(NetCounters * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        // In case of the first, we need to add
        if (Array.row(i).size() == 1u)
            return 1.0;
        
        return 
            pow(static_cast<double>(Array.row(i).size()), 1.5) -
            pow(static_cast<double>(Array.row(i).size() - 1), 1.5)
            ;
        
    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr, false, "Outdegree^(1.5)"
    );

    return;
    
}


// -----------------------------------------------------------------------------
/**@brief Sum of absolute attribute difference between ego and alter */
inline void counter_absdiff(
        NetCounters * counters,
        uint attr_id,
        double alpha = 1.0
    ) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return std::pow(std::fabs(
                Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][i] - 
                    Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][j]
        ), NET_C_DATA_NUM(0u));
        
    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (Array.D()->vertex_attr.size() == 0u)
            throw std::range_error("No attributes in the Array.");
        
        if ((NET_C_DATA_IDX(0u) != 0u) && (Array.D()->vertex_attr.size() <= (NET_C_DATA_IDX(0u) - 1u)))
            throw std::range_error("Attribute index out of range.");
        
        return 0.0;
        
    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new NetCounterData({attr_id}, {alpha}),
        true, "Absdiff"
    );
    
    return;
    
}
    
// -----------------------------------------------------------------------------
/**@brief Sum of attribute difference between ego and alter to pow(alpha)*/
inline void counter_diff(
        NetCounters * counters,
        uint attr_id,
        double alpha     = 1.0,
        double tail_head = true
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return std::pow(NET_C_DATA_NUM(1u) * (
                Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][i] - 
                    Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][j]
        ), NET_C_DATA_NUM(0u));
        
    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (Array.D()->vertex_attr.size() == 0u)
            throw std::range_error("No attributes in the Array.");
        
        if ((NET_C_DATA_IDX(0u) != 0u) && (Array.D()->vertex_attr.size() <= (NET_C_DATA_IDX(0u) - 1u)))
            throw std::range_error("Attribute index out of range.");
        
        return 0.0;
        
    };
    
    counters->add_counter(
        tmp_count, tmp_init,
        new NetCounterData({attr_id}, {alpha, tail_head ? 1.0: -1.0}),
        true, "Absdiff^(" + std::to_string(alpha) + ")"
    );
    
    return;
    
}

// Nodeicov, nodeocov, and Nodematch -------------------------------------------
NETWORK_COUNTER(init_single_attr)
{
    
    if (Array.D() == nullptr)
        throw std::logic_error("The array data has not been initialized");
    
    if (Array.D()->vertex_attr.size() == 0u)
        throw std::range_error("No attributes in the Array.");
    
    if ((NET_C_DATA_IDX(0u) != 0u) && (Array.D()->vertex_attr.size() <= (NET_C_DATA_IDX(0u) - 1u)))
        throw std::range_error("Attribute index out of range.");
    
    return 0.0;
    
}

// -----------------------------------------------------------------------------
//*@brief Attribute sum over receiver nodes */
inline void counter_nodeicov(NetCounters * counters, uint attr_id)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][j];
        
    };
    
    counters->add_counter(
        tmp_count, init_single_attr,
        new NetCounterData({attr_id}, {}),
        true, "nodeicov", "Sum of ego attribute"
    );
      
    return;

}

// -----------------------------------------------------------------------------
//*@brief Attribute sum over sender nodes */
inline void counter_nodeocov(NetCounters * counters, uint attr_id)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][i];
        
    };
    
    counters->add_counter(
        tmp_count, init_single_attr,
        new NetCounterData({attr_id}, {}),
        true, "nodeocov", "Sum of alter attribute"
    );
    
    return;

}

// -----------------------------------------------------------------------------
//*@brief Attribute sum over receiver and sender nodes */
inline void counter_nodecov(NetCounters * counters, uint attr_id)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][i] +
            Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][j];
        
    };
    
    counters->add_counter(
        tmp_count, init_single_attr,
        new NetCounterData({attr_id}, {}),
        true, "nodecov", "Sum of nodes covariates"
    );
    
    return;
}

// -----------------------------------------------------------------------------
//* @brief Number of homophililic ties */
inline void counter_nodematch(NetCounters * counters, uint attr_id)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return 
        (
                Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][i] == 
                    Array.D()->vertex_attr[NET_C_DATA_IDX(0u)][j]
        )? 1.0 : 0.0;
        
    };
    
    // Preparing the counter data and returning. We make sure that the memory is 
    // released so we set delete_data = true.
    counters->add_counter(
        tmp_count, init_single_attr,
        new NetCounterData({attr_id}, {}),
        true, "Homophily", "Number of homophilic ties"
    );
    
    return ;
    
}

// -----------------------------------------------------------------------------
/** @brief Counts number of vertices with a given in-degree */
inline void counter_idegree(
        NetCounters * counters,
        std::vector< uint > d)
{

    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        uint d = Array.col(j).size();
        if (d == NET_C_DATA_IDX(0u))
            return 1.0;
        else if (d == (NET_C_DATA_IDX(0u) + 1))
            return -1.0;
        
        return 0.0;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!Array.D()->directed)
            throw std::logic_error("-odegree- counter is only valid for directed graphs");
        
        if (NET_C_DATA_IDX(0u) == 0u)
            return static_cast<double>(Array.nrow());
        
        return 0.0;

    };
    
    for (auto iter = d.begin(); iter != d.end(); ++iter)
        counters->add_counter(
            tmp_count, tmp_init,
            new NetCounterData({*iter}, {}),
            true, "Nodes indeg " + std::to_string(*iter),
            "Number of nodes with indigree " + std::to_string(*iter)
        );
    
    return;  

}

// -----------------------------------------------------------------------------
/**@brief Counts number of vertices with a given out-degree */
inline void counter_odegree(
        NetCounters * counters,
        std::vector<uint> d
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        uint d = Array.row(i).size();
        if (d == NET_C_DATA_IDX(0u))
            return 1.0;
        else if (d == (NET_C_DATA_IDX(0u) + 1))
            return -1.0;
        
        return 0.0;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!Array.D()->directed)
            throw std::logic_error("-odegree- counter is only valid for directed graphs");
        
        if (NET_C_DATA_IDX(0u) == 0u)
            return static_cast<double>(Array.nrow());
        
        return 0.0;

    };
        
        
    for (auto iter = d.begin(); iter != d.end(); ++iter) 
        counters->add_counter(
            tmp_count, tmp_init,
            new NetCounterData({*iter}, {}),
            true, "Nodes w/ outdeg " + std::to_string(*iter),
            "Number of nodes with outdegree " + std::to_string(*iter)
        );
    
    return;  
    
}
    
// -----------------------------------------------------------------------------
/** @brief Counts number of vertices with a given out-degree */
inline void counter_degree(
        NetCounters * counters,
        std::vector<uint> d
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {
        
        uint d = Array.row(i).size();
        if (d == NET_C_DATA_IDX(0u))
            return 1.0;
        else if (d == (NET_C_DATA_IDX(0u) + 1))
            return -1.0;
        
        return 0.0;
    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init) {
        
        if (Array.D() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (Array.D()->directed)
            throw std::logic_error("-degree- counter is only valid for undirected graphs");
        
        if (NET_C_DATA_IDX(0u) == 0u)
            return static_cast<double>(Array.nrow());
        
        return 0.0;
    };
    
    
    for (auto iter = d.begin(); iter != d.end(); ++iter)
    {
        counters->add_counter(
            tmp_count, tmp_init,
            new NetCounterData({*iter}, {}),
            true
        );
    }
    
    return;  
}

#include "network-css.hpp"

///@}


/**
 * @name Rules for network models
 * @param rules A pointer to a `NetRules` object (`Rules`<`Network`, `bool`>).
 */
///@{
// -----------------------------------------------------------------------------
/**@brief Number of edges */
inline void rules_zerodiag(NetRules * rules) {
    
    NETWORK_RULE_LAMBDA(no_self_tie) {
        return i != j;
    };
    
    rules->add_rule(no_self_tie);
    
    return;
}

///@}

///@}

#undef NET_C_DATA_IDX
#undef NET_C_DATA_NUM

#endif
