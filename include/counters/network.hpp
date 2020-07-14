#include "../counters-bones.hpp"
#include "../support.hpp"
#include "../statscounter.hpp"

#ifndef BARRAY_NETWORK_H
#define BARRAY_NETWORK_H 1

/**@brief Data class for Networks.
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
  
  /**@brief Constructor using a single attribute
   * @param vertex_attr_ Double vector of length equal to the number of vertices
   * in the data.
   * @param directed_ When `true` the graph as treated as directed.
   */
  NetworkData(
    std::vector< double >  vertex_attr_,
    bool directed_ = true
  ) : directed(directed_), vertex_attr(1u, vertex_attr_) {};
  
  /**@brief Constructor using multiple attributes
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


typedef BArray<bool, NetworkData> Network;
typedef Counter<Network, NetCounterData > NetCounter;
typedef Support<Network, NetCounterData > NetSupport;
typedef StatsCounter<Network, NetCounterData> NetStatsCounter;

#define NETWORK_COUNTER(a) inline double (a) \
(const Network * Array, uint i, uint j, NetCounterData * data)

#define NETWORK_COUNTER_LAMBDA(a) Counter_fun_type<Network, NetCounterData> a = \
  [](const Network * Array, uint i, uint j, NetCounterData * data)

// Edges counter ---------------------------------------------------------------
inline NetCounter counter_edges() {
  
  NETWORK_COUNTER_LAMBDA(count_edges) {
    return 1.0;
  };
  
  NetCounter edges(count_edges);
  return edges;
}


// Isolates counter ------------------------------------------------------------
inline NetCounter counter_isolates() {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    if (i == j)
      return 0.0;
    
    double res = 0.0;
    
    // i is sending its first tie
    if (A_ROW(i).size() == 1u && A_COL(i).size() == 0u)
      res -= 1.0;
    
    // j is receiving its first tie, meaning that he
    // has no other tie but i's?
    if (A_ROW(j).size() == 0u && A_COL(j).size() == 1u)
      res -= 1.0;
    
    return res;
  };
  
  NETWORK_COUNTER_LAMBDA(tmp_init) {
    return (double) (Array->N);
  };
  
  NetCounter tmp_counter(tmp_count, tmp_init);
  return tmp_counter;
}

// Mutuals -------------------------------------------------------------------
inline NetCounter counter_mutual() {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    // Is there any tie at ji? If not, then we have a new mutual!
    // but this only makes sence if the jth row and ith column exists
    // if ((Array->N > j) && (Array->M > i)) 
    if (i == j)
      return 0.0;
    
    // printf("Checking if it is empty or not at (%i, %i)... ", i, j);
    if (!Array->is_empty(j, i, false)) {
      // printf("Yes, mutual.\n");
      return 1.0;
    }
    // printf("No, no mutual.\n");
    
    return 0.0;
  };
  
  NETWORK_COUNTER_LAMBDA(tmp_init) {
    if (Array->N != Array->M)
      throw std::logic_error("The -mutual- counter only works on square arrays.");
    
    if (!Array->data->directed)
      throw std::logic_error("The -mutual- counter only works on directed (non-symmetric) arrays.");
    
    return 0.0;
  };
  
  NetCounter tmp_counter(tmp_count, tmp_init);
  return tmp_counter;
}


// 2-istars --------------------------------------------------------------------
inline NetCounter counter_istar2() {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    // Need to check the receiving, if he/she is getting a new set of stars
    // when looking at triads
    
    if (A_COL(j).size() == 1u)
      return 0.0;
    
    return ((double) A_COL(j).size() - 1.0);
  };
  
  NetCounter tmp_counter(tmp_count);
  return tmp_counter;
}

// 2-ostars --------------------------------------------------------------------
inline NetCounter counter_ostar2() {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    // Need to check the receiving, if he/she is getting a new set of stars
    // when looking at triads
    
    if (A_ROW(i).size() == 1u)
      return 0.0;
    
    return ((double) A_ROW(i).size() - 1.0);
  };
  
  NetCounter tmp_counter(tmp_count);
  return tmp_counter;
}


// ttriads ---------------------------------------------------------------------
inline NetCounter counter_ttriads() {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    // Self ties do not count
    if (i == j)
      return 0.0;
    
    double ans = 0.0;
    
    // Case 1: i-j, i-k, j-k
    if (A_ROW(j).size() < A_ROW(i).size()) {
      
      for (auto j_row = A_ROW(j).begin(); j_row != A_ROW(j).end(); ++j_row) 
        if ((j != j_row->first) && (i != j_row->first) && !Array->is_empty(i, j_row->first, false))
          ans += 1.0;
        
    } else {
      
      for (auto i_row = A_ROW(i).begin(); i_row != A_ROW(i).end(); ++i_row) 
        if ((i != i_row->first) && (i_row->first != j) && !Array->is_empty(j, i_row->first, false))
          ans += 1.0;
        
    }
    
    // Case 2: i-j, i-k, k-j  
    if (A_ROW(i).size() > A_COL(j).size()) {
      
      for (auto j_col = A_COL(j).begin(); j_col != A_COL(j).end(); ++j_col)
        if ((j != j_col->first) && (i != j_col->first) && !Array->is_empty(i, j_col->first, false))
          ans += 1.0;
        
    } else {
      
      for (auto i_row = A_ROW(i).begin(); i_row != A_ROW(i).end(); ++i_row) 
        if ((i != i_row->first) && (j != i_row->first) && !Array->is_empty(i_row->first, j, false))
          ans += 1.0;
        
    }
    
    // Case 3: 
    if (A_COL(i).size() > A_COL(j).size()) {
      
      for (auto j_col = A_COL(j).begin(); j_col != A_COL(j).end(); ++j_col)
        if ((j != j_col->first) && (i != j_col->first) && !Array->is_empty(j_col->first, i, false))
          ans += 1.0;
        
    } else {
      
      for (auto i_col = A_COL(i).begin(); i_col != A_COL(i).end(); ++i_col) 
        if ((i != i_col->first) && (j != i_col->first) && !Array->is_empty(i_col->first, j, false))
          ans += 1.0;
        
    }
    
    // The regular counter double counts
    return ans;
  };
  
  NETWORK_COUNTER_LAMBDA(tmp_init) {
    if (!(Array->data->directed))
      throw std::invalid_argument("The ttriads counter is only valid for directed networks. This is undirected.");
    return 0.0;
  };
  
  
  
  NetCounter tmp_counter(tmp_count, tmp_init);
  return tmp_counter;
}


// Cycle triads --------------------------------------------------------------
inline NetCounter counter_ctriads() {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    if (i == j)
      return 0.0;
    
    double ans = 0.0;
    if (A_COL(i).size() < A_ROW(j).size()) {
      
      for (auto i_col = A_COL(i).begin(); i_col != A_COL(i).end(); ++i_col) 
        if ((i != i_col->first) && (j != i_col->first) && !Array->is_empty(j, i_col->first, false))
          ans += 1.0;
        
    } else {
      
      for (auto j_row = A_ROW(j).begin(); j_row != A_ROW(j).end(); ++j_row) 
        if ((j != j_row->first) && (i != j_row->first) && !Array->is_empty(j_row->first, i, false))
          ans += 1.0;
        
    }
    
    return ans;
  };
  
  NETWORK_COUNTER_LAMBDA(tmp_init) {
    if (!(Array->data->directed))
      throw std::invalid_argument("The ctriads counter is only valid for directed networks. This is undirected.");
    return 0.0;
  };
  
  NetCounter tmp_counter(tmp_count, tmp_init);
  return tmp_counter;
  
}
  
// Density --------------------------------------------------------------
inline NetCounter counter_density() {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    
    return 1.0/(Array->N * (Array->M - 1));
    
  };
  
  // Preparing the counter data and returning. We make sure that the memory is 
  // released so we set delete_data = true.
  NetCounter tmp_counter(tmp_count);
  return tmp_counter;
  
}

// idegree1.5  -------------------------------------------------------------
inline NetCounter counter_idegree15() {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    
    // In case of the first, we need to add
    if (A_COL(j).size() == 1u)
      return 1.0;
    
    return 
      pow((double) A_COL(j).size(), 1.5) - pow((double) A_COL(j).size() - 1, 1.5)
      ;
    
  };
  
  NetCounter tmp_counter(tmp_count);
  return tmp_counter;
  
}

// odegree1.5  -------------------------------------------------------------
inline NetCounter counter_odegree15() {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    
    // In case of the first, we need to add
    if (A_ROW(i).size() == 1u)
      return 1.0;
    
    return 
      pow((double) A_ROW(i).size(), 1.5) - pow((double) A_ROW(i).size() - 1, 1.5)
      ;
    
  };
  
  NetCounter tmp_counter(tmp_count);
  return tmp_counter;
  
}


// Nodematch -------------------------------------------------------------------
inline NetCounter counter_absdiff(uint attr_id, double alpha = 1.0) {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    
    return std::pow(std::abs(
        Array->data->vertex_attr[NET_C_DATA_IDX(0u)][i] - 
          Array->data->vertex_attr[NET_C_DATA_IDX(0u)][j]
    ), NET_C_DATA_NUM(0u));
    
  };
  
  NETWORK_COUNTER_LAMBDA(tmp_init) {
    
    if (Array->data == nullptr)
      throw std::logic_error("data for the array must be specified.");
    
    if (Array->data->vertex_attr.size() == 0u)
      throw std::range_error("No attributes in the Array.");
    
    if ((NET_C_DATA_IDX(0u) != 0u) && (Array->data->vertex_attr.size() <= (NET_C_DATA_IDX(0u) - 1u)))
      throw std::range_error("Attribute index out of range.");
    
    return 0.0;
    
  };
  
  NetCounter tmp_counter(tmp_count);
  tmp_counter.data = new NetCounterData({attr_id}, {alpha});
  tmp_counter.delete_data = true;
  return tmp_counter;
  
}


// Nodeicov, nodeocov, and Nodematch -------------------------------------------
NETWORK_COUNTER(init_single_attr) {
  
  if (Array->data == nullptr)
    throw std::logic_error("data for the array must be specified.");
  
  if (Array->data->vertex_attr.size() == 0u)
    throw std::range_error("No attributes in the Array.");
  
  if ((NET_C_DATA_IDX(0u) != 0u) && (Array->data->vertex_attr.size() <= (NET_C_DATA_IDX(0u) - 1u)))
    throw std::range_error("Attribute index out of range.");
  
  return 0.0;
  
}

inline NetCounter counter_nodeicov(uint attr_id) {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    
    return Array->data->vertex_attr[NET_C_DATA_IDX(0u)][j];
    
  };
  
  NetCounter tmp_counter(tmp_count, init_single_attr);
  tmp_counter.data = new NetCounterData({attr_id}, {});
  tmp_counter.delete_data = true;
  return tmp_counter;
}

inline NetCounter counter_nodeocov(uint attr_id) {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    
    return Array->data->vertex_attr[NET_C_DATA_IDX(0u)][i];
    
  };
  
  NetCounter tmp_counter(tmp_count, init_single_attr);
  tmp_counter.data = new NetCounterData({attr_id}, {});
  tmp_counter.delete_data = true;
  return tmp_counter;
}

inline NetCounter counter_nodecov(uint attr_id) {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    
    return Array->data->vertex_attr[NET_C_DATA_IDX(0u)][i] +
      Array->data->vertex_attr[NET_C_DATA_IDX(0u)][j];
    
  };
  
  NetCounter tmp_counter(tmp_count, init_single_attr);
  tmp_counter.data = new NetCounterData({attr_id}, {});
  tmp_counter.delete_data = true;
  return tmp_counter;
}

inline NetCounter counter_nodematch(uint attr_id) {
  
  NETWORK_COUNTER_LAMBDA(tmp_count) {
    
    return 
    (
        Array->data->vertex_attr[NET_C_DATA_IDX(0u)][i] == 
          Array->data->vertex_attr[NET_C_DATA_IDX(0u)][j]
    )? 1.0 : 0.0;
    
  };
  
  // Preparing the counter data and returning. We make sure that the memory is 
  // released so we set delete_data = true.
  NetCounter nodematch(tmp_count, init_single_attr);
  nodematch.data = new NetCounterData({attr_id}, {});
  nodematch.delete_data = true;
  return nodematch;
  
}

#endif
