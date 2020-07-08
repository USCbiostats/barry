#include "../counters-bones.hpp"
#include "../support.hpp"

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

typedef BArray<bool, NetworkData> Network;
typedef Counter<Network, std::vector<uint> > NetCounter;
typedef Support<Network, std::vector<uint> > NetSupport;

#define NETWORK_COUNTER(a) inline double (a) \
(const Network * Array, uint i, uint j, std::vector<uint> * data)
  

// Edges counter
NETWORK_COUNTER(count_edges) {
  return 1.0;
} 

NetCounter edges(count_edges);
 
// Isolates counter
NETWORK_COUNTER(count_isolates) {
  
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
  
}

NETWORK_COUNTER(init_isolates) {
  return (double) (Array->N);
}

NetCounter isolates(count_isolates, init_isolates);

// Mutuals -------------------------------------------------------------------
NETWORK_COUNTER(init_mutual) {
  
  if (Array->N != Array->M)
    throw std::logic_error("The -mutual- counter only works on square arrays.");

  if (!Array->data->directed)
    throw std::logic_error("The -mutual- counter only works on directed (non-symmetric) arrays.");
  
  return 0.0;
}

NETWORK_COUNTER(count_mutual) {

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
  
}

NetCounter mutual(count_mutual, init_mutual);

// 2-istars
NETWORK_COUNTER(count_istar2) {
 
  // Need to check the receiving, if he/she is getting a new set of stars
  // when looking at triads

  if (A_COL(j).size() == 1u)
    return 0.0;
  
  return ((double) A_COL(j).size() - 1.0);
 
  // return 0.0; 
}

NetCounter istar2(count_istar2);

// 2-ostars
NETWORK_COUNTER(count_ostar2) {
  
  // Need to check the receiving, if he/she is getting a new set of stars
  // when looking at triads
  
  if (A_ROW(i).size() == 1u)
    return 0.0;
  
  return ((double) A_ROW(i).size() - 1.0);
  
  // return 0.0; 
}

NetCounter ostar2(count_ostar2);

// ttriads
NETWORK_COUNTER(count_ttriads) {
  
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

}

NetCounter ttriads(count_ttriads);

// Cycle triads --------------------------------------------------------------
NETWORK_COUNTER(count_ctriads) {
  
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
  
}

NetCounter ctriads(count_ctriads);

// Density --------------------------------------------------------------
NETWORK_COUNTER (count_density) {
  return 1.0/(Array->N * (Array->M - 1));
}

NetCounter density(count_density);

// idegree1.5  -------------------------------------------------------------
NETWORK_COUNTER(count_idegree15) {
  
  // In case of the first, we need to add
  if (A_COL(j).size() == 1u)
    return 1.0;
  
  return 
    pow((double) A_COL(j).size(), 1.5) - pow((double) A_COL(j).size() - 1, 1.5)
  ;
}

NetCounter idegree15(count_idegree15);

// odegree1.5  -------------------------------------------------------------
NETWORK_COUNTER(count_odegree15) {
  
  // In case of the first, we need to add
  if (A_ROW(i).size() == 1u)
    return 1;
  
  return 
    pow((double) A_ROW(i).size(), 1.5) - pow((double) A_ROW(i).size() - 1, 1.5)
  ;
}

NetCounter odegree15(count_odegree15);

// Nodematch -------------------------------------------------------------------
NETWORK_COUNTER(count_nodematch) {
  
  return 
    (
        Array->data->vertex_attr[(*data)[0u]][i] == 
          Array->data->vertex_attr[(*data)[0u]][j]
    )? 1.0 : 0.0;
  
}

NETWORK_COUNTER(init_nodematch) {
  
  if (data == nullptr)
    throw std::logic_error("data for the counter must be specified.");
  
  if (data->size() != 1u)
    throw std::logic_error("data of the counter must be of size 1.");
  
  if (Array->data == nullptr)
    throw std::logic_error("data for the array must be specified.");
  
  if (Array->data->vertex_attr.size() == 0u)
    throw std::range_error("No attributes in the Array.");
  
  if (((*data)[0u] != 0u) && (Array->data->vertex_attr.size() <= ((*data)[0u] - 1u)))
    throw std::range_error("Attribute index out of range.");
  
  return 0.0;
  
}

NetCounter nodematch(count_nodematch, init_nodematch);

#endif
