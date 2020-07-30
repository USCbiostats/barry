#include <iostream>
#include <ostream>
#include "../include/barry.hpp"

template<typename T>
inline void print(const std::vector< T > & x) {
  std::cout << "[";
  for (auto iter = x.begin(); iter != x.end(); ++iter)
    std::cout << *iter << ", ";
  std::cout << "]" << std::endl;
  return; 
}

inline std::vector< double > keygen_phylo(
    const phylocounters::PhyloArray & Array_
  ) {
  return {
    (double) Array_.N, (double) Array_.M,
    (double) Array_.data->states[0],
    (double) Array_.data->states[1]
    };
}

typedef std::vector< unsigned int > vuint;

int main() {
  
  /**Representing the tree
   * [0]__[1]__[3]__[4]
   *   |    |    |
   *   |    |    |__[5]
   *   |    |    
   *   |    |_______[6]
   *   |
   *   |___[2]______[7]
   *         |
   *         |______[8]
   *
   * All four ancestor nodes have the same type (duplication events).
   */
  std::vector< std::pair<uint,uint> > edgelist = {
    {0,1}, {0,2},
    {1,3}, {1,6},
    {2,7}, {2,8},
    {3,4}, {3,5}
  };
  
  /**All nodes will share these four states (0,0), (1,0), (0,1), (1,1)*/
  phylocounters::PhyloArray n0(2,2);
  phylocounters::PhyloArray n1(2,2);
  phylocounters::PhyloArray n2(2,2);
  phylocounters::PhyloArray n3(2,2);
  
  n0.set_data(new phylocounters::NodeData({1u,1u}, {false,false} ), true);
  n1.set_data(new phylocounters::NodeData({1u,1u}, {true,false} ), true);
  n2.set_data(new phylocounters::NodeData({1u,1u}, {false,true} ), true);
  n3.set_data(new phylocounters::NodeData({1u,1u}, {true,true} ), true);
  
  
  // We only generate a single powerset since transitions are shared?
  phylocounters::PhyloModel model;
  model.set_keygen(keygen_phylo);
  
  // Activating the storage of powersets (because we'll need it!)
  model.store_psets();
  
  // Adding terms (gains/losses for each)
  phylocounters::counter_gains(&model.counters, {0, 1}); 
  phylocounters::counter_loss(&model.counters, {0, 1});
  
  // Now it is interesting: neofun and subfun
  phylocounters::counter_neofun(&model.counters, 0, 1);
  phylocounters::counter_subfun(&model.counters, 0, 1);
  
  // Adding the data! while forcing it to keep different counters
  std::vector< unsigned int > idx(4u);
  idx[0u] = model.add_array(n0, true);
  idx[1u] = model.add_array(n1, true);
  idx[2u] = model.add_array(n2, true);
  idx[3u] = model.add_array(n3, true);
  
  // Printing the first one
  std::cout << "pset_stat.size()      = " << model.pset_stats.size() << std::endl;
  std::cout << "pset_stat[0].size()   = " << model.pset_stats[0].size() << std::endl;
  std::cout << "pset_arrays[0].size() = " << model.pset_arrays[0].size() << std::endl;

  
  print(model.pset_stats[0][0]);
  
  std::cout << "The likelihood for model with parameters 1 equals:\n " << std::endl;
  print(idx);
  print(
    (std::vector<double>) {
      model.likelihood({1,1,1,1,1,1}, idx[0u], false),
      model.likelihood({1,1,1,1,1,1}, idx[1u], false),
      model.likelihood({1,1,1,1,1,1}, idx[2u], false),
      model.likelihood({1,1,1,1,1,1}, idx[3u], false)
    }
  );
   
  
  print(model.normalizing_constants);
  
  
  return 0;
}
 
