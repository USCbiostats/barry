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

/* This is for hashing. Ultimately, we want to know what features do we need
 * to look at when we are finding duplicates. In this case there are only
 * four:
 * 1. Number of row,
 * 2. Number of columns,
 * 3. Parent state of function 0.
 * 4. Parent state of function 1.
 */
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
   * [0]__[1]__[3]__[4]__[5]
   *   |    |    |    |
   *   |    |    |    |__[6]
   *   |    |    |
   *   |    |    |_______[7]
   *   |    |
   *   |    |____________[8]
   *   |
   *   |__[2]____________[9]
   *        |
   *        |____________[10]
   *
   * All six ancestor nodes have the same type (duplication events).
   */
  // std::vector< std::pair<uint,uint> > edgelist = {
  //   {0,1}, {0,2},
  //   {1,3}, {1,6},
  //   {2,7}, {2,8},
  //   {3,4}, {3,5}
  // };

  /**All nodes will share these five states (0,0), (1,0), (0,1), (1,1), (1,1)
   * so the arrays are of two by to (two functions x two siblings)
   */
  std::cout << "Creating arrays" << std::endl;
  phylocounters::PhyloArray n0(2,2);
  phylocounters::PhyloArray n1(2,2);
  phylocounters::PhyloArray n2(2,2);
  phylocounters::PhyloArray n3(2,2);
  phylocounters::PhyloArray n4(2,2);


  /* We now start the counter. To differentiate objects, we use the
   * keygen_phylo function defined earlier. This function receives an array
   * and returns a vector.
   */
  std::cout << "Preparing the model" << std::endl;
  phylocounters::PhyloModel model;
  model.set_keygen(keygen_phylo);

  // Activating the storage of powersets (because we'll need it!)
  model.store_psets();

  std::cout << "Adding counters" << std::endl;
  // Adding terms (gains/losses for each)
  phylocounters::counter_gains(&model.counters, {0, 1});
  phylocounters::counter_loss(&model.counters, {0, 1});

  // Now it is interesting: neofun and subfun
  phylocounters::counter_neofun(&model.counters, 0, 1);
  phylocounters::counter_subfun(&model.counters, 0, 1);

  /* We set the last argument as true so that the destructor takes care
   * of the cleaning once the arrays are deleted.
   * Branch lengths are 1
   */
  std::cout << "Adding nodes" << std::endl;
  n0.set_data(new phylocounters::NodeData({1u}, {false,false} ), true);
  n1.set_data(new phylocounters::NodeData({1u}, {true,false} ), true);
  n2.set_data(new phylocounters::NodeData({1u}, {false,true} ), true);
  n3.set_data(new phylocounters::NodeData({1u}, {true,true} ), true);
  n4.set_data(new phylocounters::NodeData({1u}, {true,true} ), true);

  // Adding the data!
  std::cout << "Adding the data" << std::endl;

  std::vector< unsigned int > idx(5u);
  idx[0u] = model.add_array(n0, false);
  idx[1u] = model.add_array(n1, false);
  idx[2u] = model.add_array(n2, false);
  idx[3u] = model.add_array(n3, false);
  idx[4u] = model.add_array(n4, false);

  // Printing the first one
  std::cout << "The number of unique statistics should equal to the number of unique supports:" << std::endl;
  std::cout << "pset_stat.size()      = " << model.pset_stats.size() << std::endl;
  std::cout << "pset_stat[0].size()   = " << model.pset_stats[0].size() << std::endl;
  std::cout << "pset_arrays[0].size() = " << model.pset_arrays[0].size() << std::endl;

  print(model.pset_stats[0][0]);

  std::vector< double > model_parameters = {.9, .8, .02, .05, .5, .7};

  std::cout << "Indices: " << std::endl;
  print(idx);
  std::cout << "The likelihood for model with parameters 1 equals: " << std::endl;
  print(
    (std::vector<double>) {
      model.likelihood(model_parameters, idx[0u], false),
      model.likelihood(model_parameters, idx[1u], false),
      model.likelihood(model_parameters, idx[2u], false),
      model.likelihood(model_parameters, idx[3u], false),
      model.likelihood(model_parameters, idx[4u], false)
    }
  );

  std::cout << "Normalizing constants: " << std::endl;
  print(model.normalizing_constants);

  for (unsigned int i = 0u; i < model.n_arrays(); ++i) {
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Looking at the support of array " << i << std::endl;
    std::cout << "gain0, gain1, loss0, loss1, neofun, subfun" << std::endl;
    model.print_stats(i);
  }

  return 0;
}
