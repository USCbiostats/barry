#include "barray.hpp"

#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H

// template <typename Cell_type>
class PhyloNode : public barray::BArray<bool> {
  
private:
  std::vector< double > lengths;
  barray::PowerSet<bool> states;
  std::vector< barray::Support< bool > > stats;
public:
  // PhyloNode() : BArray<bool>(0u, 0u) {} ;
  PhyloNode(barray::uint N_, barray::uint M_) : BArray<bool>(N_, M_), states(N_, M_) {};
  ~PhyloNode() {};
  
  void init(std::vector< barray::Counter_type< bool > > counters) {
    
    // Adding the needed counters to each state space
    for (barray::uint i = 0; i < stats.size(); ++i) {
      for (auto iter = counters.begin(); iter != counters.end(); ++iter) 
        stats[i].add_counter(*iter);
  
      // Calculagin the support
      stats[i].calc();
    }
    
    return;
  }
};

#endif
