#include <Rcpp.h>
#include "../include/barray.hpp"

using namespace Rcpp;

template <typename T>
using Vec = std::vector< T >;

// [[Rcpp::export]]
List counter_phylo(
    const LogicalVector & x,
    int nfun,
    int noffspring,
    const NumericVector & blenghts
    ) {
  
  // Initializing the node
  phylocounters::PhyloArray tree(nfun, noffspring);
  phylocounters::NodeData data(as<Vec<double>>(blenghts), as<Vec<bool>>(x)); 
  tree.data = &data;
  
  // Setting counters, one per function
  barray::Counter<phylocounters::PhyloArray, Vec<uint>> counter0 = phylocounters::gains;
  barray::Counter<phylocounters::PhyloArray, Vec<uint>> counter1 = phylocounters::gains;
  counter0.data = new Vec<uint>({0u});
  counter1.data = new Vec<uint>({1u});
  
  barray::Counter<phylocounters::PhyloArray, Vec<uint>> counter2 = phylocounters::loss;
  barray::Counter<phylocounters::PhyloArray, Vec<uint>> counter3 = phylocounters::loss;
  counter2.data = new Vec<uint>({0u});
  counter3.data = new Vec<uint>({1u});
  
  barray::Counter<phylocounters::PhyloArray, Vec<uint>> counter4 = phylocounters::subfun;
  counter4.data = new Vec<uint>({0u, 1u});
  
  barray::Counter<phylocounters::PhyloArray, Vec<uint>> counter5 = phylocounters::cogain;
  counter5.data = new Vec<uint>({0u, 1u});
  
  barray::Counter<phylocounters::PhyloArray, Vec<uint>> counter6 = phylocounters::longest;
  counter6.data = new Vec<uint>(0u);
  
  
  barray::Support<phylocounters::PhyloArray, Vec<uint>> support(&tree);
  support.add_counter(counter0);
  support.add_counter(counter1);
  support.add_counter(counter2);
  support.add_counter(counter3);
  support.add_counter(counter4);
  support.add_counter(counter5);
  support.add_counter(counter6);
  
  // Computing and retrieving
  std::vector< std::vector< double > > observed(0u);
  support.calc(0u, true, nullptr, &observed);
  
  delete counter0.data;
  delete counter1.data;
  delete counter2.data;
  delete counter3.data;
  delete counter4.data;
  delete counter5.data;
  delete counter6.data;
  counter0.data = nullptr;
  counter1.data = nullptr;
  counter2.data = nullptr;
  counter3.data = nullptr;
  counter4.data = nullptr;
  counter5.data = nullptr;
  counter6.data = nullptr;
  
  // Generating the entries
  barray::Counts_type ans = support.support.get_entries();
  
  List res(ans.size());
  for (unsigned int i = 0u; i < res.size(); ++i) {
    res[i]       = List::create(
      _["x"]     = ans.at(i).first,
      _["count"] = ans.at(i).second
    );
  }
  
  return List::create(
    _["suff_stats"] = res,
    _["pset_stats"] = wrap(observed)
  );
  
}


/***R

snames <- c(
  "gains0",
  "gains1",
  "loss0",
  "loss1",
  "subfun",
  "cogain01",
  "longest",
  "counts"
)

wrap <- function(x) {
  structure(
    do.call(rbind, lapply(x$suff_stats, unlist)),
    dimnames = list(NULL, snames)
  )
}

# Other parameters
nfun        <- 2
noffspring  <- 4
blengths    <- rep(1, noffspring)
blengths[1] <- 2

# Case 1 (0,0)
wrap(counter_phylo(c(FALSE, FALSE), nfun = nfun, noffspring = noffspring, blengths))


# Case 2 (1, 0)
wrap(counter_phylo(c(TRUE, FALSE), nfun = nfun, noffspring = noffspring, blengths))

# Case 3 (0, 1)
wrap(counter_phylo(c(FALSE, TRUE), nfun = nfun, noffspring = noffspring, blengths))

# Case 3 (1, 1)
(ans0 <- wrap(counter_phylo(c(TRUE, TRUE), nfun = nfun, noffspring = noffspring, blengths)))
sum(ans0[,"counts"])

blengths2 <- blengths
blengths2[2] <- 2
(ans1 <- wrap(counter_phylo(c(TRUE, TRUE), nfun = nfun, noffspring = noffspring, blengths2)))
sum(ans1[,"counts"])

*/

// 
// #endif
