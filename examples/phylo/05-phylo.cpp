#include <Rcpp.h>
#include "../include/barry.hpp"

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
  phylocounters::PhyloSupport support(&tree);
  
  // Adding counters
  phylocounters::counter_gains(support.counters, 0u);
  phylocounters::counter_gains(support.counters, 1u);
  phylocounters::counter_loss(support.counters, 0u);
  phylocounters::counter_loss(support.counters, 1u);
  phylocounters::counter_subfun(support.counters, 0u, 1u);
  phylocounters::counter_cogain(support.counters, 0u, 1u);
  phylocounters::counter_longest(support.counters);
  
  // Computing and retrieving
  std::vector< std::vector< double > > observed(0u);
  support.calc(0u, true, nullptr, &observed);
  
  // Generating the entries
  barry::Counts_type ans = support.get_counts();
  
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
