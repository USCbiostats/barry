#include <Rcpp.h>
#include "../include/barray.hpp"
#include "../include/lbarray-bones.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
SEXP suff_stats(const NumericMatrix & x) {
  
  Rcpp::XPtr< barray::SuffStats > xptr(
    new barray::SuffStats(),
    true
  );
  
  // Counting the stats
  for (barray::uint i = 0u; i < x.ncol(); ++i) {
    const std::vector< double > r(x.column(i).begin(), x.column(i).end());
    xptr->add(r);
  }
  
  return xptr;
}

// [[Rcpp::export]]
List get_suff_stats(SEXP x) {
  
  Rcpp::XPtr< barray::SuffStats > xptr(x);
  
  // Now, getting the data
  barray::vec_pair_dbl_uint ans = xptr->get_entries();
  
  List res(ans.size());
  for (unsigned int i = 0u; i < res.size(); ++i) {
    res[i] = List::create(
      _["x"] = ans.at(i).first,
      _["count"] = ans.at(i).second
    );
  }
  
  return res;
  
}

// A more elaborated example
// [[Rcpp::export]]
SEXP counter(
    int N, int M,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    const std::vector< double > & value
) {
  
  // Initializing the Binary array, and also the the suffstats counter
  const barray::BArray Array((uint) N, (uint) M, source, target, value);
  barray::SuffStats stats;

  
  // Creating the counter object; 
  barray::StatsCounter dat(&Array, &stats);
  
  // Adding functions
  dat.add_counter(barray::counters::edges);
  dat.add_counter(barray::counters::mutual);
  dat.add_counter(barray::counters::isolates);
  
  // Fingers crossed
  dat.count_all();
  
  barray::vec_pair_dbl_uint ans = stats.get_entries();
  
  List res(ans.size());
  for (unsigned int i = 0u; i < res.size(); ++i) {
    res[i] = List::create(
      _["x"] = ans.at(i).first,
      _["count"] = ans.at(i).second
    );
  }
  
  return res;
  
}



/***R

# Example with sufficient statistics
set.seed(123)
x <- seq(0, 10, length.out = 20)
x <- matrix(sample(x, 500*3, TRUE), ncol = 3)
xt <- t(x)

microbenchmark::microbenchmark(
  ans0 <- suff_stats(xt),
  ans1 <- table(x[,1], x[,2], x[,3])
)

ptr <- suff_stats(xt)
ans0 <- get_suff_stats(ptr)
ans1 <- t(sapply(ans0, function(a) c(a[[1]], a[[2]])))
maxcount <- ans1[order(-ans1[,4]),][1,]

nrow(ans1)

# Example with functions
set.seed(123)
N <- 100
M <- N

nedges  <- 1000
source <- sample.int(N, nedges, replace = TRUE)
target <- sample.int(M, nedges, replace = TRUE)
values <- runif(nedges)

# Removing diagonal
idx <- which(source != target)
source <- source[idx]
target <- target[idx]
values <- values[idx]

# Adding a few mutual
source <- c(source, target[1:20])
target <- c(target, source[1:20])
values <- c(values, 1:20)

el <- counter(N, M, source - 1L, target - 1L, values)

# Comparing with ergm
mat <- matrix(0, nrow = N, ncol = M)
mat[cbind(source, target)] <- 1L
microbenchmark::microbenchmark(
  # ergmito::count_stats(mat ~ edges + mutual),
  ergm::summary_formula(mat ~ edges + mutual + isolates),
  counter(N, M, source - 1L, target - 1L, values)
)


*/


