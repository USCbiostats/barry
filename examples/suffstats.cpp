#include <Rcpp.h>
#include "../include/barray.hpp"
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



/***R

# Example with sufficient statistics
set.seed(123)
x <- seq(0, 10, length.out = 20)
x <- matrix(sample(x, 50000*3, TRUE), ncol = 3)
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
*/


