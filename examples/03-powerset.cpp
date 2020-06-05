#include <Rcpp.h>
#include "../include/barray.hpp"
using namespace Rcpp;
 

// [[Rcpp::export]]
SEXP psets(int n, int m) {
    
  Rcpp::XPtr< barray::PowerSet<>> xptr(
    new barray::PowerSet<>((uint) n, (uint) m),
    true
  );
  
  // Generating the powerset 
  xptr->calc();
  // std::cout << "Size: " << sizeof(*xptr) << std::endl;
  return xptr;
}


// [[Rcpp::export]]
List get_data(SEXP x) {
  Rcpp::XPtr< barray::PowerSet<>> xptr(x);

  // Generating the data
  List ans(xptr->data.size());
  uint counter = 0u;
  for (auto iter = xptr->data.begin(); iter != xptr->data.end(); ++iter) {
    
    barray::Entries<bool> set = iter->get_entries();
    
    ans[counter++] = List::create(
      _["source"] = set.source,
      _["target"] = set.target,
      _["val"] = set.val
    );
    
  }
    
  
  return ans;
  
}

/***R
N <- 4
M <- 3

PS_2_3 <- psets(N, M)
x <- ergmito::powerset(4)

microbenchmark::microbenchmark(
  PS_2_3 <- psets(N,M),
  x <- ergmito::powerset(4)
)


stop()

PS <- psets(3,1)
PS <- get_data(PS)

PS <- lapply(PS, function(p.) {
  ans <- matrix(0, nrow = 3, ncol = 1)
  ans[cbind(p.$source, p.$target) + 1] <- p.$val
  ans
})

table(sapply(PS, paste, collapse = ""))
table(table(sapply(PS, paste, collapse = ""))) == 2 ^ (3)

*/


