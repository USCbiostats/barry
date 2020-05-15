#include <Rcpp.h>
#include "../include/barray.hpp"
using namespace Rcpp;


// [[Rcpp::export]]
List psets(int n, int m) { 
    
  Rcpp::XPtr< barray::LBArray< bool >> xptr(
    new barray::LBArray< bool >((uint) n, (uint) m),
    true
  );
  
  // Generating the powerset 
  xptr->pset();
  
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
M <- 5
PS_2_3 <- psets(N,M)

PS_2_3 <- lapply(PS_2_3, function(p.) {
  ans <- matrix(0, nrow = N, ncol = M)
  ans[cbind(p.$source, p.$target) + 1] <- p.$val
  ans
})

table(sapply(PS_2_3, paste, collapse = ""))
table(table(sapply(PS_2_3, paste, collapse = ""))) == 2 ^ (N*M)

*/


