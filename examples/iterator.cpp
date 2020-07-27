#include <Rcpp.h>
#include "../include/barry.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
SEXP new_block(const IntegerVector & row, const IntegerVector & col, int N, int M) {

  if (row.size() != col.size())
    stop("row and col should have the same size.");
    
  // Preparing the data
  typedef unsigned int uint;
  std::vector< std::pair< uint, uint > > dat;
  for (uint i = 0u; i < row.size(); ++i)
    dat.push_back({row.at(i), col.at(i)});
  
  Rcpp::XPtr< barry::CellSeq > xptr( 
    new barry::CellSeq(dat, N, M),
    true
  );
  
  return xptr;
}

// [[Rcpp::export]]
IntegerMatrix get_sequence(SEXP x) {
  
  Rcpp::XPtr< barry::CellSeq > xptr(x);
  typedef unsigned int uint;
  const std::vector< std::pair<uint,uint> > * ans = xptr->get_seq();
  
  IntegerMatrix res(ans->size(), 2u);
  int i = 0u;
  for (auto iter = ans->begin(); iter != ans->end(); ++iter) {
    res(i, 0u) = iter->first;
    res(i++, 1u) = iter->second;
  }
  
  return res;
  
}



/***R

# No diagonal
N<-M<-10
rows <- 0:(N-1)
cols <- rows

# microbenchmark::microbenchmark(
# build = {
ptr <- new_block(rows, cols, N, M)
ans <- get_sequence(ptr) + 1
# }, unit = "s")
m   <- matrix(0, N, N)
m[ans] <- 1
m

# Block diagonal
dat <- ans[ans[,1] <= ans[,2],]
ptr <- new_block(dat[,1] - 1, dat[,2] - 1, N, M)
ans <- get_sequence(ptr) + 1
m   <- matrix(0, N, N)
m[ans] <- 1
m

*/


