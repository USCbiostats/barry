#include <Rcpp.h>
#include "../include/barry.hpp"
using namespace Rcpp;

// Upper diagonal are blocked
template <typename Array_Type, typename Data_Type>
inline bool rule_lowertri(
    const Array_Type * a, barry::uint i, barry::uint j, Data_Type * d = nullptr
) {
  return i < j;
}

// Diagonal blocked
template <typename Array_Type, typename Data_Type>
inline bool rule_diag(
    const Array_Type * a, barry::uint i, barry::uint j, Data_Type * d = nullptr
) {
  return i == j;
}

// [[Rcpp::export]]
SEXP new_block(int i = 0) {

  typedef barry::Rules< barry::BArray<>, bool > rulet;
  Rcpp::XPtr< rulet > xptr( 
    new barry::Rules< barry::BArray<>, bool >(),
    true
  );
  
  // Adding two simple rules: Only lower triangle matrix with zero diagonal
  // are to be touched.
  xptr->add_rule(rule_lowertri<barry::BArray<>,bool>);
  xptr->add_rule(rule_diag<barry::BArray<>,bool>);
  
  return xptr;
}

// [[Rcpp::export]]
List get_sequence(SEXP x, int N, int K) {
  
  // Making space
  typedef unsigned int uint;
  std::vector< std::pair<uint,uint> > ans0(0);
  std::vector< std::pair<uint,uint> > ans1(0);
  
  // This function needs to be applied over an array
  Rcpp::XPtr< barry::Rules< barry::BArray<>, bool > > xptr(x);
  barry::BArray<> adjmat(N, K);
  xptr->get_seq(&adjmat, &ans0, &ans1);
  
  // Preparing the output to be exported to R
  IntegerMatrix res0(ans0.size(), 2u);
  IntegerMatrix res1(ans1.size(), 2u);
  int i = 0u;
  for (auto iter = ans0.begin(); iter != ans0.end(); ++iter) {
    res0(i, 0u) = iter->first;
    res0(i++, 1u) = iter->second;
  }
  i = 0u;
  for (auto iter = ans1.begin(); iter != ans1.end(); ++iter) {
    res1(i, 0u) = iter->first;
    res1(i++, 1u) = iter->second;
  }
  
  return List::create(_["free"] = res0, _["blocked"] = res1);
  
}



/***R

N <- 1e2
ptr <- new_block()
ans <- get_sequence(ptr, N, N)
m   <- matrix(0, N, N)
m[ans$free + 1L] <- 1

library(Matrix)
Matrix::image(as(m, "dgCMatrix"))


*/


