#include <Rcpp.h>
#include "../include/barry.hpp"
using namespace Rcpp; 

// [[Rcpp::export]]
SEXP new_Array(  
    int N, int M,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    const std::vector< double > & value
) {
  
  Rcpp::XPtr< barry::BArray<double,double> > ptr(
    new barry::BArray<double,double>((uint) N, (uint) M, source, target, value),
    true
  );
  
  return ptr;
  
}

// [[Rcpp::export]]
double get_cell(SEXP x, int i, int j) {
  Rcpp::XPtr< barry::BArray<double,double> > xptr(x);
  return xptr->get_cell(i, j);
}

// Returning the i-th row
// [[Rcpp::export]]
NumericVector get_row(SEXP x, int i) {
  
  Rcpp::XPtr< barry::BArray<double,double> > xptr(x);
  NumericVector ans(xptr->M, 0);
  const barry::Row_type<double> * m = xptr->get_row(i);
  
  for (auto row = m->begin(); row != m->end(); ++row)
    ans[row->first] = row->second.value;
  
  return ans;
  
}

// [[Rcpp::export]]
NumericVector get_col(SEXP x, int i) {
  
  Rcpp::XPtr< barry::BArray<double,double> > xptr(x);
  NumericVector ans(xptr->N, 0);
  const barry::Col_type< double > * m = xptr->get_col(i);
  
  for (auto row = m->begin(); row != m->end(); ++row) {
    ans[row->first] = row->second->value; 
  }
  
  return ans;
  
}

// [[Rcpp::export]]
List get_entries(const SEXP & x) {
  
  Rcpp::XPtr< barry::BArray<double,double> > xptr(x);
  barry::Entries<double> res = xptr->get_entries();
  
  return List::create(
    _["source"] = res.source,
    _["target"] = res.target, 
    _["val"]    = res.val
  );
  
}


/***R

set.seed(123)
N <- 1000
M <- 2000

nedges  <- 1e3
source <- sample.int(N, nedges, replace = TRUE)
target <- sample.int(M, nedges, replace = TRUE)
values <- runif(nedges)

el <- new_Array(N, M, source - 1L, target - 1L, values)
# View(cbind(cbind(source, target) - 1, values))

# First, we need to get the same values
entries <- get_entries(el)
entries <- do.call(cbind, entries)
entries <- entries[order(entries[,3]),]
range(entries - cbind(source - 1, target - 1, values)[order(values),]) # 00


ans <- sapply(1:nedges, function(i) get_cell(el, source[i]-1, target[i]-1))
range(values - ans)

# The most popular values
i_max <- table(source)
i_max <- as.integer(names(i_max)[which.max(i_max)])
j_max <- table(target)
j_max <- as.integer(names(j_max)[which.max(j_max)])

# This has 3 values
unique(get_row(el, i_max - 1))
unique(get_col(el, j_max - 1))

# Comparing with a sparse matrix
library(Matrix)
x_sparse <- sparseMatrix(i = source, j = target, x = values, dims = c(N,M))
 
# All equal?
ans <- sapply(1:N - 1, get_row, x = el)
range(t(ans) - x_sparse)
mean(ans)

# Is it worth it?
microbenchmark::microbenchmark(
  Matrix = x_sparse[i_max,],
  Array  = get_row(el, i_max - 1)
)
# 
microbenchmark::microbenchmark(
  Matrix = x_sparse[,j_max],
  Array  = get_col(el, j_max - 1)
)
range(x_sparse[j_max,] - get_row(el, j_max - 1))

*/


