#include <Rcpp.h>
#include "array.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
SEXP new_EdgeList(
    int N, int M,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    const std::vector< double > & value
) {
  
  Rcpp::XPtr< EdgeList > ptr(
    new EdgeList((uint) N, (uint) M, source, target, value),
    true
  );
  
  return ptr;
  
}

// [[Rcpp::export]]
double get_cell(SEXP x, int i, int j) {
  Rcpp::XPtr< EdgeList > xptr(x);
  return xptr->get_cell(i, j);
}

// Returning the i-th row
// [[Rcpp::export]]
NumericVector get_row(SEXP x, int i) {
  
  Rcpp::XPtr< EdgeList > xptr(x);
  NumericVector ans(xptr->M, 0);
  const umap_int_cell * m = xptr->get_row(i);
  
  for (auto row = m->begin(); row != m->end(); ++row)
    ans[row->first] = row->second.value;
  
  return ans;
  
}

// [[Rcpp::export]]
NumericVector get_col(SEXP x, int i) {
  
  Rcpp::XPtr< EdgeList > xptr(x);
  NumericVector ans(xptr->N, 0);
  const umap_int_cell_ptr * m = xptr->get_col(i);
  
  for (auto row = m->begin(); row != m->end(); ++row)
    ans[row->first] = row->second->value;
  
  return ans;
  
}

// [[Rcpp::export]]
int rm_cell(SEXP x, int i, int j) {
  
  Rcpp::XPtr< EdgeList > xptr(x);
  xptr->rm_cell(i, j);
  return 0; 
  
}

// [[Rcpp::export]]
int insert_cell(SEXP x, int i, int j, double v) {
  Rcpp::XPtr< EdgeList > xptr(x);
  xptr->insert_cell(i, j, v);
  return 0;
}

/***R

set.seed(123)
N <- 1000
M <- 2000

nedges  <- 1e3
source <- sample.int(N, nedges, replace = TRUE)
target <- sample.int(M, nedges, replace = TRUE)
values <- runif(nedges)

el <- new_EdgeList(N, M, source - 1L, target - 1L, values)
View(cbind(cbind(source, target) - 1, values))

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
library(microbenchmark)
microbenchmark(
  Matrix = x_sparse[i_max,],
  Array  = get_row(el, i_max - 1)
)

microbenchmark(
  Matrix = x_sparse[,j_max],
  Array  = get_col(el, j_max - 1)
)
# range(x_sparse[7425,] - get_row(el, 7424))

# Checking removing cells
set.seed(123)
M <- N <- 10
source <- sample.int(N, 5)
target <- sample.int(N, 5)
values <- runif(5)

el2 <- new_EdgeList(N, M, source - 1L, target - 1L, values)

# Getting matrix
mat0 <- sparseMatrix(i = source, j = target, x = values, dims = c(N, M))
mat1 <- as(t(sapply(1:N - 1, get_row, x = el2)), "dgCMatrix")

# Removing one
rm_cell(el2, 5, 0)
as(t(sapply(1:N - 1, get_row, x = el2)), "dgCMatrix")

# Adding two
insert_cell(el2, 0, 0, 1)
insert_cell(el2, 1, 1, 2)
as(t(sapply(1:N - 1, get_row, x = el2)), "dgCMatrix")

*/
