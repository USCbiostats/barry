#include <Rcpp.h>
#include "array.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
SEXP new_Array(
    int N, int M,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    const std::vector< double > & value
) {
  
  Rcpp::XPtr< Array > ptr(
    new Array((uint) N, (uint) M, source, target, value),
    true
  );
  
  return ptr;
  
}

// [[Rcpp::export]]
double get_cell(SEXP x, int i, int j) {
  Rcpp::XPtr< Array > xptr(x);
  return xptr->get_cell(i, j);
}

// Returning the i-th row
// [[Rcpp::export]]
NumericVector get_row(SEXP x, int i) {
  
  Rcpp::XPtr< Array > xptr(x);
  NumericVector ans(xptr->M, 0);
  const umap_int_cell * m = xptr->get_row(i);
  
  for (auto row = m->begin(); row != m->end(); ++row)
    ans[row->first] = row->second.value;
  
  return ans;
  
}

// [[Rcpp::export]]
NumericVector get_col(SEXP x, int i) {
  
  Rcpp::XPtr< Array > xptr(x);
  NumericVector ans(xptr->N, 0);
  const umap_int_cell_ptr * m = xptr->get_col(i);
  
  for (auto row = m->begin(); row != m->end(); ++row)
    ans[row->first] = row->second->value;
  
  return ans;
  
}

// [[Rcpp::export]]
int rm_cell(SEXP x, int i, int j) {
  
  Rcpp::XPtr< Array > xptr(x);
  xptr->rm_cell(i, j);
  return 0; 
  
}

// [[Rcpp::export]]
int insert_cell(SEXP x, int i, int j, double v) {
  Rcpp::XPtr< Array > xptr(x);
  xptr->insert_cell(i, j, v);
  return 0;
}

// [[Rcpp::export]]
int swap_cells(SEXP x, int i0, int j0, int i1, int j1) {
  Rcpp::XPtr< Array > xptr(x);
  xptr->swap_cells(i0, j0, i1, j1);
  return 0;
}

// [[Rcpp::export]]
int swap_rows(SEXP x, int i0, int i1) {
  Rcpp::XPtr< Array > xptr(x);
  xptr->swap_rows(i0, i1);
  return 0;
}

// [[Rcpp::export]]
int swap_cols(SEXP x, int j0, int j1) {
  Rcpp::XPtr< Array > xptr(x);
  xptr->swap_cols(j0, j1);
  return 0;
}

// [[Rcpp::export]]
int transpose(SEXP x) {
  Rcpp::XPtr< Array > xptr(x);
  xptr->transpose();
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

el <- new_Array(N, M, source - 1L, target - 1L, values)
# View(cbind(cbind(source, target) - 1, values))

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
# library(microbenchmark)
# microbenchmark(
#   Matrix = x_sparse[i_max,],
#   Array  = get_row(el, i_max - 1)
# )
# 
# microbenchmark(
#   Matrix = x_sparse[,j_max],
#   Array  = get_col(el, j_max - 1)
# )
# range(x_sparse[7425,] - get_row(el, 7424))

# Checking removing cells
set.seed(123)
M <- N <- 1000
nedges <- 500
source <- sample.int(N, nedges)
target <- sample.int(N, nedges)
values <- runif(nedges)

el2 <- new_Array(N, M, source - 1L, target - 1L, values)

# Getting matrix
mat0 <- sparseMatrix(i = source, j = target, x = values, dims = c(N, M))
mat1 <- as(t(sapply(1:N - 1, get_row, x = el2)), "dgCMatrix")

# Removing one
rm_cell(el2, 5, 0)
as(t(sapply(1:N - 1, get_row, x = el2)), "dgCMatrix")[1:5, 1:5]

# Adding two
insert_cell(el2, 0, 0, 1)
insert_cell(el2, 1, 1, 2)
as(t(sapply(1:N - 1, get_row, x = el2)), "dgCMatrix")[1:5, 1:5]

# Swapping cells
swap_cells(el2, 0, 0, 1, 1)
swap_cells(el2, N-1, N-1, N-2, N-2)
swap_cells(el2, 9, 3, 8, 2)
swap_cells(el2, 9, 3, 8, 2)
# insert_cell(el2, 1, 1, 2)
as(t(sapply(1:N - 1, get_row, x = el2)), "dgCMatrix")[1:5, 1:5]

# microbenchmark(
#   Array  = swap_rows(el2, 0, 1),
#   Matrix = {mat0[c(1,2),] <- mat0[c(2,1),]}
# )
# 
# microbenchmark(
#   Array  = swap_cols(el2, 0, 1),
#   Matrix = {mat0[,c(1,2)] <- mat0[,c(2,1)]}
# )

swap_rows(el2, 0, 1)
as(t(sapply(1:N - 1, get_row, x = el2)), "dgCMatrix")[1:5, 1:5]
as(sapply(1:M - 1, get_col, x = el2), "dgCMatrix")[1:5, 1:5] # Should be equivalent

# Checking transpose -----------------------------------------------------------
set.seed(123)
M <- N <- 2000
M <- M/2
ncells <- M
source <- sample.int(N, ncells)
target <- sample.int(M, ncells)
values <- runif(ncells)

el3 <- new_Array(N, M, source - 1L, target - 1L, values)
(m0 <- as(sapply(1:M - 1, get_col, x = el3), "dgCMatrix"))[1:5, 1:5]

# transpose(el3)
# (m1 <- as(sapply(1:N - 1, get_col, x = el3), "dgCMatrix"))


# Comparing 
microbenchmark::microbenchmark(
  m0 <<- t(m0),
  transpose(el3), times = 1e3,
  unit = "relative"
)

*/
