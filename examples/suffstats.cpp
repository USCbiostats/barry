#include <Rcpp.h>
#include "../include/barray.hpp"
#include "../include/lbarray-bones.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
SEXP suff_stats(const NumericMatrix & x) {
  
  Rcpp::XPtr< barray::StatsDB > xptr( 
    new barray::StatsDB(),
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
  
  Rcpp::XPtr< barray::StatsDB > xptr(x);
  
  // Now, getting the data
  barray::Counts_type ans = xptr->get_entries();
  
  List res(ans.size());
  for (unsigned int i = 0u; i < res.size(); ++i) {
    res[i] = List::create(
      _["x"]     = ans.at(i).first,
      _["count"] = ans.at(i).second
    );
  }
  
  return res;
  
}

// A more elaborated example
// [[Rcpp::export]]
NumericVector counter(
    int N, int M,
    const std::vector< uint > & source,
    const std::vector< uint > & target
) {
  
  // Initializing the Binary array, and also the the suffstats counter
  barray::BArray<bool> Array((uint) N, (uint) M, source, target);

  // Array.meta.set("undirected", true);
  
  // Creating the counter object; 
  barray::StatsCounter<bool> dat(&Array);
  
  // Adding functions 
  dat.add_counter(barray::counters::edges);
  dat.add_counter(barray::counters::mutual);
  dat.add_counter(barray::counters::isolates);
  dat.add_counter(barray::counters::istar2);
  dat.add_counter(barray::counters::ostar2);
  dat.add_counter(barray::counters::ttriads);
  dat.add_counter(barray::counters::ctriads);
  dat.add_counter(barray::counters::density);
  dat.add_counter(barray::counters::idegree15);
  dat.add_counter(barray::counters::odegree15);
  
  // Fingers crossed
  std::vector< double > ans = dat.count_all();
  
  return wrap(ans);
  
}

// To get the support
// [[Rcpp::export]] 
List support (
    int N, int M
) {
  
  // Initializing the Binary array, and also the the suffstats counter
  barray::Support<bool> dat(N, M);
  
  // Adding functions
  dat.add_counter(barray::counters::edges);
  dat.add_counter(barray::counters::mutual);
  dat.add_counter(barray::counters::isolates);
  dat.add_counter(barray::counters::istar2);
  dat.add_counter(barray::counters::ostar2);
  dat.add_counter(barray::counters::ttriads);
  dat.add_counter(barray::counters::ctriads);
  dat.add_counter(barray::counters::density);
  dat.add_counter(barray::counters::idegree15);
  dat.add_counter(barray::counters::odegree15);
  
  // Generating the data
  dat.calc();
  
  // Generating the entries
  barray::Counts_type ans = dat.support.get_entries();
  
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

nedges  <- N*20
source <- sample.int(N, nedges, replace = TRUE)
target <- sample.int(M, nedges, replace = TRUE)
values <- runif(nedges)

# Removing diagonal
idx <- which(source != target)
source <- source[idx]
target <- target[idx]
values <- values[idx]

# Adding a few mutual
# source <- c(source, target[1:20])
# target <- c(target, source[1:20])
# values <- c(values, 1:20)

el <- counter(N, M, source - 1L, target - 1L)

# Comparing with ergm
net <- network::as.edgelist(cbind(source, target), n = N)
net <- network::as.network(net)
microbenchmark::microbenchmark(
  # ergmito::count_stats(mat ~ edges + mutual + istar2 + ostar2 + ttriad),
  ergm = ergm::summary_formula(net ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad + ctriad + density + idegree1.5 + odegree1.5),
  barray = counter(N, M, source - 1L, target - 1L),
  unit = "relative",
  times = 100
)


ergm::summary_formula(
  net ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad + ctriad +
    density + idegree1.5 + odegree1.5
  ) - 
counter(N, M, source - 1L, target - 1L)


# stop()
set.seed(123)
N <- 4
M <- N

nedges  <- 1
source <- sample.int(N, nedges, replace = TRUE)
target <- sample.int(M, nedges, replace = TRUE)
values <- runif(nedges)

# Removing diagonal
idx <- which(source != target)
source <- source[idx]
target <- target[idx]
values <- values[idx]

mat <- matrix(0, nrow = N, ncol = M)
mat[cbind(source, target)] <- 1L
mat <- network::as.network(mat)

ans0 <- support(N, M)
ans0 <- t(sapply(ans0, function(i) c(i$x, i$count)))

microbenchmark::microbenchmark(
  barray = support(N, M),
  ergm   = ergm::ergm.allstats(mat ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad  + ctriad + density + idegree1.5 + odegree1.5, zeroobs = FALSE),
  times = 100,
  unit = "relative"
)

ans1 <- ergm::ergm.allstats(mat ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad + ctriad + density + idegree1.5 + odegree1.5, zeroobs = FALSE)
ans1 <- cbind(ans1$statmat, w = ans1$weights)

# Sorting equally
sort_all <- function(x) {
  x[do.call(order, lapply(1:ncol(x), function(i) x[, i])),,drop=FALSE]
  
}
ans1 <- sort_all(ans1)
ans0 <- sort_all(ans0)
colnames(ans0) <- colnames(ans1)

# Checking as support
set.seed(123)
pars <- runif(ncol(ans1) - 1)
log(exp(pars %*% t(ans1[, -ncol(ans1)])) %*% cbind(ans1[,ncol(ans1)]))
log(exp(pars %*% t(ans0[, -ncol(ans0)])) %*% cbind(ans0[,ncol(ans0)]))

# Checking range
range(ans1 - ans0)


ps <- ergmito::powerset(3)
ps <- sapply(ps, function(ps.) {
  d <- which(ps. != 0, arr.ind = TRUE) - 1
  counter(3, 3, d[,1], d[,2])
})
ps <- t(ps)
unique(ps[,-ncol(ps)])

colnames(ps) <- c(
  "edges",
  "mutual",
  "isolates",
  "istar2",
  "ostar2",
  "ttriad",
  "ctriad",
  # "count",
  "density",
  "idegree1.5",
  "odegree1.5"
)
ps
*/


