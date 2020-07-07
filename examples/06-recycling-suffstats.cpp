#include <Rcpp.h>
#include "../include/barry.hpp"
// #include "../include/lbarray-bones.hpp"
using namespace Rcpp; 


typedef std::vector< unsigned int > vuint;
 

// To get the support
// [[Rcpp::export]] 
List support (
    std::vector< int > N,
    std::vector< int >  M,
    std::vector< std::vector< double > > gender
  ) { 
  
  if (N.size() != M.size())
    throw std::length_error("N and M should have the same length");
  
  // Initializing the Binary array, and also the the suffstats counter
  
  List res(N.size());
  
  // Preparing model
  barry::Support<netcounters::Network, vuint> dat(N[0], M[0]);
  dat.add_counter(netcounters::edges);
  dat.add_counter(netcounters::mutual);
  dat.add_counter(netcounters::isolates);
  dat.add_counter(netcounters::istar2);
  dat.add_counter(netcounters::ostar2);
  dat.add_counter(netcounters::ttriads);
  dat.add_counter(netcounters::ctriads);
  dat.add_counter(netcounters::density);
  dat.add_counter(netcounters::idegree15);
  dat.add_counter(netcounters::odegree15);
  
  // Adding functions
  netcounters::NetCounter nodematchfem = netcounters::nodematch;
  nodematchfem.data = new vuint({0u});
  dat.add_counter(nodematchfem);  
  
  for (unsigned int i = 0u; i < N.size(); ++i) {
    
    netcounters::Network net(N[i],M[i]);
    net.data = new netcounters::NetworkData(gender[i]);
    
    // Need to pass the data evey time
    dat.reset(&net);
    
    // Generating the data
    dat.calc(0u, false); 
    
    // Generating the entries
    barry::Counts_type ans = dat.support.get_entries();
    
    List res_tmp(ans.size());
    for (unsigned int j = 0u; j < res_tmp.size(); ++j) {
      res_tmp[j] = List::create(
        _["x"] = ans.at(j).first,
        _["count"] = ans.at(j).second
      );
    }
    
    delete net.data;
    res[i] = clone(res_tmp);
    
  }
  
  // Final cleanup
  delete nodematchfem.data;
  
  return res;
}


/***R

set.seed(123)
N <- c(5, 4)
M <- N
gender <- lapply(N, function(i) sample.int(2, i, TRUE))

nets <- vector("list", length(N))
for (i in seq_along(N)) {
  
  nedges  <- 1
  source <- sample.int(N[i], nedges, replace = TRUE)
  target <- sample.int(M[i], nedges, replace = TRUE)
  values <- runif(nedges)
  
  # Removing diagonal
  idx <- which(source != target)
  source <- source[idx]
  target <- target[idx]
  values <- values[idx]
  
  mat <- matrix(0, nrow = N[i], ncol = M[i])
  mat[cbind(source, target)] <- 1L
  mat <- network::as.network(mat)
  network::set.vertex.attribute(mat, "gender", gender)
  
  nets[[i]] <- mat
  
}

ans <- support(N, M, gender)

stop()

ans0_5 <- t(sapply(ans[[1]], function(i) c(i$x, i$count)))
ans0_4 <- t(sapply(ans[[2]], function(i) c(i$x, i$count)))
str(ans0_4)
str(ans0_5)

stopifnot(sum(ans0_4[,ncol(ans0_4)]) == 2^(12))
stopifnot(sum(ans0_5[,ncol(ans0_5)]) == 2^(20))

ans1_05 <- ergm::ergm.allstats(
  nets[[1]] ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad + ctriad +
    density + idegree1.5 + odegree1.5 + nodematch("gender"), zeroobs = FALSE,
  force = TRUE, maxNumChangeStatVectors = 2^20)
ans1_04 <- ergm::ergm.allstats(
  nets[[2]] ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad + ctriad +
    density + idegree1.5 + odegree1.5 + nodematch("gender"), zeroobs = FALSE,
  force = TRUE, maxNumChangeStatVectors = 2^20)

ans1_04 <- cbind(ans1_04$statmat, w = ans1_04$weights)
ans1_05 <- cbind(ans1_05$statmat, w = ans1_05$weights)

# Sorting equally
sort_all <- function(x) {
  x[do.call(order, lapply(1:ncol(x), function(i) x[, i])),,drop=FALSE]
  
}
ans1 <- sort_all(ans1)
ans0 <- sort_all(ans0)
colnames(ans0) <- colnames(ans1)

# Checking as support
set.seed(123)
pars <- runif(ncol(ans1_04) - 1)
log(exp(pars %*% t(ans1_04[, -ncol(ans1_04)])) %*% cbind(ans1_04[,ncol(ans1_04)]))
log(exp(pars %*% t(ans0_4[, -ncol(ans0_4)])) %*% cbind(ans0_4[,ncol(ans0_4)]))

log(exp(pars %*% t(ans1_05[, -ncol(ans1_05)])) %*% cbind(ans1_05[,ncol(ans1_05)]))
log(exp(pars %*% t(ans0_5[, -ncol(ans0_5)])) %*% cbind(ans0_5[,ncol(ans0_5)]))



*/


