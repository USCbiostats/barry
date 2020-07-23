#include <Rcpp.h>
#include "../include/barry.hpp"
// #include "../include/lbarray-bones.hpp"
using namespace Rcpp; 


typedef std::vector< unsigned int > vuint;
 
// A more elaborated example
// [[Rcpp::export]]
SEXP new_model(
   std::vector< int > N,
   std::vector< int > M,
   const std::vector< std::vector< uint >> & source,
   const std::vector< std::vector< uint >> & target,
   const std::vector< std::vector< double >> & gender
) {
 
  // Checking sizes
  if (N.size() != M.size())
   throw std::length_error("N and M should have the same length");
  if (N.size() != gender.size())
   throw std::length_error("N and gender should have the same length");
  if (N.size() != source.size())
   throw std::length_error("N and source should have the same length");
  if (N.size() != target.size())
   throw std::length_error("N and target should have the same length");

  
  // Building a model
  Rcpp::XPtr< netcounters::NetModel > model(
    new netcounters::NetModel,
    true
  );
  // netcounters::NetModel model;
  
  netcounters::counter_edges(&model->counters);
  netcounters::counter_mutual(&model->counters);
  netcounters::counter_isolates(&model->counters);
  netcounters::counter_istar2(&model->counters);
  netcounters::counter_ostar2(&model->counters);
  netcounters::counter_ttriads(&model->counters);
  netcounters::counter_ctriads(&model->counters);
  netcounters::counter_density(&model->counters);
  netcounters::counter_idegree15(&model->counters);
  netcounters::counter_odegree15(&model->counters);
  netcounters::counter_nodematch(&model->counters, 0u);
  
  // Adding the arrays
  std::vector< unsigned int > ids(gender.size());
  for (unsigned int i = 0u; i < gender.size(); ++i) {
    
    // Creating the network object
    netcounters::Network net((uint) N[i], (uint) M[i], source[i], target[i]);
    net.set_data(new netcounters::NetworkData(gender[i]), true);
    
    // For now, we are just forcing the thing
    ids[i] = model->add_array(net, true);

    
  }

  return model;
 
}

// [[Rcpp::export]]
double loglik(SEXP model, const std::vector<double> & par, bool as_log = true) {

  Rcpp::XPtr< netcounters::NetModel > ptr(model);
  return ptr->likelihood_total(par, as_log);  
}

 

/***R

# Checking the counters --------------------------------------------------------
set.seed(173)
# adjm <- ergmito::rbernoulli(c(4, 4), .5)
adjm <- list(structure(c(0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 
                 0), .Dim = c(4L, 4L)), structure(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 
                                                    1, 0, 0, 0, 1, 1, 0), .Dim = c(4L, 4L)))
# gender <- list(gender=sample.int(2, 4, TRUE), gender=sample.int(2, 4, TRUE))
gender <- list(gender = c(1L, 2L, 1L, 1L), gender = c(2L, 1L, 1L, 2L))

el <- lapply(adjm, function(i) which(i!=0, arr.ind = TRUE) - 1L)
param <- runif(11)
model <- new_model(
  N = c(4L, 4L),
  M = c(4L, 4L),
  source = lapply(el, "[", i=, j=1),
  target = lapply(el, "[", i=, j=2),
  gender = gender
  )

ans1 <- loglik(model, param)

networks <- Map(function(a,b) network::network(a, vertex.attr = list(gender = b)), a = adjm, b=gender)

f <- ergmito::ergmito_formulae(
  networks ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad + ctriad +
  density + idegree1.5 + odegree1.5 + nodematch("gender")
)

f$loglik(param) - ans1

# How faster/slower?
microbenchmark::microbenchmark(
  f$loglik(runif(11)),
  loglik(model, runif(11)), times = 1000, unit = "relative"
)

# Consistency
nchecks <- 1e3
consistency <- replicate(
  nchecks, {
    p <- runif(11)
    f$loglik(p) - loglik(model, p)
  })

hist(consistency, xlim = c(-0.01, 0.01))

*/


