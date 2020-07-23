#include <vector>
#include "../include/barry.hpp"
#include "catch.hpp"
#include "tests.h"

TEST_CASE("Computing support for networks (with Model)", "[support w model]") {
  
  // Reading large network
  /**
   set.seed(123)
   (x <- ergmito::rbernoulli(4, .2))
   
   
   gender <- c(0,0,1,0)
   x <- network::network(x, vertex.attr = list(gender = gender))
   ans <- ergm::ergm.allstats(x ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad +
   ctriad + density + idegree1.5 + odegree1.5 + nodematch("gender"),maxNumChangeStatVectors = 2^20,
   zeroobs = FALSE)
   ans <- cbind(ans$statmat, w = ans$weights)
   set.seed(9988)
   p0 <- runif(11)
   p1 <- runif(11)
   
   ans_ergmito <- ergmito::ergmito_formulae(
   x ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad +
   ctriad + density + idegree1.5 + odegree1.5 + nodematch("gender")
   )
   
   sufstat <- ergm::summary_formula(x ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad +
   ctriad + density + idegree1.5 + odegree1.5 + nodematch("gender"))
   
   l0 <- p0 %*% sufstat - log(exp(p0 %*% t(ans[, -ncol(ans)])) %*% cbind(ans[,ncol(ans)]))
   l1 <- p1 %*% sufstat - log(exp(p1 %*% t(ans[, -ncol(ans)])) %*% cbind(ans[,ncol(ans)]))
   cat(sprintf("%.5f", p0), sep = ", ")
   cat(sprintf("%.5f", p1), sep = ", ")
   cat(sprintf("%.5f", c(l0, l1)), sep = ", ")
   cat(sprintf("%.5f", colMeans(ans[, -12])), sep = ", ")
   */
  
  std::vector< double > p0 = {0.60959, 0.01940, 0.90534, 0.62935, 0.01958, 0.97016, 0.51485, 0.47980, 0.18479, 0.87739, 0.02286};
  std::vector< double > p1 = {0.34609, 0.81370, 0.79881, 0.96398, 0.48765, 0.13675, 0.47716, 0.11797, 0.13809, 0.69155, 0.07703};
  std::vector< double > logs_expected = {-61.79280, -48.59951};
  
  // phylocounters::PhyloArray node(4, 2);
  netcounters::Network net(4, 4, {2}, {3});
  net.set_data(new netcounters::NetworkData({0,0,1,0}), true);
  
  netcounters::NetModel model;
  
  // Preparing model  
  
  netcounters::counter_edges(&model.counters);
  netcounters::counter_mutual(&model.counters);
  netcounters::counter_isolates(&model.counters);
  netcounters::counter_istar2(&model.counters);
  netcounters::counter_ostar2(&model.counters);
  netcounters::counter_ttriads(&model.counters);
  netcounters::counter_ctriads(&model.counters);
  netcounters::counter_density(&model.counters);
  netcounters::counter_idegree15(&model.counters);
  netcounters::counter_odegree15(&model.counters);
  netcounters::counter_nodematch(&model.counters,0u);
  
  // Adding the network to the model
  model.add_array(net);
  
  std::vector< double > logs = {0.0, 0.0};
  logs[0u] = model.likelihood_total(p0, true);
  logs[1u] = model.likelihood_total(p1, true);
  
  std::vector< double > margin = {0.00001, 0.00001};
  print(logs);
  print(logs_expected);
  REQUIRE(vabsdiff(logs, logs_expected) < margin);
  
}

