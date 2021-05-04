// #include <vector>
// #include "../include/barry.hpp"
// #include "catch.hpp"
// #include "tests.h"

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
  model.store_psets(); // Need this for sampling
  
  // Preparing model  
  
  netcounters::counter_edges(model.get_counters());
  netcounters::counter_mutual(model.get_counters());
  netcounters::counter_isolates(model.get_counters());
  netcounters::counter_istar2(model.get_counters());
  netcounters::counter_ostar2(model.get_counters());
  netcounters::counter_ttriads(model.get_counters());
  netcounters::counter_ctriads(model.get_counters());
  netcounters::counter_density(model.get_counters());
  netcounters::counter_idegree15(model.get_counters());
  netcounters::counter_odegree15(model.get_counters());
  netcounters::counter_nodematch(model.get_counters(), 0u);
  
  // Adding rules
  netcounters::rules_zerodiag(model.get_rules());
  
  // Adding the network to the model
  model.add_array(net);
  
  model.set_seed(1231);
  std::vector< double > p2 = {-1.0, 2.0, 0.0, 0.0, 0.0, 1.0*0, 0.0, 0.0, 0.0, 0.0, 0.0};
  netcounters::Network net0 = model.sample(0, p2);
  netcounters::Network net1 = model.sample(0, p2);
  netcounters::Network net2 = model.sample(0, p2);

  std::cout << "Printing networks" << std::endl;
  std::cout << "Net 0" << std::endl;
  net0.print();
  std::cout << "Net 1" << std::endl;
  net1.print();
  std::cout << "Net 2" << std::endl;
  net2.print();


  std::vector< double > logs0(2);
  std::vector< double > logs1(2);
  logs0[0u] = model.likelihood_total(p0, true); 
  logs0[1u] = model.likelihood_total(p1, true);
  logs1[0u] = model.likelihood(p0, net, 0, true); 
  logs1[1u] = model.likelihood(p1, net, 0, true);

  std::vector< double > margin = {0.00001, 0.00001};

  std::cout << "Printing logs: " << std::endl;
  print(logs0);
  print(logs1);
  print(logs_expected);
  REQUIRE(vabsdiff(logs0, logs_expected) < margin);
  REQUIRE(vabsdiff(logs1, logs_expected) < margin);
  
}

