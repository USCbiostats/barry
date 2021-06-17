// #include <vector>
// #include "../include/barry.hpp"
// #include "catch.hpp"
// #include "tests.h"

TEST_CASE("Computing support for networks", "[support]") {
  
  // Reading large network
  /**
   set.seed(123)
   x <- ergmito::rbernoulli(4, .2)
   gender <- c(0,0,1,0)
   x <- network::network(x, vertex.attr = list(gender = gender))
   ans <- ergm::ergm.allstats(x ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad +
   ctriad + density + idegree1.5 + odegree1.5 + nodematch("gender"),maxNumChangeStatVectors = 2^20,
   zeroobs = FALSE)
   ans <- cbind(ans$statmat, w = ans$weights)
   set.seed(9988)
   p0 <- runif(11)
   p1 <- runif(11)
   
   l0 <- log(exp(p0 %*% t(ans[, -ncol(ans)])) %*% cbind(ans[,ncol(ans)]))
   l1 <- log(exp(p1 %*% t(ans[, -ncol(ans)])) %*% cbind(ans[,ncol(ans)]))
   cat(sprintf("%.5f", p0), sep = ", ")
   cat(sprintf("%.5f", p1), sep = ", ")
   cat(sprintf("%.5f", c(l0, l1)), sep = ", ")
   cat(sprintf("%.5f", colMeans(ans[, -12])), sep = ", ")
   */
  
  std::vector< double > p0 = {0.60959, 0.01940, 0.90534, 0.62935, 0.01958, 0.97016, 0.51485, 0.47980, 0.18479, 0.87739, 0.02286};
  std::vector< double > p1 = {0.34609, 0.81370, 0.79881, 0.96398, 0.48765, 0.13675, 0.47716, 0.11797, 0.13809, 0.69155, 0.07703};
  std::vector< double > logs_expected = {65.31523, 51.38268};
  
  using namespace barry::counters::network;

  Network net(4, 4);
  net.set_data(new NetworkData({0,0,1,0}), true);
  NetSupport support(net); 
  
  // Preparing model  
  
  counter_edges(support.get_counters());
  counter_mutual(support.get_counters());
  counter_isolates(support.get_counters());
  counter_istar2(support.get_counters());
  counter_ostar2(support.get_counters());
  counter_ttriads(support.get_counters());
  counter_ctriads(support.get_counters());
  counter_density(support.get_counters());
  counter_idegree15(support.get_counters());
  counter_odegree15(support.get_counters());
  counter_nodematch(support.get_counters(),0u);
  
  rules_zerodiag(support.get_rules());
  
  // Getting the full support
  support.calc(); 
  
  // Generating the entries
  barry::Counts_type ans = support.get_counts();
  
  // log(exp(p0 %*% t(ans[, -ncol(ans)])) %*% cbind(ans[,ncol(ans)]))
  std::vector< double > logs = {0.0, 0.0};
  auto nnets = 0u;

  for (auto n = 0u; n < ans.size(); ++n) {
    
    // // Printing the first 10
    // if (n < 5)
    //   print(ans[n].first);
    
    double tmp0 = 0.0;
    double tmp1 = 0.0;
    nnets += ans[n].second;
    
    // Now iterating through the parameters
    for (auto j = 0u; j < p0.size(); ++j) {
      tmp0 += p0[j] * ans[n].first[j];
      tmp1 += p1[j] * ans[n].first[j];
    }
    
    // Exp and weights
    logs[0u] += exp(tmp0) * ans[n].second;
    logs[1u] += exp(tmp1) * ans[n].second;
  }
  
  logs[0u] = log(logs[0u]);
  logs[1u] = log(logs[1u]);
  
  
  // std::vector< double > margin = {std::fabs(logs_expected[0u] * .001), std::fabs(logs_expected[1u] * .001)};
  std::cout << nnets << " networks." << std::endl;
  print(logs);
  print(logs_expected);
  REQUIRE_THAT(logs, Catch::Approx(logs_expected).epsilon(0.001));
  
}

