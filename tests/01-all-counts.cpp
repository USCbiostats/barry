// #include <vector>
// #include "../include/barry.hpp"
// #include "catch.hpp"
#include "tests.hpp"

using namespace barry::counters::network;

BARRY_TEST_CASE("Network counts work", "[counts]") {
  
  // Reading large network
  /**
   set.seed(123); x <- ergmito::rbernoulli(20, .2)
   x[c(10, 12), ] <- 0
   x[, c(10, 12)] <- 0
   idx <- which(x != 0, arr.ind = TRUE) - 1
   # Ego
   cat(idx[,1], sep = ",")
   cat(idx[,2], sep = ",")
   
   set.seed(2244)
   gender <- sample.int(2, 20, TRUE)
   age    <- runif(20)
   x <- network::network(x, vertex.attr = list(gender = gender, age = age))
   cat(gender, sep=",", '\n')
   cat(sprintf("%.4f", age), sep=",", '\n')
   library(ergm)
   ans <- ergm::summary_formula(x ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad +
   ctriad + density + idegree1.5 + odegree1.5 + nodematch("gender") +
   nodeicov("gender") + nodeocov("gender") + nodecov("gender") +
   absdiff("age") + diff("age", pow=1) + diff("age", pow=2) + diff("age", pow=3) +
   idegree(1:4) + odegree(1:4))
   cat(sprintf("%.4f", unname(ans)), sep = ", ")
   */
  
  std::vector< double > gender = {2,2,2,2,1,1,2,1,1,2,2,2,2,1,1,1,1,1,2,2};
  std::vector< double > age    = {0.1476,0.2609,0.6445,0.0231,0.1808,0.7323,0.1051,0.7529,0.1255,0.1288,0.1147,0.6299,0.5036,0.9889,0.3125,0.6167,0.7872,0.3611,0.2492,0.0615};
  std::vector< size_t > source = {5,14,17,14,0,4,5,10,13,16,1,13,19,10,15,17,12,15,7,8,2,3,14,8,14,2,16,8,12,18,1,4,6,8,12,17,18,1,6,7,13,2,4,5,7,1,10,13,0,1,6,7,6,8,14,17,18};
  std::vector< size_t > target = {0,0,0,1,2,2,2,2,2,2,3,3,3,4,4,4,5,5,6,6,7,7,7,10,10,12,12,13,13,13,14,14,14,14,14,14,14,15,15,15,15,16,16,16,16,17,17,17,18,18,18,18,19,19,19,19,19};

  // As computed by ergm  
  std::vector< double > expected_counts = {
    57.0000, 2.0000, 2.0000, 86.0000, 75.0000, 25.0000, 11.0000, 0.1500,
    111.8882, 107.2032, 24.0000, 85.0000, 82.0000, 167.0000, 15.6926,
    -0.7202, 7.4380, -0.4130,
    1.0000, 4.0000, 6.0000, 3.0000, 2.0000, 3.0000, 6.0000, 4.0000
    };
  
  std::cout << "Creating the network...";
  Network net(20, 20, source, target);
  std::cout << "done." << std::endl;
  
  std::cout << "Preparating the model...";
  std::vector< std::vector< double > > vattrs(2);
  vattrs[0] = gender;
  vattrs[1] = age;
  net.set_data(new NetworkData(vattrs), true);
  
  NetStatsCounter<Network> counter(&net); 

  counter_edges<Network>(counter.get_counters());
  counter_mutual<Network>(counter.get_counters());
  counter_isolates<Network>(counter.get_counters());
  counter_istar2<Network>(counter.get_counters());
  counter_ostar2<Network>(counter.get_counters());
  counter_ttriads<Network>(counter.get_counters());
  counter_ctriads<Network>(counter.get_counters());
  counter_density<Network>(counter.get_counters());
  counter_idegree15<Network>(counter.get_counters());
  counter_odegree15<Network>(counter.get_counters());
  counter_nodematch<Network>(counter.get_counters(), 0u);
  counter_nodeicov<Network>(counter.get_counters(), 0u);
  counter_nodeocov<Network>(counter.get_counters(), 0u);
  counter_nodecov<Network>(counter.get_counters(), 0u);
  counter_absdiff<Network>(counter.get_counters(), 1u, 1.0);
  counter_diff<Network>(counter.get_counters(), 1u, 1.0);
  counter_diff<Network>(counter.get_counters(), 1u, 2.0);
  counter_diff<Network>(counter.get_counters(), 1u, 3.0);
  counter_idegree<Network>(counter.get_counters(), {1u,2u,3u,4u});
  counter_odegree<Network>(counter.get_counters(), {1u,2u,3u,4u});
  std::cout << "done." << std::endl;
  
  std::cout << "Starting the count...";
  std::vector< double > ans = counter.count_all();
  std::cout << "done." << std::endl;
  
  // delete net.data;
  // std::vector<double> margin(ans.size(), 0.001);
  std::cout << "Observed counts: " << std::endl;
  // counter.EmptyArray.print();
  print(ans);
  
  std::cout << "Expected counts: " << std::endl;
  // net.print();
  print(expected_counts);
  
  #ifdef CATCH_CONFIG_MAIN
  REQUIRE_THAT(expected_counts, Catch::Approx(ans).epsilon(.001));
  #endif
}

