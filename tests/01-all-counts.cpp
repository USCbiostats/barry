// #include <vector>
// #include "../include/barry.hpp"
// #include "catch.hpp"

TEST_CASE("Network counts work", "[counts]") {
  
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
   nodeicov("gender") + nodeocov("gender") + nodecov("gender") + absdiff("age"))
   cat(sprintf("%.4f", unname(ans)), sep = ", ")
   */
  
  std::vector< double > gender = {2,2,2,2,1,1,2,1,1,2,2,2,2,1,1,1,1,1,2,2};
  std::vector< double > age    = {0.1476,0.2609,0.6445,0.0231,0.1808,0.7323,0.1051,0.7529,0.1255,0.1288,0.1147,0.6299,0.5036,0.9889,0.3125,0.6167,0.7872,0.3611,0.2492,0.0615};
  std::vector< unsigned int > source = {5,14,17,14,0,4,5,10,13,16,1,13,19,10,15,17,12,15,7,8,2,3,14,8,14,2,16,8,12,18,1,4,6,8,12,17,18,1,6,7,13,2,4,5,7,1,10,13,0,1,6,7,6,8,14,17,18};
  std::vector< unsigned int > target = {0,0,0,1,2,2,2,2,2,2,3,3,3,4,4,4,5,5,6,6,7,7,7,10,10,12,12,13,13,13,14,14,14,14,14,14,14,15,15,15,15,16,16,16,16,17,17,17,18,18,18,18,19,19,19,19,19};

  // As computed by ergm  
  std::vector< double > expected_counts = {57.0000, 2.0000, 2.0000, 86.0000, 75.0000, 25.0000, 11.0000, 0.1500, 111.8882, 107.2032, 24.0000, 85.0000, 82.0000, 167.0000, 15.6926};
  
  netcounters::Network net(20, 20, source, target);
  std::vector< std::vector< double > > vattrs(2);
  vattrs[0] = gender;
  vattrs[1] = age;
  net.data = new netcounters::NetworkData(vattrs);
  
  netcounters::NetStatsCounter counter(&net); 
  
  counter.add_counter(netcounters::counter_edges()); 
  counter.add_counter(netcounters::counter_mutual()); 
  counter.add_counter(netcounters::counter_isolates());
  counter.add_counter(netcounters::counter_istar2());
  counter.add_counter(netcounters::counter_ostar2());
  counter.add_counter(netcounters::counter_ttriads());
  counter.add_counter(netcounters::counter_ctriads());
  counter.add_counter(netcounters::counter_density());
  counter.add_counter(netcounters::counter_idegree15());
  counter.add_counter(netcounters::counter_odegree15()); 
  counter.add_counter(netcounters::counter_nodematch(0u));
  counter.add_counter(netcounters::counter_nodeicov(0u));
  counter.add_counter(netcounters::counter_nodeocov(0u));
  counter.add_counter(netcounters::counter_nodecov(0u));
  counter.add_counter(netcounters::counter_absdiff(1u, 1.0));
  
  
  std::vector< double > ans = counter.count_all();
  
  delete net.data;
  
  std::vector<double> margin(ans.size(), 0.01);
  print(ans);
  print(expected_counts);
  REQUIRE(vabsdiff(ans, expected_counts) < margin);
}

