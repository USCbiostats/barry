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
   x <- network::network(x, vertex.attr = list(gender = gender))
   cat(gender, sep=",", '\n')
   library(ergm)
   ans <- ergm::summary_formula(x ~ edges + mutual + isolates + istar(2) + ostar(2) + ttriad +
   ctriad + density + idegree1.5 + odegree1.5 + nodematch("gender"))
   cat(sprintf("%.4f", unname(ans)), sep = ", ")
   */
  
  std::vector< double > gender = {2,2,2,2,1,1,2,1,1,2,2,2,2,1,1,1,1,1,2,2};
  std::vector< unsigned int > source = {5,14,17,14,0,4,5,10,13,16,1,13,19,10,15,17,12,15,7,8,2,3,14,8,14,2,16,8,12,18,1,4,6,8,12,17,18,1,6,7,13,2,4,5,7,1,10,13,0,1,6,7,6,8,14,17,18};
  std::vector< unsigned int > target = {0,0,0,1,2,2,2,2,2,2,3,3,3,4,4,4,5,5,6,6,7,7,7,10,10,12,12,13,13,13,14,14,14,14,14,14,14,15,15,15,15,16,16,16,16,17,17,17,18,18,18,18,19,19,19,19,19};

  // As computed by ergm  
  std::vector< double > expected_counts = {57.0000, 2.0000, 2.0000, 86.0000, 75.0000, 25.0000, 11.0000, 0.1500, 111.8882, 107.2032, 24.0000};
  
  netcounters::Network net(20, 20, source, target);
  net.data = new netcounters::NetworkData(gender);
  
  barry::StatsCounter<netcounters::Network, std::vector< unsigned int >> counter(&net);
  
  counter.add_counter(&netcounters::edges);
  counter.add_counter(&netcounters::mutual);
  counter.add_counter(&netcounters::isolates);
  counter.add_counter(&netcounters::istar2);
  counter.add_counter(&netcounters::ostar2);
  counter.add_counter(&netcounters::ttriads);
  counter.add_counter(&netcounters::ctriads);
  counter.add_counter(&netcounters::density);
  counter.add_counter(&netcounters::idegree15);
  counter.add_counter(&netcounters::odegree15);
  
  netcounters::NetCounter nodematchfem = netcounters::nodematch;
  nodematchfem.data = new std::vector< unsigned int >({0u});
  counter.add_counter(&nodematchfem);
  
  std::vector< double > ans = counter.count_all();
  
  delete net.data;
  delete nodematchfem.data;
  
  std::vector<double> margin(ans.size(), 0.01);
  print(ans);
  print(expected_counts);
  REQUIRE(vabsdiff(ans, expected_counts) < margin);
}

