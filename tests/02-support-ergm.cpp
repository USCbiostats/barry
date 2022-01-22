#include "tests.hpp"

BARRY_TEST_CASE("Computing support for networks", "[support]")
{
    
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
    NetSupport<> support(net); 
    
    // Preparing model  
    
    counter_edges<>(support.get_counters());
    counter_mutual<>(support.get_counters());
    counter_isolates<>(support.get_counters());
    counter_istar2<>(support.get_counters());
    counter_ostar2<>(support.get_counters());
    counter_ttriads<>(support.get_counters());
    counter_ctriads<>(support.get_counters());
    counter_density<>(support.get_counters());
    counter_idegree15<>(support.get_counters());
    counter_odegree15<>(support.get_counters());
    counter_nodematch<>(support.get_counters(),0u);
    
    rules_zerodiag<>(support.get_rules());
    
    // Getting the full support
    support.calc(); 
    
    // Generating the entries
    auto ans = support.get_counts();
    
    // log(exp(p0 %*% t(ans[, -ncol(ans)])) %*% cbind(ans[,ncol(ans)]))
    std::vector< double > logs = {0.0, 0.0};
    auto nnets = 0u;
    size_t n_counters = support.get_counters()->size();
    size_t n_unique   = ans.size() / (n_counters + 1);

    for (auto n = 0u; n < n_unique; ++n)
    {
        
        // // Printing the first 10
        // if (n < 5)
        //   print(ans[n].first);
        
        double tmp0 = 0.0;
        double tmp1 = 0.0;
        nnets += ans[n * (n_counters + 1u)];
        
        // Now iterating through the parameters
        for (auto j = 0u; j < p0.size(); ++j)
        {

            tmp0 += p0[j] * ans[n * (n_counters + 1u) + j + 1u];
            tmp1 += p1[j] * ans[n * (n_counters + 1u) + j + 1u];

        }
        
        // Exp and weights
        logs[0u] += exp(tmp0) * ans[n * (n_counters + 1u)];
        logs[1u] += exp(tmp1) * ans[n * (n_counters + 1u)];

    }
    
    logs[0u] = log(logs[0u]);
    logs[1u] = log(logs[1u]);
    
    
    // std::vector< double > margin = {std::fabs(logs_expected[0u] * .001), std::fabs(logs_expected[1u] * .001)};
    std::cout << nnets << " networks." << std::endl;
    print(logs);
    print(logs_expected);

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE_THAT(logs, Catch::Approx(logs_expected).epsilon(0.001));
    #endif

    // Checking change statistics
    // - Triangle: 0-1-2-0
    // - Isolates: 3
    // - Extra dyad: 0-2
    // 0 1 1 0
    // 0 0 1 0
    // 1 0 0 0
    // 0 0 0 0
    net.clear();
    net(0, 1) = 1;
    net(0, 2) = 1;
    net(1, 2) = 1;
    net(2, 0) = 1;

    // Preparing the model
    NetModel<Network> model2;
    model2.set_counters(support.get_counters());
    model2.add_array(net);
    NetStatsCounter<> counter2(&net);
    counter2.set_counters(support.get_counters());

    // Counting stats
    auto counts0 = counter2.count_all();
    net(3, 0) = 1;
    auto counts1 = counter2.count_all();

    std::vector< double > diff(counts1.size(), 0.0);
    for (unsigned int i = 0; i < counts0.size(); ++i)
        diff[i] = counts1[i] - counts0[i];

    double logistic_prob0 = 1.0/
        (1.0 + std::exp(-barry::vec_inner_prod<double>(p0, diff)));

    double logistic_prob1 = model2.conditional_prob(net, p0, 3, 0);

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(logistic_prob0 == Approx(logistic_prob1).epsilon(0.001));
    #endif

    
}

