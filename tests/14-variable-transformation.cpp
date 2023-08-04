#include "tests.hpp"

BARRY_TEST_CASE("Transformation of models", "[transformation]")
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

    // Reverse version (for checking transformation)
    std::vector< double > p0_rev = p0; 
    std::reverse(p0_rev.begin(), p0_rev.end());

    std::vector< double > p1_rev = p1; 
    std::reverse(p1_rev.begin(), p1_rev.end());

    std::vector< double > logs_expected = {65.31523, 51.38268};
    
    using namespace barry::counters::network;

    NetworkDense net(4, 4);
    net.set_data(new NetworkData({0,0,1,0}), true);
    NetSupport<NetworkDense> support(net); 
    
    // Preparing model  
    counter_edges<NetworkDense>(support.get_counters());
    counter_mutual<NetworkDense>(support.get_counters());
    counter_isolates<NetworkDense>(support.get_counters());
    counter_istar2<NetworkDense>(support.get_counters());
    counter_ostar2<NetworkDense>(support.get_counters());
    counter_ttriads<NetworkDense>(support.get_counters());
    counter_ctriads<NetworkDense>(support.get_counters());
    counter_density<NetworkDense>(support.get_counters());
    counter_idegree15<NetworkDense>(support.get_counters());
    counter_odegree15<NetworkDense>(support.get_counters());
    counter_nodematch<NetworkDense>(support.get_counters(),0u);
    
    rules_zerodiag<NetworkDense>(support.get_rules());
    
    // Getting the full support
    support.calc(); 
    
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

    // Preparing the model -----------------------------------------------------
    NetModel<NetworkDense> model2;
    model2.set_counters(support.get_counters());
    model2.set_rules(support.get_rules());
    model2.add_array(net);

    // This transformation reverses the order of the terms
    std::function<std::vector<double>(double *,size_t)> tfun = []
    (double * dat, size_t k) {

        std::vector< double > v(k);

        for (auto i = 0u; i < k; ++i)
            v[k - i - 1] = (*(dat + i));

        return v;

    };

    // Same variables
    std::vector< std::string > newnames = {
        "New nodematch",
        "New odegree15",
        "New idegree15",
        "New density",
        "New ctriads",
        "New ttriads",
        "New ostar2",
        "New istar2",
        "New isolates",
        "New mutual",
        "New edges"
    };

    model2.print();

    auto loglik0 = model2.likelihood_total(p0, true);
    auto loglik1 = model2.likelihood_total(p1, true);

    model2.set_transform_model(tfun, newnames);

    auto loglik0_rev = model2.likelihood_total(p0_rev, true);
    auto loglik1_rev = model2.likelihood_total(p1_rev, true);

    model2.print();

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(loglik0 == Approx(loglik0_rev).epsilon(0.00001));
    REQUIRE(loglik1 == Approx(loglik1_rev).epsilon(0.00001));
    #endif

    // Preparing the model Setting zeros and other mults -----------------------
    NetModel<NetworkDense> model3;
    model3.set_counters(support.get_counters());
    model3.set_rules(support.get_rules());
    model3.add_array(net);

    // This transformation reverses the order of the terms
    std::function<std::vector<double>(double *,size_t)> tfun2 = []
    (double * dat, size_t k) {

        // Removing the edge variable
        auto k_new = k - 1;
        std::vector< double > v(k_new);

        for (auto i = 0u; i < k_new; ++i)
            v[k_new - i - 1u] = (*(dat + i + 1u));

        // Rescaling mutual
        v[9u] /= 1.25;

        return v;

    };

    // Same variables
    std::vector< std::string > newnames2 = {
        "New nodematch",
        "New odegree15",
        "New idegree15",
        "New density",
        "New ctriads",
        "New ttriads",
        "New ostar2",
        "New istar2",
        "New isolates",
        "New mutual"
    };

    model3.print();

    // Updating the likelihoods to reveal the changes
    // - edges: 0
    // - mutuals: /1.25
    p0[0u] *= 0.0;
    p0[1u] /= 1.25;
    p1[0u] *= 0.0;
    p1[1u] /= 1.25;

    p0_rev.pop_back();
    p1_rev.pop_back();

    auto loglik0b = model3.likelihood_total(p0, true);
    auto loglik1b = model3.likelihood_total(p1, true);

    auto target = *model3.get_stats_target();

    model3.set_transform_model(tfun2, newnames2);

    auto loglik0b_rev = model3.likelihood_total(p0_rev, true);
    auto loglik1b_rev = model3.likelihood_total(p1_rev, true);

    auto targetb = *model3.get_stats_target();

    model3.print();

    #ifdef CATCH_CONFIG_MAIN
    REQUIRE(loglik0b == Approx(loglik0b_rev).epsilon(0.00001));
    REQUIRE(loglik1b == Approx(loglik1b_rev).epsilon(0.00001));
    #endif
   
}

