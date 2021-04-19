// #include "/opt/intel/oneapi/advisor/2021.2.0/include/advisor-annotate.h"
#include "../../include/barry/barry.hpp"
#include <chrono>

typedef std::chrono::time_point<std::chrono::system_clock> clocktime;

#define TOKENPASTE(a,b) a ## b
#define TIME_START(a) clocktime TOKENPASTE(a,_start) = std::chrono::system_clock::now();
#define TIME_END(a) clocktime TOKENPASTE(a,_end) = std::chrono::system_clock::now(); \
    std::chrono::duration< double, std::milli > \
    TOKENPASTE(a,_diff) = TOKENPASTE(a,_end) - TOKENPASTE(a,_start);
#define ELAPSED(a) TOKENPASTE(a,_diff).count()


int main(int argc, char* argv[]) {

    int nsamples;
    if (argc == 1u) {
        nsamples = 1;
    } else {
        nsamples = static_cast< int >(std::strtol(argv[1u], nullptr, 10));
    }

    using namespace barry::counters;

    auto params = {-.5, .5, .2, 2.0};
    double avg_build   = 0.0;
    double avg_loglike = 0.0;


    TIME_START(total)
    for (int i = 0u; i < nsamples; ++i) {

        TIME_START(build)
        network::Network a(5, 5);

        a.set_data(new network::NetworkData({1.0, 1.5, 0.0, 2.0}, true), true);

        network::NetModel netmod(1);
        
        network::rules_zerodiag(&netmod.rules);
        network::counter_edges(&netmod.counters);
        network::counter_ttriads(&netmod.counters);
        network::counter_idegree(&netmod.counters, {0});
        network::counter_nodeicov(&netmod.counters, 0u);
        // network::counter_absdiff(&netmod.counters, 0u);

        netmod.add_array(a);
        TIME_END(build)
        avg_build += (ELAPSED(build) / nsamples);
    
        TIME_START(likelihood)
        double ll = netmod.likelihood_total(params, true);
        TIME_END(likelihood)

        avg_loglike += (ELAPSED(likelihood) / nsamples);


        if ((i + 1) == nsamples) {

            #ifdef BARRY_USE_OMP
            printf("With OMP\n");
            #pragma omp parallel
            {
                #pragma omp single
                {
                    printf("Running with %i threads\n", omp_get_num_threads());
                    printf("Max number of threads: %i\n", omp_get_max_threads());
                }
            }
            #else
            printf("Without OMP\n");
            #endif

            printf(
                "Likelihood: %.4f, Time(build): %.6f, Time(ll): %.6f\n",
                ll, avg_build, avg_loglike
                );


        }

    }
    TIME_END(total)
    printf("\nElapsed time (%i reps): %.4f\n", nsamples, ELAPSED(total));

    return 0;

}