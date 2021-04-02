#ifndef GEESE_FLOCK_MEET_HPP 
#define GEESE_FLOCK_MEET_HPP 1

#include "flock-bones.hpp"

inline unsigned int Flock::add_data(
    std::vector< std::vector<unsigned int> > & annotations,
    std::vector< unsigned int > &              geneid,
    std::vector< int > &                       parent,
    std::vector< bool > &                      duplication
) {

    // Setting up the model
    if (dat.size() == 0u) {

        support.set_rengine(&this->rengine, false);
        support.set_keygen(keygen_full);
        support.store_psets();

    }

    // Generating the Geese object
    dat.push_back(Geese(annotations, geneid, parent, duplication));
       
    return dat.size() - 1u;

}

inline void Flock::init() {

    

    // For some strange reason, pointing to support during
    // the add_data function changes addresses once its out.
    for (auto& a : dat) {
        if (a.delete_support)
            delete a.support;

        a.support         = &support;
        a.delete_support  = false;

        if (a.delete_rengine)
            delete a.rengine;

        a.rengine         = &rengine;
        a.delete_rengine  = false;

        if (a.delete_counters)
            delete a.counters;

        a.counters        = &support.counters;
        a.delete_counters = false;
    }

    // Initializing the models.
    for (auto& d : dat) 
        d.init();
    
}

inline phylocounters::PhyloCounters * Flock::counters_ptr() {
    if (dat.size() == 0u)
        throw std::logic_error("The flock has no data yet.");

    return &this->support.counters;
}

inline double Flock::likelihood_joint(
    const std::vector< double > & par,
    bool as_log,
    bool use_likelihood_sequence
    ) {

    double ans = as_log ? 0.0: 1.0;
    if (as_log) {
        for (auto& d : this->dat) {
            ans += std::log(d.likelihood(par, use_likelihood_sequence));
        }
    } else {
        for (auto& d : this->dat) {
            ans *= d.likelihood(par, use_likelihood_sequence);
        }
    }
    
    return ans;

}
#endif