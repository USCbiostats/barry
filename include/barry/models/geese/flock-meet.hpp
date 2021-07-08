#ifndef GEESE_FLOCK_MEET_HPP 
#define GEESE_FLOCK_MEET_HPP 1

// #include "flock-bones.hpp"

inline unsigned int Flock::add_data(
    std::vector< std::vector<unsigned int> > & annotations,
    std::vector< unsigned int > &              geneid,
    std::vector< int > &                       parent,
    std::vector< bool > &                      duplication
) {

    // Setting up the model
    if (dat.size() == 0u) {

        model.set_rengine(&this->rengine, false);
        model.set_keygen(keygen_full);
        model.store_psets();

    } else {

        if (annotations[0u].size() != nfuns())
            throw std::length_error("The number of functions in the new set of annotations does not match that of the first Geese.");

    }

    // Generating the Geese object
    dat.push_back(Geese(annotations, geneid, parent, duplication));

    if (dat.size() == 1u)
        this->nfunctions = dat[0].nfuns();
       
    return dat.size() - 1u;

}

inline void Flock::set_seed(const unsigned int & s) {
    this->rengine.seed(s);
}

inline void Flock::init(bool verb) {

    // For some strange reason, pointing to model during
    // the add_data function changes addresses once its out.
    for (auto& a : dat) {

        if (a.delete_support)
            delete a.model;

        a.model          = &model;
        a.delete_support = false;

        if (a.delete_rengine)
            delete a.rengine;

        a.rengine         = &rengine;
        a.delete_rengine  = false;
        
    }

    // Initializing the models.
    for (auto& d : dat) 
        d.init(verb);

    this->initialized = true;
    
}

inline phylocounters::PhyloCounters * Flock::get_counters() {

    if (dat.size() == 0u)
        throw std::logic_error("The flock has no data yet.");

    return this->model.get_counters();

}

inline phylocounters::PhyloSupport *  Flock::get_support() {
    return this->model.get_support();
}

inline phylocounters::PhyloModel *  Flock::get_model() {
    return &this->model;
}

inline double Flock::likelihood_joint(
    const std::vector< double > & par,
    bool as_log,
    bool use_reduced_sequence
) {

    INITIALIZED()

    double ans = as_log ? 0.0: 1.0;
    if (as_log) {

        for (auto& d : this->dat) 
            ans += d.likelihood(par, as_log, use_reduced_sequence);

    } else {

        for (auto& d : this->dat) 
            ans *= d.likelihood(par, as_log, use_reduced_sequence);
            
    }
    
    return ans;

}

inline unsigned int Flock::nfuns() const noexcept {

    return this->nfunctions;

}

inline unsigned int Flock::ntrees() const noexcept {

    return this->dat.size();

}

inline std::vector< unsigned int > Flock::nnodes() const noexcept {

    std::vector< unsigned int > res;
    res.reserve(this->ntrees());

    for (const auto& d : dat)
        res.push_back(d.nnodes());

    return res;
}

inline std::vector< unsigned int > Flock::nleafs() const noexcept {

    std::vector< unsigned int > res;
    res.reserve(this->ntrees());

    for (const auto& d : dat)
        res.push_back(d.nleafs());

    return res;

}

inline unsigned int Flock::nterms() const {

    INITIALIZED()
    return model.nterms() + this->nfuns();

}

inline unsigned int Flock::support_size() const noexcept {

    return this->model.support_size();

}

inline std::vector< std::string > Flock::colnames() const {

    return this->model.colnames();

}

inline unsigned int Flock::parse_polytomies(bool verb) const noexcept {

    unsigned int ans = 0;
    int i = 0;
    for (const auto & d : dat) {

        if (verb)
            printf_barry("Checking tree %i\n", i);

        unsigned int tmp = d.parse_polytomies(verb);

        if (tmp > ans)
            ans = tmp;
    }

    return ans;

}

inline Geese* Flock::operator()(unsigned int i, bool check_bounds)
{

    if (check_bounds && i >= ntrees())
        throw std::logic_error("Geese not found in the flock (out of range).");

    return &dat[i];

}

#endif