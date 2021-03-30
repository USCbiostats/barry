#ifndef GEESE_FLOCK_MEET_HPP 
#define GEESE_FLOCK_MEET_HPP 1

inline unsigned int Flock::add_data(
    std::vector< std::vector<unsigned int> > & annotations,
    std::vector< unsigned int > &              geneid,
    std::vector< int > &                       parent,
    std::vector< bool > &                      duplication
) {

    // Generating the Geese object
    dat.push_back(Geese(annotations, geneid, parent, duplication));
    unsigned int i = dat.size() - 1;

    // We proceede depending on whether the support has already been initialized
    if (i == 0u) {
        dat[i].set_support(new phylocounters::PhyloModel(), true);
    } else
        dat[i].inherit_support(dat[0u], false);
    
    // Copying random number generation
    if (i != 0u) {
        dat[i].rengine = dat[0u].rengine;
        dat[i].delete_rengine = false;
    }

    return i;

}

inline phylocounters::PhyloCounters * Flock::counters_ptr() {
    if (dat.size() == 0u)
        throw std::logic_error("The flock has no data yet.");

    return &(dat[0u].counters);
}

inline double Flock::likelihood_joint(const std::vector< double > & par, bool as_log) {

    double ans = as_log ? 1.0: 0.0;
    if (as_log) {
        for (auto& d : this->dat) {
            ans += std::log(d.likelihood(par));
        }
    } else {
        for (auto& d : this->dat) {
            ans *= d.likelihood(par);
        }
    }
    
    return ans;

}
#endif