#ifndef GEESE_FLOCK_BONES_HPP
#define GEESE_FLOCK_BONES_HPP 1

class Geese;

class Flock {
public:

    std::vector< Geese > dat;
    
    // Common components
    std::mt19937 *                 rengine  = nullptr;
    phylocounters::PhyloCounters * counters = nullptr;
    phylocounters::PhyloModel *    support  = nullptr;

    Flock() {};
    ~Flock() {};

    unsigned int add_data(
        std::vector< std::vector<unsigned int> > & annotations,
        std::vector< unsigned int > &              geneid,
        std::vector< int > &                       parent,
        std::vector< bool > &                      duplication
    );

    void init();
    
    // void add_geese(Geese x);
    phylocounters::PhyloCounters * counters_ptr();

    double likelihood_joint(const std::vector< double > & par, bool as_log = false);

};

#endif