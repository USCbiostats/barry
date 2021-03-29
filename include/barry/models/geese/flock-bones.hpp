#ifndef GEESE_FLOCK_BONES_HPP
#define GEESE_FLOCK_BONES_HPP 1

class Geese;

class Flock {
public:

    std::vector< Geese > dat;

    Flock() {};
    ~Flock() {};

    unsigned int add_data(
        std::vector< std::vector<unsigned int> > & annotations,
        std::vector< unsigned int > &              geneid,
        std::vector< int > &                       parent,
        std::vector< bool > &                      duplication
    );
    
    // void add_geese(Geese x);
    phylocounters::PhyloCounters * counters_ptr();

};

#endif