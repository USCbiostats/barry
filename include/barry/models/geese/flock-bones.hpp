#ifndef GEESE_FLOCK_BONES_HPP
#define GEESE_FLOCK_BONES_HPP 1

class Geese;

/**
 * @brief A Flock is a group of Geese
 * @details This object buils a model with multiple trees (Geese objects),
 * with all of these using the same PhyloModel object.
 * 
 */
class Flock {
public:

    std::vector< Geese > dat;
    
    // Common components
    std::mt19937              rengine;
    phylocounters::PhyloModel support = phylocounters::PhyloModel();

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

    double likelihood_joint(
        const std::vector< double > & par,
        bool as_log = false,
        bool use_likelihood_sequence = true
        );

};

#endif