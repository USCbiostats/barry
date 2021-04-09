#ifndef GEESE_FLOCK_BONES_HPP
#define GEESE_FLOCK_BONES_HPP 1

class Geese;

/**
 * @ingroup stat-models
 * @brief A Flock is a group of Geese
 * @details This object buils a model with multiple trees (Geese objects),
 * with all of these using the same PhyloModel object. Available counters
 * (terms) can be found in \ref counter-phylo.
 * 
 */
class Flock {
public:

    std::vector< Geese > dat;
    unsigned int         nfunctions  = 0u;
    bool                 initialized = false;
    
    // Common components
    std::mt19937              rengine;
    phylocounters::PhyloModel support = phylocounters::PhyloModel();

    Flock() {};
    ~Flock() {};

    /**
     * @brief Add a tree to the flock
     * 
     * @param annotations see Geese::Geese.
     * @param geneid see Geese.
     * @param parent see Geese.
     * @param duplication see Geese.
     * @return unsigned int The number of tree in the model (starting from zero).
     */
    unsigned int add_data(
        std::vector< std::vector<unsigned int> > & annotations,
        std::vector< unsigned int > &              geneid,
        std::vector< int > &                       parent,
        std::vector< bool > &                      duplication
    );

    /**
     * @brief Set the seed of the model
     * 
     * @param s Passed to the `rengine.seed()` member object.
     */
    void set_seed(const unsigned int & s);

    void init();
    
    // void add_geese(Geese x);
    phylocounters::PhyloCounters * counters_ptr();

    /**
     * @brief Returns the joint likelihood of the model
     * 
     * @param par Vector of model parameters.
     * @param as_log When `true` it will return the value as log.
     * @param use_reduced_sequence When `true` (default) will compute the
     * likelihood using the reduced sequence, which is faster.
     * @return double 
     */
    double likelihood_joint(
        const std::vector< double > & par,
        bool as_log = false,
        bool use_reduced_sequence = true
    );

    /**
     * @name Information about the model 
     */
    ///@{
    unsigned int nfuns() const noexcept;
    unsigned int ntrees() const noexcept;
    std::vector< unsigned int > nnodes() const noexcept;
    std::vector< unsigned int > nleafs() const noexcept;
    unsigned int nterms() const;
    ///@}

};

#endif