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
    size_t nfunctions  = 0u;
    bool initialized = false;
    
    // Common components
    std::mt19937 rengine;
    PhyloModel model = PhyloModel();

    Flock() {};
    ~Flock() {};

    /**
     * @brief Add a tree to the flock
     * 
     * @param annotations see Geese::Geese.
     * @param geneid see Geese.
     * @param parent see Geese.
     * @param duplication see Geese.
     * @return size_t The number of tree in the model (starting from zero).
     */
    size_t add_data(
        std::vector< std::vector<size_t> > & annotations,
        std::vector< size_t > &              geneid,
        std::vector< int > &                 parent,
        std::vector< bool > &                duplication
    );

    /**
     * @brief Set the seed of the model
     * 
     * @param s Passed to the `rengine.seed()` member object.
     */
    void set_seed(const size_t & s);

    void init(size_t bar_width = BARRY_PROGRESS_BAR_WIDTH);
    
    // void add_geese(Geese x);
    PhyloCounters * get_counters();
    PhyloSupport *  get_support_fun();
    std::vector< std::vector< double > > * get_stats_support();
    std::vector< std::vector< double > > * get_stats_target();
    PhyloModel *  get_model();

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
        bool use_reduced_sequence = true,
        size_t ncores = 1u
    );

    /**
     * @name Information about the model 
     */
    ///@{
    size_t nfuns() const noexcept;
    size_t ntrees() const noexcept;
    std::vector< size_t > nnodes() const noexcept;
    std::vector< size_t > nleafs() const noexcept;
    size_t nterms() const;
    size_t support_size() const noexcept;
    std::vector< std::string > colnames() const;
    size_t parse_polytomies(
        bool verb = true,
        std::vector< size_t > * dist = nullptr
        ) const noexcept;  ///< Check polytomies and return the largest.
    void print() const;
    ///@}

    /**
     * @brief Access the i-th geese element
     * 
     * @param i Element to access
     * @param check_bounds When true, it will check bounds.
     * @return Geese*
     */
    Geese * operator()(size_t i, bool check_bounds = true);

};

#endif