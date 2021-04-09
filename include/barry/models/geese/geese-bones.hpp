#ifndef GEESE_BONES_HPP
#define GEESE_BONES_HPP 1

// #include <vector>
// #include <algorithm>
// #include <random>
// #include <stdexcept>

template<typename Ta, typename Tb>
inline std::vector< Ta > vector_caster(const std::vector< Tb > & x) {
    std::vector< Ta > ans;
    ans.reserve(x.size());
    for (auto i = x.begin(); i != x.end(); ++i)
        ans.push_back(static_cast< Ta >(*i));
    return ans;
}

#define INITIALIZED() if (!this->initialized) \
    throw std::logic_error("The model has not been initialized yet.");

// The same need to be locked
RULE_FUNCTION(rule_empty_free) {
    if (Array.get_cell(i, j) == 9u)
        return false;
    return true;
}



// Hasher
inline std::vector< double > keygen_full(
    const phylocounters::PhyloArray & array
    ) {

    // Baseline data: nrows and columns
    std::vector< double > dat = {
        (double) array.nrow(), (double) array.ncol()
    };

    // State of the parent
    for (bool i : array.data->states) {
        dat.push_back(i ? 1.0 : 0.0);
    }

    // Type of the parent
    dat.push_back(array.data->duplication ? 1.0 : 0.0);

    return dat;
}

inline bool vec_diff(
    const std::vector< unsigned int > & s,
    const std::vector< unsigned int > & a
) {

    for (unsigned int i = 0u; i < a.size(); ++i)
        if ((a.at(i) != 9u) && (a.at(i) != s.at(i)))
            return true;

    return false;
}

/**
 * @ingroup stat-models
 * @brief Annotated Phylo Model
 * @details A list of available terms for this model can be found in the
 * \ref counters-phylo section.
 *
 */
class Geese {
public:

    /**
     * @name Shared objects within a `Geese`
     * @details
     * Since users may start adding counters before initializing the PhyloModel
     * object, the object `counter` is initialized first.
     * 
     * While the member `support` has an `rengine`, since `Geese` can sample trees,
     * we have the option to keep it separate.
     * 
     */
    ///@{
    std::mt19937 *                     rengine  = nullptr;
    phylocounters::PhyloCounters *     counters = nullptr;
    phylocounters::PhyloModel *        support  = nullptr;
    std::vector< std::vector< bool > > states;
    ///@}

    // Data
    unsigned int                       nfunctions;
    std::map< unsigned int, Node >     nodes;
    barry::MapVec_type< unsigned int > map_to_nodes;

    // Tree-traversal sequence
    std::vector< unsigned int > sequence;
    std::vector< unsigned int > reduced_sequence;  

    // Admin-related objects
    bool initialized     = false;
    bool delete_rengine  = false;
    bool delete_counters = false;
    bool delete_support  = false;

    /**
     * @name Construct a new Geese object
     *
     * The model includes a total of `N + 1` nodes, the `+ 1` beign
     * the root node.
     *
     * @param annotations A vector of vectors with annotations. It should be of
     * length `k` (number of functions). Each vector should be of length `N`
     * (equal to the number of nodes, including interior). Possible values are
     * 0, 1, and 9.
     * @param geneid Id of the gene. It should be of length `N`.
     * @param parent Id of the parent gene. Also of length `N`
     */
    ///@{
    Geese();

    Geese(
        std::vector< std::vector<unsigned int> > & annotations,
        std::vector< unsigned int > &              geneid,
        std::vector< int> &                        parent,
        std::vector< bool > &                      duplication
        );

    // Copy constructor
    Geese(const Geese & model_, bool copy_data = true);
    
    // Constructor move
    Geese(Geese && x) noexcept;

    // Copy assignment
    Geese & operator=(const Geese & model_) = delete;

    // // Move assignment
    Geese & operator=(Geese && model_) noexcept = delete;

    ///@}

    ~Geese();

    void init();

    void inherit_support(const Geese & model_, bool delete_support_ = false);

    // Node * operator()(unsigned int & nodeid);
    void calc_sequence(Node * n = nullptr);
    void calc_reduced_sequence();

    double likelihood(
        const std::vector< double > & par,
        bool as_log = false,
        bool use_reduced_sequence = true
        );

    double likelihood_exhaust(const std::vector< double > & par);

    std::vector< double > get_probabilities() const;

    void set_seed(const unsigned int & s);
    std::vector< std::vector< unsigned int > > simulate(
        const std::vector< double > & par
        );

    /**
     * @name Information about the model 
     * 
     */
    ///@{
    unsigned int nfuns() const noexcept;
    unsigned int nnodes() const noexcept;
    unsigned int nleafs() const noexcept;
    unsigned int nterms() const;
    ///@}

    std::vector< std::vector<double> > observed_counts();
    void print_observed_counts();

    /**
     * @brief Calculate the conditional probability
     * 
     * @param p Vector of parameters.
     * @param res_prob Vector indicating each nodes' state probability.
     * 
     * @details When `res_prob` is specified, the function will attach
     * the member vector `probabilities` from the `Node`s objects. This
     * contains the probability that the ith node has either of the
     * possible states.
     * 
     * @return std::vector< double > Returns the posterior probability
     */
    std::vector< std::vector< double > > predict(
        const std::vector< double > & p,
        std::vector< std::vector< double > > * res_prob = nullptr,
        bool use_reduced_sequence = true,
        bool only_annotated       = false
        );

    void init_node(Node & n);
    void update_annotations(
        unsigned int nodeid,
        std::vector< unsigned int > newann
    );  

};

#endif
