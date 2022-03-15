#ifndef DEFM_BONES_HPP
#define DEFM_BONES_HPP 1

// #include <vector>
// #include <algorithm>
// #include <random>
// #include <stdexcept>

#define INITIALIZED() if (!this->initialized) \
    throw std::logic_error("The model has not been initialized yet.");

// The same need to be locked
RULE_FUNCTION(rule_empty_free) {

    return Array(i, j) == 9u;
    
}

// Hasher
inline std::vector< double > keygen_full(
    const phylocounters::PhyloArray & array
    ) {

    // Baseline data: nrows and columns
    std::vector< double > dat = {
        static_cast<double>(array.nrow()) * 100000 +
         static_cast<double>(array.ncol()),
         0.0
    };

    // State of the parent
    unsigned int count = 0u;
    for (bool i : array.D()->states) {
        dat[dat.size() - 1u] += (i ? 1.0 : 0.0) * pow(10, static_cast<double>(count));
        count++;
    }

    // Type of the parent
    dat.push_back(array.D()->duplication ? 1.0 : 0.0);

    return dat;
}

inline bool vec_diff(
    const std::vector< unsigned int > & s,
    const std::vector< unsigned int > & a
) {

    for (unsigned int i = 0u; i < a.size(); ++i)
        if ((a[i] != 9u) && (a[i] != s[i]))
            return true;

    return false;
}

class Flock;

/**
 * @ingroup stat-models
 * @brief Annotated Phylo Model
 * @details A list of available terms for this model can be found in the
 * \ref counters-phylo section.
 *
 */
class Geese {
    friend Flock;
private:

    /**
     * @name Shared objects within a `Geese`
     * @details
     * Since users may start adding counters before initializing the PhyloModel
     * object, the object `counter` is initialized first.
     * 
     * While the member `model` has an `rengine`, since `Geese` can sample trees,
     * we have the option to keep it separate.
     * 
     */
    ///@{
    std::mt19937 *                     rengine = nullptr;
    phylocounters::PhyloModel *        model   = nullptr;
    std::vector< std::vector< bool > > states;
    unsigned int n_zeros       = 0u; ///< Number of zeros
    unsigned int n_ones        = 0u; ///< Number of ones
    unsigned int n_dupl_events = 0u; ///< Number of duplication events
    unsigned int n_spec_events = 0u; ///< Number of speciation events
    ///@}

public:

    // Data
    unsigned int                       nfunctions;
    std::map< unsigned int, Node >     nodes;
    barry::MapVec_type< unsigned int > map_to_nodes;
    std::vector< std::vector< std::vector< size_t > > > pset_loc;    ///< Locations of columns

    // Tree-traversal sequence
    std::vector< unsigned int > sequence;
    std::vector< unsigned int > reduced_sequence;  

    // Admin-related objects
    bool initialized     = false;
    bool delete_rengine  = false;
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
     * @param duplication Logical scalar indicating the type of event (true:
     * duplication, false: speciation.)
     * 
     * @details 
     * The ordering of the entries does not matter. Passing the nodes in post
     * order or not makes no difference to the constructor.
     */
    ///@{
    Geese();

    Geese(
        std::vector< std::vector<unsigned int> > & annotations,
        std::vector< unsigned int > &              geneid,
        std::vector< int > &                       parent,
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

    void init(unsigned int bar_width = BARRY_PROGRESS_BAR_WIDTH);

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
     * @param verb When `true` it will print out information about the encountered
     * polytomies.
     */
    ///@{
    unsigned int nfuns() const noexcept;             ///< Number of functions analyzed
    unsigned int nnodes() const noexcept;            ///< Number of nodes (interior + leaf)
    unsigned int nleafs() const noexcept;            ///< Number of leaf
    unsigned int nterms() const;                     ///< Number of terms included
    unsigned int support_size() const noexcept;      ///< Number of unique sets of sufficient stats.
    std::vector< unsigned int > nannotations() const noexcept;      ///< Number of annotations.
    std::vector< std::string > colnames() const;     ///< Names of the terms in the model.
    unsigned int parse_polytomies(
        bool verb = true,
        std::vector< size_t > * dist = nullptr
        ) const noexcept;  ///< Check polytomies and return the largest.

    ///@}

    std::vector< std::vector<double> > observed_counts();
    void print_observed_counts();

    /**
     * @brief Prints information about the DEFM
     */
    void print() const;

    /**
     * @name Geese prediction
     * @brief Calculate the conditional probability
     * 
     * @param par Vector of parameters (terms + root).
     * @param res_prob Vector indicating each nodes' state probability.
     * @param leave_one_out When `true`, it will compute the predictions using
     * leave-one-out, thus the prediction will be repeated nleaf times.
     * @param only_annotated When `true`, it will make the predictions only
     * on the induced sub-tree with annotated leafs.
     * @param use_reduced_sequence  Passed to the `likelihood` method.
     * @param preorder For the tree traversal.
     * 
     * @details When `res_prob` is specified, the function will attach
     * the member vector `probabilities` from the `Node`s objects. This
     * contains the probability that the ith node has either of the
     * possible states.
     * 
     * @return std::vector< double > Returns the posterior probability
     */
    ///@{
    std::vector< std::vector< double > > predict(
        const std::vector< double > & par,
        std::vector< std::vector< double > > * res_prob = nullptr,
        bool leave_one_out        = false,
        bool only_annotated       = false,
        bool use_reduced_sequence = true
        );
    
    std::vector< std::vector<double> > predict_backend(
        const std::vector< double > & par,
        bool use_reduced_sequence,
        const std::vector< uint > & preorder
        );

    std::vector< std::vector< double > > predict_exhaust_backend(
        const std::vector< double > & par,
        const std::vector< uint > & preorder
        );

    std::vector< std::vector< double > > predict_exhaust(
        const std::vector< double > & par
        );

    std::vector< std::vector< double > > predict_sim(
        const std::vector< double > & par,
        bool only_annotated       = false,
        unsigned int nsims        = 10000u
        );
    ///@}

    void init_node(Node & n);
    void update_annotations(
        unsigned int nodeid,
        std::vector< unsigned int > newann
    );

    /**
     * @name Non-const pointers to shared objects in `Geese`
     * 
     * @details These functions provide direct access to some member
     * objects that are shared by the nodes within `Geese`.
     * 
     * @return `get_rengine()` returns the Pseudo-RNG engine used.
     * @return `get_counters()` returns the vector of counters used.
     * @return `get_model()` returns the `Model` object used.
     * @return `get_support_fun()` returns the computed support of the model.
     */
    ///@{
    std::mt19937 *                     get_rengine();
    phylocounters::PhyloCounters *     get_counters();
    phylocounters::PhyloModel *        get_model();
    phylocounters::PhyloSupport *      get_support_fun();
    ///@}
    
    /**
     * @brief Powerset of a gene's possible states
     * @details This list of vectors is used throughout `Geese`. It lists
     * all possible combinations of functional states for any gene. Thus,
     * for `P` functions, there will be `2^P` possible combinations.
     * 
     * @return std::vector< std::vector< bool > > of length `2^P`.
     */
    std::vector< std::vector< bool > > get_states() const;  
    std::vector< unsigned int >        get_annotated_nodes() const; ///< Returns the ids of the nodes with at least one annotation

};

#endif
