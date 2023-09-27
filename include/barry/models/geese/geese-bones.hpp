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

    return Array(i, j) == 9u;
    
}



// Hasher

inline std::vector< double > keygen_full(
    const PhyloArray & array,
    const PhyloCounterData * d
    ) {

    // Baseline data: nrows and columns
    std::vector< double > dat = {
        static_cast<double>(array.nrow()) * 100000 +
         static_cast<double>(array.ncol()),
         1000000.0, // state of the parent
         array.D_ptr()->duplication ? 1.0 : 0.0 // type of the parent
    };

    // State of the parent
    double pow10 = 1.0;
    for (bool i : array.D_ptr()->states) {
        dat[1u] += (i ? 1.0 : 0.0) * pow10;
        pow10 *= 10.0;
    }

    return dat;
    
}

inline bool vec_diff(
    const std::vector< size_t > & s,
    const std::vector< size_t > & a
) {

    for (size_t i = 0u; i < a.size(); ++i)
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
/**
 * @brief Class representing a phylogenetic tree model with annotations.
 * 
 * The `Geese` class represents a phylogenetic tree model with annotations. It
 * includes a total of `N + 1` nodes, the `+ 1` being the root node. The class
 * provides methods for initializing the model, calculating the likelihood,
 * simulating trees, and making predictions. 
 * 
 * The class includes shared objects within a `Geese` object, such as `rengine`,
 * `model`, `states`, `n_zeros`, `n_ones`, `n_dupl_events`, and `n_spec_events`.
 * It also includes information about the type of event, such as `etype_default`,
 * `etype_speciation`, `etype_duplication`, and `etype_either`.
 * 
 * The class provides constructors, a destructor, and methods for initializing
 * the model, inheriting support, calculating the sequence, calculating the
 * reduced sequence, calculating the likelihood, calculating the likelihood
 * exhaustively, getting probabilities, setting the seed, simulating trees,
 * parsing polytomies, getting observed counts, printing observed counts,
 * printing information about the GEESE, and making predictions.
 * 
 * @see Flock
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
    PhyloModel *        model   = nullptr;
    std::vector< std::vector< bool > > states;
    size_t n_zeros       = 0u; ///< Number of zeros
    size_t n_ones        = 0u; ///< Number of ones
    size_t n_dupl_events = 0u; ///< Number of duplication events
    size_t n_spec_events = 0u; ///< Number of speciation events
    ///@}

public:

    // Data
    size_t                       nfunctions;
    std::map< size_t, Node >     nodes;
    
    barry::MapVec_type< size_t > map_to_state_id;
    std::vector< std::vector< std::vector< size_t > > > pset_loc;    ///< Locations of columns

    // Tree-traversal sequence
    std::vector< size_t > sequence;
    std::vector< size_t > reduced_sequence;  

    // Admin-related objects
    bool initialized     = false;
    bool delete_rengine  = false;
    bool delete_support  = false;

    // Information about the type of event
    
    /***
     * @name Information about the type of event
     * @details
     * The type of event is stored in the `etype` member. The possible values
     * are `etype_default`, `etype_speciation`, `etype_duplication`, and
     * `etype_either`.
     * 
    */
    ///@{
    static const size_t etype_default     = 1ul;
    static const size_t etype_speciation  = 0ul;
    static const size_t etype_duplication = 1ul;
    static const size_t etype_either      = 2ul;
    ///@}

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
        std::vector< std::vector<size_t> > & annotations,
        std::vector< size_t > &              geneid,
        std::vector< int > &                 parent,
        std::vector< bool > &                duplication
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

    void init(size_t bar_width = BARRY_PROGRESS_BAR_WIDTH);

    void inherit_support(const Geese & model_, bool delete_support_ = false);

    // Node * operator()(size_t & nodeid);
    void calc_sequence(Node * n = nullptr);
    void calc_reduced_sequence();

    double likelihood(
        const std::vector< double > & par,
        bool as_log = false,
        bool use_reduced_sequence = true,
        BARRY_NCORES_ARG(= 1)
        );

    double likelihood_exhaust(const std::vector< double > & par);

    std::vector< double > get_probabilities() const;

    void set_seed(const size_t & s);
    std::vector< std::vector< size_t > > simulate(
        const std::vector< double > & par
        );

    /**
     * @name Information about the model 
     * @param verb When `true` it will print out information about the encountered
     * polytomies.
     */
    ///@{
    size_t nfuns() const noexcept;             ///< Number of functions analyzed
    size_t nnodes() const noexcept;            ///< Number of nodes (interior + leaf)
    size_t nleafs() const noexcept;            ///< Number of leaf
    size_t nterms() const;                     ///< Number of terms included
    size_t support_size() const noexcept;      ///< Number of unique sets of sufficient stats.
    std::vector< size_t > nannotations() const noexcept;      ///< Number of annotations.
    std::vector< std::string > colnames() const;     ///< Names of the terms in the model.
    size_t parse_polytomies(
        bool verb = true,
        std::vector< size_t > * dist = nullptr
        ) const noexcept;  ///< Check polytomies and return the largest.

    ///@}

    std::vector< std::vector<double> > observed_counts();
    void print_observed_counts();

    /**
     * @brief Prints information about the GEESE
     */
    void print() const;
    void print_nodes() const;


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
        const std::vector< size_t > & preorder
        );

    std::vector< std::vector< double > > predict_exhaust_backend(
        const std::vector< double > & par,
        const std::vector< size_t > & preorder
        );

    std::vector< std::vector< double > > predict_exhaust(
        const std::vector< double > & par
        );

    std::vector< std::vector< double > > predict_sim(
        const std::vector< double > & par,
        bool only_annotated       = false,
        size_t nsims        = 10000u
        );
    ///@}

    void init_node(Node & n);
    void update_annotations(
        size_t nodeid,
        std::vector< size_t > newann
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
    std::mt19937 *  get_rengine();
    PhyloCounters * get_counters();
    PhyloModel *    get_model();
    PhyloSupport *  get_support_fun();
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
    std::vector< size_t > get_annotated_nodes() const; ///< Returns the ids of the nodes with at least one annotation
    std::vector< size_t > get_annotations() const; ///< Returns the annotations of the nodes with at least one annotation

};

#endif
