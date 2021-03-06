#ifndef APHYLOMODEL_BONES_HPP
#define APHYLOMODEL_BONES_HPP 1

#define INITIALIZED() if (!this->initialized) \
    throw std::logic_error("The model has not been initialized yet.");

// The same need to be locked
RULE_FUNCTION(rule_empty_free) {
    if (Array->get_cell(i, j) == 9u)
        return false;
    return true;
}

// using namespace phylocounters;


template<typename T1, typename T2>
std::vector< T1 > caster(const std::vector< T2 > & vec) {

    std::vector< T1 > ans;
    ans.reserve(vec.size());

    for (auto& i : vec) {
        ans.push_back(static_cast< T1 >(i));
    }

    return ans;

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
 * @brief A single node for the model
 *
 * Each node contains all the information to compute the conditional probability
 * of the pruning algorithm at that node.
 *
 */
class Node {
public:
    unsigned int                 id;
    phylocounters::PhyloArray array;
    std::vector< unsigned int >              annotations;         ///< Observed annotations (only defined for APhyloModel)
    bool                                     duplication;

    std::vector< phylocounters::PhyloArray > arrays    = {};      ///< Arrays given all possible states
    Node *                                   parent    = nullptr; ///< Parent node
    std::vector< Node* >                     offspring = {};      ///< Offspring nodes
    std::vector< unsigned int >              idx_full  = {};
    bool                                     visited   = false;

    std::vector< double > probabilities; ///< The probability of observing each state

    Node() {};
    Node(unsigned int id_, bool duplication_) : id(id_), duplication(duplication_) {};
    Node(unsigned int id_, std::vector< unsigned int > annotations_, bool duplication_) :
        id(id_), annotations(annotations_), duplication(duplication_) {};
    ~Node() {};

    int get_parent() const {
        if (parent == nullptr)
            return -1;
        else
            return static_cast<int>(parent->id);
    };

    bool is_leaf() const {
        return offspring.size() == 0u;
    };

};

/**
 * @brief Annotated Phylo Model
 *
 */
class APhyloModel {
public:

    std::mt19937                       rengine;
    phylocounters::PhyloModel          model_full;
    unsigned int                       nfunctions;
    barry::Map< unsigned int, Node >   nodes;
    std::vector< unsigned int >        sequence;
    std::vector< bool >                visited;
    std::vector< std::vector< bool > > states;
    barry::MapVec_type< unsigned int > map_to_nodes;
    phylocounters::PhyloCounters       counters;
    bool                               initialized = false;

    APhyloModel();

    /**
     * @brief Construct a new APhyloModel object
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
     * @param counters An object of class `PhyloCounters`.
     */
    APhyloModel(
        std::vector< std::vector<unsigned int> > & annotations,
        std::vector< unsigned int > &              geneid,
        std::vector< int> &                        parent,
        std::vector< bool > &                      duplication
        );

    ~APhyloModel() {};

    void init();

    double operator()(std::vector< double > & par, unsigned int & i);
    void calc_sequence(Node * n = nullptr);

    double likelihood(const std::vector< double > & par);
    double likelihood_exact(const std::vector< double > & par);

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
    unsigned int nfuns() const;
    unsigned int nnodes() const;
    unsigned int nleafs() const;
    unsigned int nterms() const;
    ///@}

};

#endif
