#include "../../include/barry/barry.hpp"

#ifndef APHYLOMODEL_BONES_HPP
#define APHYLOMODEL_BONES_HPP 1

// The same need to be locked
RULE_FUNCTION(rule_blocked) {
    if (Array->get_cell(i, j) == 9u)
        return false;
    return true;
}

using namespace phylocounters;


template<typename T1, typename T2>
std::vector< T1 > caster(const std::vector< T2 > & vec) {

    std::vector< T1 > ans;
    ans.reserve(vec.size());

    for (auto &i : vec) {
        ans.push_back(*i);
    }

    return ans;

}

// Hasher
inline std::vector< double > keygen_const(const PhyloArray & array) {

    // Baseline data: nrows and columns
    std::vector< double > dat = {
        (double) array.nrow(), (double) array.ncol()
    };

    // State of the parent
    for (bool i : array.data->states) {
        dat.push_back(i ? 1.0 : 0.0);
    }

    // Free cells
    for (auto i = 0u; i < array.nrow(); ++i)
        for (auto j = 0u; j < array.ncol(); ++j)
            dat.push_back((double) array.get_cell(i, j, false));

    return dat;
}

// Hasher
inline std::vector< double > keygen_full(const PhyloArray & array) {

    // Baseline data: nrows and columns
    std::vector< double > dat = {
        (double) array.nrow(), (double) array.ncol()
    };

    // State of the parent
    for (bool i : array.data->states) {
        dat.push_back(i ? 1.0 : 0.0);
    }

    return dat;
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
    unsigned int id;
    PhyloArray array;

    std::vector< unsigned int > annotations;         ///< Observed annotations (only defined for APhyloModel)
    std::vector< PhyloArray >   arrays    = {};      ///< Arrays given all possible states
    Node *                      parent    = nullptr; ///< Parent node
    std::vector< Node* >        offspring = {};      ///< Offspring nodes
    std::vector< unsigned int > idx_cons  = {};      ///< Id of the constrained support.
    std::vector< unsigned int > idx_full  = {};
    bool                        visited   = false;

    std::vector< double > probabilities; ///< The probability of observing each state

    Node() {};
    Node(unsigned int id_) : id(id_) {};
    Node(unsigned int id_, std::vector< unsigned int > annotations_) :
        id(id_), annotations(annotations_) {};
    ~Node() {};

    int get_parent() const {
        if (parent == nullptr)
            return -1;
        else
            return (int) parent->id;
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
    PhyloModel                         model_const;
    PhyloModel                         model_full;
    unsigned int                       nfuns;
    barry::Map< unsigned int, Node >   nodes;
    std::vector< unsigned int >        sequence;
    std::vector< bool >                visited;
    std::vector< std::vector< bool > > states;
    barry::MapVec_type< unsigned int > map_to_nodes;
    PhyloCounters                      counters;

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
        std::vector< unsigned int > & geneid,
        std::vector< unsigned int> & parent
        );

    ~APhyloModel() {};

    void init();

    double operator()(std::vector< double > & par, unsigned int & i);
    void calc_sequence(Node * n = nullptr);
    double likelihood(const std::vector< double > & par);

    std::vector< double > get_probabilities() const;

    void set_seed(const unsigned int & s);
    std::vector< std::vector< bool > > simulate(
        const std::vector< double > & par
        );
};

#endif