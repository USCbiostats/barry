#ifndef GEESE_NODE_BONES
#define GEESE_NODE_BONES 1

/**
 * @brief A single node for the model
 *
 * Each node contains all the information to compute the conditional probability
 * of the pruning algorithm at that node.
 *
 */
class Node {
public:

    size_t id; ///< Id of the node (as specified in the input)
    size_t ord; ///< Order in which the node was created

    PhyloArray array;       ///< Array of the node
    std::vector< size_t >     annotations; ///< Observed annotations (only defined for Geese)
    bool                      duplication;

    std::vector< PhyloArray > arrays = {}; ///< Arrays given all possible states

    Node *                parent    = nullptr; ///< Parent node
    std::vector< Node* >  offspring = {};      ///< Offspring nodes
    std::vector< size_t > narray    = {};      ///< ID of the array in the model
    bool                  visited   = false;

    std::vector< double > subtree_prob; ///< Induced subtree probabilities
    std::vector< double > probability;  ///< The probability of observing each state
    
    /**
     * @name Construct a new Node object
     * 
     */
    ///@{
    
    Node() : ord(std::numeric_limits< size_t >::max()) {};
    Node(size_t id_, size_t ord_, bool duplication_);
    Node(size_t id_, size_t ord_, std::vector< size_t > annotations_, bool duplication_);
    
    // Move constructor
    Node(Node && x) noexcept;

    // Copy constructor
    Node(const Node & x);
    ///@}

    ~Node() {};

    int get_parent() const;

    size_t noffspring() const noexcept;
    bool is_leaf() const noexcept;

};

inline Node::Node(size_t id_, size_t ord_, bool duplication_)
    : id(id_), ord(ord_), duplication(duplication_) {

    return;
}

inline Node::Node(
    size_t id_,
    size_t ord_,
    std::vector< size_t > annotations_,
    bool duplication_
    ) : id(id_), ord(ord_), annotations(annotations_), duplication(duplication_) {}

inline Node::Node(Node && x) noexcept :
    id(x.id), ord(x.ord), array(std::move(x.array)),
    annotations(std::move(x.annotations)),
    duplication(x.duplication), arrays(std::move(x.arrays)),
    parent(std::move(x.parent)),
    offspring(std::move(x.offspring)),
    narray(std::move(x.narray)),
    visited(x.visited),
    subtree_prob(std::move(x.subtree_prob)),
    probability(std::move(x.probability)) {

    return;
    
}

inline Node::Node(const Node & x) :
    id(x.id), ord(x.ord), array(x.array), 
    annotations(x.annotations),
    duplication(x.duplication), arrays(x.arrays),
    parent(x.parent),
    offspring(x.offspring),
    narray(x.narray),
    visited(x.visited),
    subtree_prob(x.subtree_prob),
    probability(x.probability) {

        return;
    
}

inline int Node::get_parent() const {
    if (parent == nullptr)
            return -1;
    else
        return static_cast<int>(parent->id);
}
inline size_t Node::noffspring() const noexcept {

    return this->offspring.size();

}

inline bool Node::is_leaf() const noexcept {
    return offspring.size() == 0u;
}

#endif