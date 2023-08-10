
#ifndef GEESE_HPP
#define GEESE_HPP 1

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-types.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef GEESE_TYPES_HPP
#define GEESE_TYPES_HPP
/**
 * @name Convenient typedefs for Node objects.
 * */
/**
 * @brief Data definition for the `PhyloArray` class.
 * 
 * This holds basic information about a given node.
 * 
 */
class NodeData {
public:
  
    /**
     * Branch length.
     */
    std::vector< double > blengths = {};
    
    /**
     * State of the parent node.
     */
    std::vector< bool > states = {};
    
    /**
     * 
     */
    bool duplication = true;
    
    // NodeData() : blengths(0u), states(0u) {};
    
    NodeData(
        const std::vector< double > & blengths_,
        const std::vector< bool > & states_,
        bool duplication_ = true
    ) : blengths(blengths_), states(states_), duplication(duplication_) {};
    
    // ~NodeData() {};
  
};

class PhyloCounterData {
private:
    std::vector< size_t > data;
    std::vector< double > * counters;

public:
    PhyloCounterData(
        std::vector< size_t > data_,
        std::vector< double > * counters_ = nullptr
        ) : data(data_), counters(counters_) {};

    PhyloCounterData() : data(0u) {};

    size_t at(size_t d) {return data.at(d);};
    size_t operator()(size_t d) {return data.at(d);};
    size_t operator[](size_t d) {return data[d];};
    void reserve(size_t x) {return data.reserve(x);};
    void push_back(size_t x) {return data.push_back(x);};
    void shrink_to_fit()  {return data.shrink_to_fit();};
    size_t size() {return data.size();};

    std::vector< size_t >::iterator begin() {return data.begin();};
    std::vector< size_t >::iterator end() {return data.end();};

    bool empty() {return data.empty();};
    std::vector< double > * get_counters() {return counters;};

};

class PhyloRuleDynData {
public:
    const std::vector< double > * counts;
    size_t pos;
    size_t lb;
    size_t ub;
    size_t duplication;

    PhyloRuleDynData(
        const std::vector< double > * counts_,
        size_t pos_,
        size_t lb_,
        size_t ub_,
        size_t duplication_
        ) :
        counts(counts_), pos(pos_), lb(lb_), ub(ub_), duplication(duplication_) {};

    const double operator()() const
    {
        return (*counts)[pos];
    }
    
    ~PhyloRuleDynData() {};
    
};


typedef std::vector< std::pair< size_t, size_t > > PhyloRuleData;

///@{
typedef barry::BArrayDense<size_t, NodeData> PhyloArray;
typedef barry::Counter<PhyloArray, PhyloCounterData > PhyloCounter;
typedef barry::Counters< PhyloArray, PhyloCounterData> PhyloCounters;

typedef barry::Rule<PhyloArray,PhyloRuleData> PhyloRule;
typedef barry::Rules<PhyloArray,PhyloRuleData> PhyloRules;

typedef barry::Rule<PhyloArray,PhyloRuleDynData> PhyloRuleDyn;
typedef barry::Rules<PhyloArray,PhyloRuleDynData> PhyloRulesDyn;

typedef barry::Support<PhyloArray, PhyloCounterData, PhyloRuleData, PhyloRuleDynData > PhyloSupport;
typedef barry::StatsCounter<PhyloArray, PhyloCounterData> PhyloStatsCounter;
typedef barry::Model<PhyloArray, PhyloCounterData, PhyloRuleData, PhyloRuleDynData > PhyloModel;
typedef barry::PowerSet<PhyloArray, PhyloRuleData> PhyloPowerSet;
///@}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-types.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-node-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


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

    std::vector< bool > arrays_valid = {}; ///< Whether the arrays are valid according to the rules of the model.

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
    arrays_valid(std::move(x.arrays_valid)), 
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
    arrays_valid(x.arrays_valid),
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
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-node-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


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
         static_cast<double>(array.ncol())
    };

    // State of the parent
    dat.push_back(0.0);
    size_t count = 0u;
    for (bool i : array.D_ptr()->states) {
        dat[dat.size() - 1u] += (i ? 1.0 : 0.0) * pow(10, static_cast<double>(count));
        count++;
    }

    // Type of the parent
    dat.push_back(array.D_ptr()->duplication ? 1.0 : 0.0);

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
    barry::MapVec_type< size_t > map_to_nodes;
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
        bool use_reduced_sequence = true
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
    std::mt19937 *                     get_rengine();
    PhyloCounters *     get_counters();
    PhyloModel *        get_model();
    PhyloSupport *      get_support_fun();
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
    std::vector< size_t >        get_annotated_nodes() const; ///< Returns the ids of the nodes with at least one annotation

};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include "geese-bones.hpp"

#ifndef GEESE_MEAT_HPP
#define GEESE_MEAT_HPP 1

inline void Geese::init_node(Node & n)
{

    // Creating the phyloarray, nfunctions x noffspring
    n.array = PhyloArray(nfunctions, n.offspring.size());

    std::vector< bool > tmp_state = vector_caster<bool,size_t>(n.annotations);

    std::vector< double > blen(n.offspring.size(), 1.0);

    n.array.set_data(
        new NodeData(blen, tmp_state, n.duplication),
        true
    );

    // We initialize all with a zero since, if excluded from the pruning process,
    // We need to set it to one (as the result of the full integration).
    n.subtree_prob.resize(states.size(), 1.0);

    // Adding the data, first through functions
    for (size_t k = 0u; k < nfunctions; ++k)
    {

        // Then through the offspring
        size_t j = 0;
        for (auto& o : n.offspring)
        {

            // If leaf, then it may have an annotation
            if (o->is_leaf())
            {

                if (o->annotations[k] != 0)
                    n.array.insert_cell(k, j, o->annotations[k], false, false);

            }
            else
            {
                // [2022-02-11]: (IMPORTANT COMMENT!)
                // Otherwise, we fill it with a 0 so the support works correctly.
                // When adding an array from the interior, we don't need to deal
                // with the actual value as it is the powerset that matters. Using
                // nine instead will block the cell and stop the routine for computing
                // the values correctly
                n.array.insert_cell(k, j, 9u, false, false);

            }

            ++j;

        }

    }

    // We then need to set the powerset
    if (n.arrays.size() != states.size())
    {

        n.arrays.resize(states.size());
        n.narray.resize(states.size());
        n.arrays_valid.resize(states.size(), false);

    }
    
    // Here we have an issue: Some transitions may not be right
    // under the dynamic rules. So not all states can be valid.
    // The arrays and narrays need to be updated once the model
    // is initialized.
    //
    // The later is especially true for leaf nodes, where the
    // limitations are not known until the model is initialized.
    PhyloStatsCounter stats_counter;
    stats_counter.set_counters(model->get_counters());
    for (size_t s = 0u; s < states.size(); ++s)
    {

        n.arrays[s] = PhyloArray(n.array, false);

        n.arrays[s].set_data(
            new NodeData(blen, states[s], n.duplication),
            true
        );


        // Checking the rule. We need to make sure the counts match
        // the counts of the current array.
        if (model->get_rules_dyn() != nullptr)
        {
            // Once the array is ready, we can add it to the model
            stats_counter.reset_array(&n.arrays[s]);
            auto counts = stats_counter.count_all();

            PhyloRulesDyn dyn_rule(*model->get_rules_dyn());
            for (auto & r : dyn_rule)
                r.D().counts = &counts;

            // Finally, we can check if it can bee added. If not,
            // then we need to skip it.
            if (!dyn_rule(n.arrays[s], 0u, 0u))
                continue;
        }

        // Use try catch to run the following lines of code
        // only if the array is valid.
        try
        {
            n.narray[s] = model->add_array(n.arrays[s]);
        }
        catch (const std::exception & e)
        {
            auto err = std::string(e.what());

            err = "Array " + std::to_string(n.id) +
                " cannot be added to the model with error:\n" + err +
                "\n. This is likely due to a dynamic rule. " +
                "The array to be added was in the following state:";
                
            std::string state_str = "";
            for (const auto & ss : states[s])
                state_str += std::to_string(ss) + " ";

            err += state_str + "\n";

            throw std::runtime_error(err);
            
        }

        // n.narray[s] = model->add_array(n.arrays[s]);

        if (model->get_pset(n.narray[s])->size() != 0u)
            n.arrays_valid[s] = true;

    }

    return;

}

inline Geese::~Geese() {

    if (delete_support)
        delete model;

    if (delete_rengine)
        delete rengine;

    return;

}

inline void Geese::init(size_t bar_width) {

    // Initializing the model, if it is null
    if (this->model == nullptr)
    {

        this->model = new PhyloModel();

        this->delete_support = true;
        this->model->add_hasher(keygen_full);

        this->model->store_psets();

    }

    // Checking rseed, this is relevant when dealing with a flock. In the case of
    // flock, both model and rengine are shared.
    if (this->model->get_rengine() == nullptr) 
        this->model->set_rengine(this->rengine, false);

    // All combinations of the function
    PhyloPowerSet pset(nfunctions, 1u);

    pset.calc();

    states.reserve(pset.data.size());

    size_t i = 0u;

    for (auto& iter : pset.data)
    {

        states.push_back(std::vector< bool >(nfunctions, false));
        
        for (auto j = 0u; j < nfunctions; ++j)
        {

            if (!iter.is_empty(j, 0u, false))
                states[i][j] = true;

        }

        // Adding to map so we can look at it later on
        map_to_nodes.insert({iter.get_col_vec(0u, false), i});

        i++;

    }

    if (bar_width > 0u)
    {
        printf_barry("Initializing nodes in Geese (this could take a while)...\n");

        barry::Progress prog_bar(this->nnodes(), bar_width);

        // Iterating throught the nodes
        for (auto& iter : nodes)
        {

            // Only parents get a node
            if (!iter.second.is_leaf())
                this->init_node(iter.second); 
                
            prog_bar.next();
            
        }

        prog_bar.end();


    }
    else
    {

        // Iterating throught the nodes
        for (auto& iter : nodes)
        {

            // Only parents get a node
            if (!iter.second.is_leaf())
                this->init_node(iter.second); 
            
        }

    }

    // Resetting the sequence
    for (auto& n: this->nodes)
        n.second.visited = false;

    // The first time it is called, it need to generate the corresponding
    // hashes of the columns so it is fast to access then (saves time
    // hashing and looking in the map.)
    auto sup_arrays = model->get_pset_arrays();

    pset_loc.resize(sup_arrays->size());
    std::vector< size_t > tmpstate(nfunctions);

    for (auto s = 0u; s < sup_arrays->size(); ++s)
    {

        auto sup_array = sup_arrays->operator[](s);
        pset_loc[s].resize(sup_array.size());

        for (auto a = 0u; a < sup_array.size(); ++a)
        {

            for (auto o = 0u; o < sup_array[a].ncol(); ++o)
            {

                sup_array[a].get_col_vec(&tmpstate, o, false);
                pset_loc[s][a].push_back(map_to_nodes[tmpstate]);
                
            }   

        }

    }
    
    // So that others now know it was initialized
    initialized = true;

    return;

}

inline void Geese::inherit_support(const Geese & model_, bool delete_support_)
{
    
    if (this->model != nullptr)
        throw std::logic_error(
            "There is already a -model- in this Geese. Cannot set a -model- after one is present."
            );

    this->model = model_.model;

    this->delete_support = delete_support_;

    // And random number generation
    if (this->delete_rengine)
    {

        delete this->rengine;

        this->delete_rengine = false;

    }
    
    this->rengine = model_.rengine;
    
    return;

}

inline void Geese::update_annotations(
    size_t nodeid,
    std::vector< size_t > newann
) {

    // This can only be done if it has been initialized
    INITIALIZED()

    // Is this node present?
    if (nodes.find(nodeid) == nodes.end())
        throw std::length_error("The requested node is not present.");

    if (nodes[nodeid].annotations.size() != newann.size())
        throw std::length_error("Incorrect length of the new annotations.");

    // Resetting the annotations, and updating the stats from the
    // parent node
    nodes[nodeid].annotations = newann;

    // This only makes sense (for now) if it is a tip 
    if (!nodes[nodeid].is_leaf())
        return;

    init_node(*nodes[nodeid].parent);

    return;

}

inline void Geese::calc_sequence(Node * n)
{

    if (sequence.size() == nodes.size())
        return;

    // First iteration
    if (n == nullptr)
        n = &(nodes.begin()->second);

    // Here before?
    if (n->visited)
        return;

    n->visited = true;

    if (!n->is_leaf())
    {

        // iterating over its offspring, only if not there before
        for (auto& it : n->offspring)
        {

            if (!it->visited)
                calc_sequence(it);

        }

    }

    // Now, adding to the list and going to its parent
    sequence.push_back(n->id);

    if (n->parent == nullptr)
        return;

    // Go to the parent iff not visited
    if (!n->parent->visited)
        calc_sequence(n->parent);

    return;

}

inline void Geese::calc_reduced_sequence()
{

    // The criteria, if none of its decendants is annotated, then we can remove
    // the node from the model
    std::vector< bool > includeit(nodes.size(), false);

    for (auto& i : sequence)
    {

        Node & n = nodes[i];

        // We will count this at the end
        if (n.is_leaf())
        {

            for (size_t k = 0u; k < nfuns(); ++k)
                if (n.annotations[k] != 9u)
                {

                    includeit[n.ord] = true;
                    reduced_sequence.push_back(i);
                    break;

                }

        }
        else
        {

            // Checking, am I including any of my offspring?
            for (auto& o : n.offspring) 

                if (includeit[o->ord])
                {
                    
                    includeit[n.ord] = true;
                    reduced_sequence.push_back(i);
                    break;

                }

        }

    }

}

inline std::vector< double > Geese::get_probabilities() const
{

    std::vector< double > res;

    res.reserve(
        this->states.size() * nodes.size()
        );
    
    for (auto& i : sequence)
    {

        for (auto& p : this->nodes.at(i).subtree_prob)
            res.push_back(p);

    }

    return res;
    
}

inline size_t Geese::nfuns() const noexcept
{

    return this->nfunctions;

}

inline size_t Geese::nnodes() const noexcept
{

    return this->nodes.size();

}

inline size_t Geese::nleafs() const noexcept
{

    size_t n = 0u;

    for (auto& iter : this->nodes)
        if (iter.second.is_leaf())
            n++;

    return n;
}

inline size_t Geese::nterms() const
{

    INITIALIZED()
    return model->nterms() + this->nfuns();

}

inline size_t Geese::support_size() const noexcept
{

    if (model == nullptr)
        return 0u;

    return model->support_size();
    
}

inline std::vector< size_t > Geese::nannotations() const noexcept
{

    std::vector< size_t > ans = {this->n_zeros, this->n_ones};

    return ans;

}

inline std::vector< std::string > Geese::colnames() const
{

    return this->model->colnames();

}

inline size_t Geese::parse_polytomies(
    bool verb,
    std::vector< size_t > * dist
) const noexcept
{

    size_t largest = 0u;
    for (const auto& n : this->nodes)
    {

        if (n.second.is_leaf())
            continue;

        size_t noff = n.second.noffspring();

        if (dist)
            dist->push_back(noff);

        if (noff > 2u)
        {

            if (verb)
                printf_barry("Node id: %li has polytomy size %li\n", n.second.id, noff);
                
        }

        if (noff > largest)
            largest = noff;

    }

    return largest;

}

inline std::vector< std::vector<double> > Geese::observed_counts()
{

    // Making room for the output
    std::vector<std::vector<double>> ans;

    ans.reserve(nnodes());

    // Creating counter
    PhyloStatsCounter tmpcount;

    tmpcount.set_counters(this->model->get_counters());

    // Iterating through the nodes
    for (auto& n : nodes)
    {

        if (n.second.is_leaf())
        {

            ans.push_back({});
            continue;

        }

        PhyloArray tmparray(nfuns(), n.second.offspring.size());

        size_t j = 0u;

        for (auto& o : n.second.offspring)
        {

            for (size_t k = 0u; k < nfuns(); ++k)
            {

                if (o->annotations.at(k) != 0)
                {

                    tmparray.insert_cell(
                        k, j, o->annotations.at(k), false, false
                        );

                }

            }

            ++j;

        }

        std::vector< bool > tmp_state = vector_caster<bool,size_t>(
            n.second.annotations
            );

        std::vector< double > blen(n.second.offspring.size(), 1.0);

        tmparray.set_data(
            new NodeData(blen, tmp_state, n.second.duplication),
            true
        );

        tmpcount.reset_array(&tmparray);

        ans.push_back(tmpcount.count_all());

    }

    return ans;

}

inline void Geese::print_observed_counts()
{

    // Making room for the output
    std::vector<std::vector<double>> ans;
    ans.reserve(nnodes());

    // Creating counter
    PhyloStatsCounter tmpcount;
    tmpcount.set_counters(this->model->get_counters());

    // Iterating through the nodes
    for (auto& n : nodes) {

        if (n.second.is_leaf()) {
            ans.push_back({});
            continue;
        }

        PhyloArray tmparray(nfuns(), n.second.offspring.size());

        size_t j = 0u;
        for (auto& o : n.second.offspring) {
            for (size_t k = 0u; k < nfuns(); ++k) {
                if (o->annotations.at(k) != 0) {
                    tmparray.insert_cell(
                        k, j, o->annotations.at(k), false, false
                        );
                }
            }
            ++j;
        }

        std::vector< bool > tmp_state =vector_caster<bool,size_t>(n.second.annotations);
        std::vector< double > blen(n.second.offspring.size(), 1.0);
        tmparray.set_data(
            new NodeData(blen, tmp_state, n.second.duplication),
            true
        );

        tmpcount.reset_array(&tmparray);
        std::vector< double > counts = tmpcount.count_all();

        // Printing
        auto dpl = n.second.duplication ? "duplication" : "speciation";
        printf_barry("----------\n");
        printf_barry("nodeid: % 3li (%s)\nstate: [", n.second.id, dpl);
        for (size_t f = 0u; f < nfuns(); ++f)
            printf_barry("%i, ", (tmparray.D_ptr()->states[f] ? 1 : 0));

        printf_barry("]; Array:\n");
        tmparray.print();
        printf_barry("Counts: ");
        for (auto& c : counts)
            printf_barry("%.2f, ", c);
        printf_barry("\n");

    }

    return;

}

inline void Geese::print() const
{

    // Information about the tree:
    // - Number of functions
    // - Number of nodes and leafs
    // - Number of annotated leafs (0/1)
    printf_barry("GEESE\nINFO ABOUT PHYLOGENY\n");
    printf_barry("# of functions           : %li\n", this->nfuns());
    printf_barry("# of nodes [int; leaf]   : [%li; %li]\n", this->nnodes() - this->nleafs(), this->nleafs());
    printf_barry("# of ann. [zeros; ones]  : [%li; %li]\n", this->n_zeros, this->n_ones);
    printf_barry("# of events [dupl; spec] : [%li; %li]\n", this->n_dupl_events, this->n_spec_events);
    printf_barry("Largest polytomy         : %li\n", parse_polytomies(false));
    printf_barry("\nINFO ABOUT THE SUPPORT\n");
    this->model->print();

}

inline std::mt19937 * Geese::get_rengine()
{
    return this->rengine;
}

inline PhyloCounters * Geese::get_counters()
{
    return this->model->get_counters();
}

inline PhyloModel * Geese::get_model() {
    return this->model;
}

inline PhyloSupport * Geese::get_support_fun() {
    return this->model->get_support_fun();
}

inline std::vector< std::vector< bool > > Geese::get_states() const {
    return this->states;
}

inline std::vector< size_t > Geese::get_annotated_nodes() const {

    std::vector< size_t > ids(0u);
    for (auto & n : nodes)
    {

        // Counting non-9 annotations
        for (size_t f = 0u; f < nfuns(); ++f)
        {
            // If it has one non-9, then add it to the list
            // and continue to the next node.
            if (n.second.annotations[f] != 9u) {
                ids.push_back(n.second.id);
                break;
            }
        }

    }

    return ids;

}


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-meat-constructors.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include "geese-bones.hpp"

#ifndef GEESE_MEAT_CONSTRUCTORS_HPP
#define GEESE_MEAT_CONSTRUCTORS_HPP 1

inline Geese::Geese() {

    // In order to start...
    this->rengine         = new std::mt19937;
    this->delete_rengine  = true;
    this->model           = new PhyloModel();
    this->delete_support  = true;

    this->model->add_hasher(keygen_full);
    this->model->store_psets();

    return;
}

inline Geese::Geese(
    std::vector< std::vector<size_t> > & annotations,
    std::vector< size_t > &              geneid,
    std::vector< int > &                 parent,
    std::vector< bool > &                duplication
) {

    // In order to start...
    this->rengine         = new std::mt19937;
    this->delete_rengine  = true;
    this->model           = new PhyloModel();
    this->delete_support  = true;

    this->model->add_hasher(keygen_full);
    this->model->store_psets();

    // Check the lengths
    if (annotations.size() == 0u)
        throw std::logic_error("Annotations is empty");

    nfunctions = annotations.at(0u).size();

    // size_t n = annotations.size();
    for (auto& iter : annotations)
    {

        if (iter.size() != nfunctions)
            throw std::length_error(
                "Not all the annotations have the same length"
                );

    }

    // Grouping up the data by parents -----------------------------------------
    for (size_t i = 0u; i < geneid.size(); ++i)
    {

        // Temp vector with the annotations
        std::vector< size_t > & funs(annotations.at(i));

        // Case 1: Not the root node, and the parent does not exists
        if ((parent.at(i) >= 0) && (nodes.find(parent.at(i)) == nodes.end()))
        {

            // Adding parent
            auto key_par = nodes.insert({
                parent.at(i),
                Node(parent.at(i), std::numeric_limits< size_t >::max(), true)
            });

            // Case 1a: i does not exists
            if (nodes.find(geneid.at(i)) == nodes.end())
            {

                auto key_off = nodes.insert({
                    geneid.at(i),
                    Node(geneid.at(i), i, funs, duplication.at(i))
                    });

                // Adding the offspring to the parent
                key_par.first->second.offspring.push_back(
                    &key_off.first->second
                );

                // Adding the parent to the offspring
                key_off.first->second.parent = &key_par.first->second;

            } else { // Case 1b: i does exists (we saw it earlier)

                // We just need to make sure that we update it!
                nodes[geneid.at(i)].duplication = duplication.at(i);
                nodes[geneid.at(i)].annotations = funs;
                nodes[geneid.at(i)].parent      = &nodes[parent.at(i)];
                nodes[geneid.at(i)].ord         = i;

                nodes[parent.at(i)].offspring.push_back(
                    &nodes[geneid.at(i)]
                );

            }

        } else { // Case 2: Either this is the root, or the parent does exists

            // Case 2a: i does not exists (but its parent does)
            if (nodes.find(geneid.at(i)) == nodes.end())
            {

                // Adding i
                auto key_off = nodes.insert({
                    geneid.at(i),
                    Node(geneid.at(i), i, funs, duplication.at(i))
                    });

                // We only do this if this is not the root
                if (parent.at(i) >= 0)
                {

                    nodes[parent.at(i)].offspring.push_back(
                        &key_off.first->second
                    );

                    // Adding the parent to the offspring
                    key_off.first->second.parent = &nodes[parent.at(i)];

                }

            } else { // Case 2b: i does exists (and so does its parent)

                // We just need to make sure that we update it!
                nodes[geneid.at(i)].duplication = duplication.at(i);
                nodes[geneid.at(i)].annotations = funs;
                nodes[geneid.at(i)].ord         = i;

                if (parent.at(i) >= 0)
                {

                    nodes[geneid.at(i)].parent = &nodes[parent.at(i)];
                    nodes[parent.at(i)].offspring.push_back(
                        &nodes[geneid.at(i)]
                    );

                }

            }
        }

    }

    // Verifying that all have the variable ord
    for (auto& n : nodes)
    {

        Node & node = n.second;

        // Checking variable
        if (node.ord == std::numeric_limits< size_t >::max())
        {

            const char *fmt = "Node id %i was not included in geneid.";
            int sz = std::snprintf(nullptr, 0, fmt, node.id);
            std::vector<char> buf(sz + 1);
            std::snprintf(&buf[0], buf.size(), fmt, node.id);
            throw std::logic_error(&buf[0]);

        }

        // Checking duplication
        if (node.duplication != duplication[node.ord])
        {

            const char *fmt = "Node id %i's duplication was not properly recorded.";
            int sz = std::snprintf(nullptr, 0, fmt, node.id);
            std::vector<char> buf(sz + 1);
            std::snprintf(&buf[0], buf.size(), fmt, node.id);
            throw std::logic_error(&buf[0]);

        }

        // Counting the type of annotations
        if (node.is_leaf())
        {

            for (const auto & a : node.annotations)
            {

                if (a == 1u)
                    this->n_ones++;
                else if (a == 0u)
                    this->n_zeros++;

            }

        } else {

            if (node.duplication)
                this->n_dupl_events++;
            else
                this->n_spec_events++;

        }

    }


    // Computing the pruning sequence.
    calc_sequence();
    calc_reduced_sequence();

    // Are the sequences OK?
    if (this->sequence.size() != this->nnodes())
        throw std::logic_error("The pruning sequence's length is different from nnodes(). This should not happen! (contact the developers).");

    return;

}

inline Geese::Geese(const Geese & model_, bool copy_data) : 
    states(model_.states),
    n_zeros(model_.n_zeros),
    n_ones(model_.n_ones),
    n_dupl_events(model_.n_dupl_events),
    n_spec_events(model_.n_spec_events),
    nfunctions(model_.nfunctions),
    nodes(model_.nodes),
    map_to_nodes(model_.map_to_nodes),
    pset_loc(model_.pset_loc),
    sequence(model_.sequence),
    reduced_sequence(model_.reduced_sequence),
    initialized(model_.initialized) {

    
    // Replicating -------------------------------------------------------------
    if (copy_data)
    {

        if (model_.rengine != nullptr)
        {
            rengine = new std::mt19937(*(model_.rengine));
            delete_rengine = true;
        }

        if (model_.model != nullptr)
        {
            model = new PhyloModel(*(model_.model));
            delete_support = true;
        }

    } else {
        
        if (model_.rengine != nullptr)
        {
            rengine = model_.rengine;
            delete_rengine = false;
        }

        if (model_.model != nullptr)
        {
            model = model_.model;
            delete_support = false;
        }

    }

    // These should not be necesary as they are already initialized.
    // this->model->set_keygen(keygen_full);
    // this->model->store_psets();

    // Dealing with the nodes is a bit different -------------------------------
    auto revseq = this->sequence;
    std::reverse(revseq.begin(), revseq.end());

    for (auto& i : revseq)
    {

        // Leaf do not have offspring
        if (this->nodes[i].is_leaf())
            continue;

        // Clearing offspring
        this->nodes[i].offspring.clear();

        // I cannot directly access the node since, if non existent, it will 
        // create an entry with it (alegedly).
        auto n = model_.nodes.find(i);

        for (const auto& off : n->second.offspring)
            this->nodes[i].offspring.push_back(&this->nodes[off->id]);

    }

    return;
  
}

// Constructor move
inline Geese::Geese(Geese && x) noexcept :
    rengine(nullptr),
    model(nullptr),
    states(std::move(x.states)),
    n_zeros(std::move(x.n_zeros)),
    n_ones(std::move(x.n_ones)),
    n_dupl_events(std::move(x.n_dupl_events)),
    n_spec_events(std::move(x.n_spec_events)),
    nfunctions(x.nfunctions),
    nodes(std::move(x.nodes)),
    map_to_nodes(std::move(x.map_to_nodes)),
    pset_loc(std::move(x.pset_loc)),
    sequence(std::move(x.sequence)),
    reduced_sequence(std::move(x.reduced_sequence)),
    initialized(x.initialized)
{

    if (x.delete_rengine)
    {

        rengine = new std::mt19937(*x.rengine);
        delete_rengine = true;

    } else {

        rengine = x.rengine;
        delete_rengine = false;

    }

    if (x.delete_support)
    {

        model = new PhyloModel(*x.model);
        delete_support = true;

    } else {

        model = x.model;
        delete_support = false;
        
    }

    // Figuring out if model needs to be updated
    if ((model != nullptr) && (x.delete_support | x.delete_rengine))
        model->set_rengine(rengine, false);

    return;

}



#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-meat-constructors.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-meat-likelihood.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef GEESE_MEAT_LIKELIHOOD_HPP
#define GEESE_MEAT_LIKELIHOOD_HPP 1

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models//geese/geese-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


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
         static_cast<double>(array.ncol())
    };

    // State of the parent
    dat.push_back(0.0);
    size_t count = 0u;
    for (bool i : array.D_ptr()->states) {
        dat[dat.size() - 1u] += (i ? 1.0 : 0.0) * pow(10, static_cast<double>(count));
        count++;
    }

    // Type of the parent
    dat.push_back(array.D_ptr()->duplication ? 1.0 : 0.0);

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
    barry::MapVec_type< size_t > map_to_nodes;
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
        bool use_reduced_sequence = true
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
    std::mt19937 *                     get_rengine();
    PhyloCounters *     get_counters();
    PhyloModel *        get_model();
    PhyloSupport *      get_support_fun();
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
    std::vector< size_t >        get_annotated_nodes() const; ///< Returns the ids of the nodes with at least one annotation

};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models//geese/geese-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



inline double Geese::likelihood(
    const std::vector< double > & par,
    bool as_log,
    bool use_reduced_sequence
) {

    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    double ll = 0.0;

    Node * n_off;

    // Following the prunning sequence
    std::vector< size_t > * preseq;

    if (use_reduced_sequence)
    {

        preseq = &this->reduced_sequence;

    }
    else
    {   

        preseq = &this->sequence;

    }

    // The first time it is called, it need to generate the corresponding
    // hashes of the columns so it is fast to access then (saves time
    // hashing and looking in the map.)
    auto arrays2support = model->get_arrays2support();

    for (auto& i : *preseq)
    {

        // We cannot compute probability at the leaf, we need to continue
        if (this->nodes[i].is_leaf())
            continue;

        // Since we are using this a lot...
        Node & node = nodes[i];

        // Iterating through states
        for (size_t s = 0u; s < states.size(); ++s)
        {

            // Starting the prob
            double totprob = 0.0;

            // If the transition doesn't make sense, then we skip it.
            // This is determined during the construction of the node, when
            // the rule_dyn is called.
            if (!node.arrays_valid[s])
            {
                node.subtree_prob[s] = 0.0;
                continue;
            }

            // Retrieving the sets of arrays
            const std::vector< PhyloArray > * psets =
                model->get_pset(node.narray[s]);

            const std::vector<double> * psets_stats =
                model->get_pset_stats(node.narray[s]);

            std::vector< std::vector< size_t > > & locations = pset_loc[
                arrays2support->operator[](node.narray[s])
                ];
            
            // Summation over all possible values of X
            size_t nstate = 0u;
            size_t narray = 0u;
            for (auto x = psets->begin(); x != psets->end(); ++x)
            {

                if (!x->is_dense())
                    throw std::logic_error("This is only supported for dense arrays.");

                std::vector< size_t > & location_x = locations[narray++];

                // Extracting the possible values of each offspring
                double off_mult = 1.0;

                for (auto o = 0u; o < x->ncol(); ++o)
                {

                    // Setting the node
                    n_off = node.offspring[o];
                    
                    // In the case that the offspring is a leaf, then we need to
                    // check whether the state makes sense.
                    if (n_off->is_leaf())
                    {
                        for (auto f = 0u; f < nfunctions; ++f)
                        {
                            if (n_off->annotations[f] != 9u)
                            {

                                if (x->operator()(f, o) != n_off->annotations[f])
                                {

                                    off_mult = -1.0;
                                    break;

                                }
                                
                            }

                        }

                        // Going out
                        if (off_mult < 0)
                            break;
                
                        continue;

                    }

                    // Retrieving the location to the respective set of probabilities
                    off_mult *= node.offspring[o]->subtree_prob[location_x[o]];

                }

                // Is this state valid?
                if (off_mult < 0.0)
                {

                    ++nstate;
                    continue;
                    
                }

                // Multiplying by P(x|x_n), the transition probability
                std::vector< double > temp_stats(par0.size(), 0.0);
                for (auto p = 0u; p < par0.size(); ++p)
                    temp_stats[p] = psets_stats->operator[](par0.size() * nstate + p);

                nstate++;

                // Use try catch in the following line
                try {
                    off_mult *= model->likelihood(
                        par0,
                        temp_stats,
                        node.narray[s]
                    );
                } catch (std::exception & e) {

                    auto err = std::string(e.what());

                    std::string state_str = "";
                    for (const auto & ss : states[s])
                        state_str += std::to_string(ss) + " ";

                    err = "Error computing the likelihood at node " +
                        std::to_string(node.id) + " with state " + state_str +
                        ". Error message:\n" +
                        err;

                    throw std::runtime_error(err);
                    
                }

                // Adding to the total probabilities
                totprob += off_mult;

            }

            // Setting the probability at the node
            node.subtree_prob[s] = totprob;

        }

        // All probabilities should be completed at this point
        if (node.parent == nullptr)
        {

            for (size_t s = 0u; s < states.size(); ++s)
            {

                double tmpll = 1.0;

                for (auto k = 0u; k < nfunctions; ++k)
                {

                    tmpll *= states[s][k] ? par_root[k] : (1 - par_root[k]);

                }

                ll += tmpll * node.subtree_prob[s];

            }
        }

    }

    // In the case that the sequence is empty, then it means
    // that we are looking at a completely unnanotated tree,
    // thus the likelihood should be one
    if (preseq->size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 1.0;


    return as_log ? std::log(ll) : ll;

}
#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-meat-likelihood.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-meat-likelihood_exhaust.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



#ifndef GEESE_MEAT_LIKELIHOOD_EXHAUST_HPP
#define GEESE_MEAT_LIKELIHOOD_EXHAUST_HPP 1
// #include "../../barry.hpp"
// #include "geese-bones.hpp" 

inline double Geese::likelihood_exhaust(
    const std::vector< double > & par
)
{

    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    // This is only worthwhile if the number of nodes is small
    if (this->nnodes() > 6)
        throw std::overflow_error("Too many nodes! Exhaust calculation of likelihood cannot be done for such cases.");

    if (this->nfuns() > 3)
        throw std::overflow_error("Too many functions! Exhaust calculation of likelihood cannot be done for such cases.");

    // Computing all combinations ----------------------------------------------
    PhyloArray base(nfuns(), nnodes());
    for (auto& n : nodes)
    {

        for (size_t i = 0u; i < nfuns(); ++i)
            base(i, n.second.ord) = n.second.annotations[i];
            
    }

    PhyloPowerSet pset(base);//this->nfuns(), this->nnodes());
    pset.add_rule(
            rule_empty_free<PhyloArray,PhyloRuleData>,
            PhyloRuleData()
            );
    pset.calc();

    // Inverse sequence
    std::vector< size_t > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    double totprob = 0.0;
    
    // This vector says whether the probability has to be included in 
    // the final likelihood or not.
    for (size_t p = 0u; p < pset.size(); ++p)
    {
        
        // ith state
        const PhyloArray * s = &pset[p];
        
        // Following the sequence
        double prob = 1.0;
        std::vector< size_t > tmpstates(this->nfuns());

        Node * node;
        for (auto& i : preorder)
        {

            node = &nodes[i];
            std::fill(tmpstates.begin(), tmpstates.end(), 0u);
            s->get_col_vec(&tmpstates, node->ord, false);

            // Root node first
            if (node->parent == nullptr)
            {               
                // Since it is the root, the first probability is computed using
                // the root only
                for (auto k = 0u; k < this->nfuns(); ++k)
                    prob *= tmpstates[k] == 1u ? par_root[k] : (1.0 - par_root[k]);

            }
            else if (node->is_leaf())
                continue;

            // Computing the transition
            PhyloArray transition(nfuns(), node->offspring.size());

            std::vector< double > bl(node->offspring.size(), 1.0);

            std::vector< bool > sl = vector_caster<bool,size_t>(tmpstates);

            transition.set_data(
                new NodeData(bl, sl, node->duplication),
                true
            );

            // Filling the array
            for (size_t a = 0u; a < nfuns(); ++a)
            {

                for (size_t o = 0u; o < node->offspring.size(); ++o)
                {

                    if (s->get_cell(a, node->offspring[o]->id) == 1u)
                        transition(a, o) = 1u;

                }

            }

            prob *= this->model->likelihood(
                par0,
                transition,
                node->narray[this->map_to_nodes[tmpstates]],
                false
                );

        }

        totprob += prob;
    }

    return totprob;

}
#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-meat-likelihood_exhaust.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-meat-simulate.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef GEESE_MEAT_SIMULATE_HPP
#define GEESE_MEAT_SIMULATE_HPP 1

inline void Geese::set_seed(const size_t & s) {
    rengine->seed(s);
}

inline std::vector< std::vector< size_t > > Geese::simulate(
    const std::vector< double > & par
    ) {

    INITIALIZED()

    // Splitting the probabilities
    std::vector< double > par0(par.begin(), par.end() - nfunctions);
    std::vector< double > par_root(par.end() - nfunctions, par.end());

    // Scaling root
    for (auto& p : par_root) {
        p = std::exp(p)/(std::exp(p) + 1);
    }

    // Making room 
    std::vector< std::vector< size_t > > res(nodes.size());

    // Inverse sequence
    std::vector< size_t > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    // Generating probabilities at the root-level (root state)
    std::vector< double > rootp(states.size(), 1.0);
    const Node & rootnode = nodes[preorder[0u]];
    for (size_t i = 0u; i < rootp.size(); ++i)
    {

        if (!rootnode.arrays_valid[i])
        {
            rootp[i] = 0.0;
            continue;
        }

        for (size_t j = 0u; j < nfuns(); ++j)
            rootp[i] *= states[i][j] ? par_root[j] : (1.0 - par_root[j]);
    }

    // Preparing the random number generator
    std::uniform_real_distribution<> urand(0, 1);
    double r         = urand(*rengine);
    size_t idx = 0u;
    double cumprob = rootp[idx];
    while ((idx < rootp.size()) && (cumprob < r))
    {
        cumprob += rootp[++idx];
    }

    #ifdef BARRY_DEBUG
    
    // auto totprob = std::accumulate(rootp.begin(), rootp.end(), 0.0);
    // if (totprob < 0.9999999999999999 || totprob > 1.0000000000000001)
    //     throw std::runtime_error("Root probabilities do not sum to 1!"
    //         " (totprob = " + std::to_string(totprob) + ")");
    
    #endif
    
    // We now know the state of the root
    res[nodes[preorder[0u]].ord] =
        vector_caster< size_t, bool>(states[idx]);

    // Going in the opposite direction
    for (auto& i : preorder)
    {

        if (nodes[i].is_leaf())
            continue;

        const Node & n = nodes[i];

        // Getting the state of the node      
        size_t l = map_to_nodes[res[n.ord]];

        // Given the state of the current node, sample the state of the
        // offspring, all based on the current state
        // auto z = n.narray;
        auto tmp = model->sample(n.narray[l], par0);

        // Iterating through the offspring to assign the state
        for (size_t j = 0u; j < n.offspring.size(); ++j)
            res[n.offspring[j]->ord] = tmp.get_col_vec(j, false);

    }

    return res;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-meat-simulate.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-meat-predict.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include "geese-bones.hpp"

#ifndef GEESE_MEAT_PREDICT_HPP
#define GEESE_MEAT_PREDICT_HPP 1

inline std::vector< std::vector<double> > Geese::predict_backend(
    const std::vector< double > & par,
    bool use_reduced_sequence,
    const std::vector< size_t > & preorder
)
{

    // Splitting the probabilities
    std::vector< double > par_terms(par.begin(), par.end() - nfuns());
    std::vector< double > par_root(par.end() - nfuns(), par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    // Generating probabilities at the root-level (root state)
    std::vector< double > rootp(this->states.size(), 1.0);
    for (size_t s = 0u; s < rootp.size(); ++s)
    {

        for (size_t f = 0u; f < nfuns(); ++f)
            rootp[s] *= states[s][f] ? par_root[f] : (1.0 - par_root[f]);
        
    }

    // Making room 
    std::vector< std::vector<double> > res(
        nnodes(), std::vector<double>(nfuns())
        );

    // Step 1: Computing the probability at the root node
    std::vector< double > tmp_prob(nfuns(), 0.0);
    size_t root_id = preorder[0u];
    Node * tmp_node      = &nodes[root_id];
    tmp_node->probability.resize(states.size(), 0.0);
    double tmp_likelihood = likelihood(par, false, use_reduced_sequence);

    for (size_t s = 0u; s < states.size(); ++s)
    {

        // This is mostly relevant when the root is connected
        // to a leaf node, where the array may not be valid 
        // according to the dynamic rules, i.e., gains/loses.
        // If 2 genes have a function, the root state is zero,
        // and the rule is no more than one gain, then the
        // state is not valid.
        if (!tmp_node->arrays_valid[s])
        {
            tmp_node->probability[s] = 0.0;
            continue;
        }

        // Overall state probability P(x_s | D)
        tmp_node->probability[s] = tmp_node->subtree_prob[s] * rootp[s] /
            tmp_likelihood;

        // Marginalizing the probabilities P(x_sf | D)
        for (size_t f = 0u; f < nfuns(); ++f)
        {

            // Since the probability, the expected value, is for
            // observing an x = 1, then we need to make sure that we
            // are multiplying by the corresponding state
            if (states[s][f])
                tmp_prob[f] += tmp_node->probability[s];

        }
        

    }

    // Storing the final prob
    res[nodes[preorder[0u]].ord] = tmp_prob;
    size_t n_pars = par_terms.size();

    for (auto & i : preorder)
    {

        // Leafs have nothing to do here
        Node & parent = nodes[i];
        if (parent.is_leaf())
            continue;

        // Creating space.
        std::vector< std::vector< double > > everything_below(states.size());
        std::vector< std::vector< double > > everything_above(states.size());
        std::vector< std::vector< PhyloArray > > psets(states.size());

        // Making space for the offspring
        for (auto & off : parent.offspring)
        {
            off->probability.resize(states.size(), 0.0);
            std::fill(off->probability.begin(), off->probability.end(), 0.0);
        }

        // Iterating through the parent states
        for (size_t s = 0u; s < states.size(); ++s)
        {

            // Only if it is valid: See the same block applied to the root
            if (!parent.arrays_valid[s])
            {
                parent.probability[s] = 0.0;
                continue;
            }

            // Retrieving powerset of stats and arrays
            // it is not const since we will flip the states back and forth
            // to generate the key
            const auto & pset_arrays = model->get_pset(parent.narray[s]);
            const std::vector<double> * pset_target = model->get_pset_stats(parent.narray[s]);

            for (size_t p = 0u; p < pset_arrays->size(); ++p)
            {

                // Corresponding graph and target stats
                const PhyloArray & array_p = pset_arrays->at(p);
                std::vector<double> target_p(n_pars, 0.0);
                for (size_t par_i = 0u; par_i < target_p.size(); ++par_i)
                    target_p[par_i] = pset_target->operator[](p * n_pars + par_i);

                PhyloArray tmp_array(nfuns(), array_p.ncol());
                tmp_array += array_p;

                // Adding to the map, we only do this during the first run,
                // afterwards, we need to actually look for the array.
                bool in_the_set = true; /// < True if the array belongs to the set
                
                // Everything below just need to be computed only once
                // and thus, if already added, no need to go through all of this!
                double everything_below_p = 1.0;
                for (size_t off = 0u; off < parent.offspring.size(); ++off)
                {

                    // Below leafs, the everything below is 1.
                    if (parent.offspring[off]->is_leaf())
                    {

                        // But we can only includ it if the current state actually
                        // matches the leaf data (otherwise the prob is 0)
                        const auto & off_ann = parent.offspring[off]->annotations;
                        for (size_t f = 0u; f < nfuns(); ++f)
                        {

                            if ((off_ann[f] != 9u) && (off_ann[f] != array_p(f, off)))
                            {
                                in_the_set = false;
                                break;
                            }
                                
                        }

                        if (!in_the_set)
                            break;

                        continue;

                    } else {

                        // Getting the offspring state, and how it maps, only
                        // if it is not an offspring
                        const auto & off_state = array_p.get_col_vec(off);
                        size_t loc = this->map_to_nodes[off_state];

                        everything_below_p *= parent.offspring[off]->subtree_prob[loc];

                    }

                }

                // If it is not in the set, then continue to the next array
                if (!in_the_set)
                    continue;

                psets[s].push_back(array_p); // Generating a copy
                everything_below[s].push_back(everything_below_p);

                // The first run, we only need to grow the list
                everything_above[s].push_back(
                    model->likelihood(
                        par_terms, target_p, parent.narray[s], false
                    ) *  parent.probability[s] / parent.subtree_prob[s]
                );


            } // end for psets
            
        } // end for states

        // Marginalizing at the state level
        for (size_t s = 0u; s < states.size(); ++s)
        {

            // Only if it is valid
            if (!parent.arrays_valid[s])
                continue;

            for (size_t p = 0u; p < everything_above[s].size(); ++p)
            {

                // p-th pset
                const auto & pset_p = psets[s][p];

                // Updating the probability (it is the product)
                everything_above[s][p] *= everything_below[s][p];

                for (size_t off = 0u; off < parent.offspring.size(); ++off)
                {

                    // Figuring out the state of the offspring
                    size_t off_s = this->map_to_nodes[pset_p.get_col_vec(off)];
                    parent.offspring[off]->probability[off_s] += everything_above[s][p];


                }

            }
        }

        // Finally, we can marginalize the values at the 
        // gene function level.
        for (const auto & off : parent.offspring)
        {
            for (size_t s = 0u; s < states.size(); ++s)
            {

                // Only if it is valid
                if (!parent.arrays_valid[s])
                    continue;

                for (size_t f = 0u; f < nfuns(); ++f)
                    if (states[s][f]) 
                        res[off->ord][f] += off->probability[s];

            }

            // Checking that probabilities add up to one
            for (size_t f = 0u; f < nfuns(); ++f)
            {
                if ((res[off->ord][f] > 1.00001) || (res[off->ord][f] < -.0000009))
                {
                    auto msg = "[geese] Out-of-range probability for node.id " +
                        std::to_string(off->id) + " for function " +
                        std::to_string(f) + ": " +
                        std::to_string(res[off->ord][f]);

                    throw std::logic_error(msg);
                    
                } 

                if (res[off->ord][f] > 1.0)
                    res[off->ord][f] = 1.0;
                else if (res[off->ord][f] < 0.0)
                    res[off->ord][f] = 0.0;

            }
   

        }

    } // end for over preorder
        
    return res;

}

inline std::vector< std::vector<double> > Geese::predict(
    const std::vector< double > & par,
    std::vector< std::vector< double > > * res_prob,
    bool leave_one_out,
    bool only_annotated,
    bool use_reduced_sequence
)
{

    INITIALIZED()

    // Inverse sequence
    std::vector< size_t > preorder;
    if (only_annotated)
        preorder = this->reduced_sequence;
    else
        preorder = this->sequence;

    std::reverse(preorder.begin(), preorder.end());

    // Full prediction (first run, right now I am doing this
    // twice. Need to fix in the future)
    std::vector< std::vector<double> > res = predict_backend(
        par, use_reduced_sequence, preorder
        );

    // If the user requires the probability matrix per state
    if (res_prob != nullptr)
    {

        res_prob->resize(nnodes());
        for (auto& i : sequence)
            res_prob->at(nodes[i].ord) = nodes[i].probability;

    }


    // In this case, we need to update the predictions, mostly of the annotated
    // leaf nodes. Because of 
    if (leave_one_out)
    {

        std::vector< size_t > default_empty(nfuns(), 9u);
        for (auto& n : nodes)
        {

            if (n.second.is_leaf()) {

                Node & ntmp = n.second;

                // We only make the changes if it is not all missing
                bool use_it = false;
                for (auto& n_state : ntmp.annotations)
                    if (n_state != 9u)
                    {

                        use_it = true;
                        break;

                    }
                

                if (!use_it)
                    continue;

                // Recording the original annotation
                auto old_ann = ntmp.annotations;

                // Removing the entire gene
                update_annotations(ntmp.id, default_empty);

                // Making the prediction
                res[ntmp.ord] = (
                    predict_backend(par, use_reduced_sequence, preorder)
                )[ntmp.ord];

                // Restoring the gene
                update_annotations(ntmp.id, old_ann);

                if (res_prob != nullptr)
                    res_prob->operator[](ntmp.ord) = ntmp.probability;


            }

        }

    }

    
    return res;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-meat-predict.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-meat-predict_exhaust.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



#ifndef GEESE_MEAT_PREDICT_EXHAUST_HPP
#define GEESE_MEAT_PREDICT_EXHAUST_HPP 1

inline std::vector< std::vector<double> > Geese::predict_exhaust(
    const std::vector< double > & par
) {

    INITIALIZED()

    // This is only worthwhile if the number of nodes is small
    if (this->nnodes() > 6)
        throw std::overflow_error("Too many nodes! Exhaust calculation of likelihood cannot be done for such cases.");

    if (this->nfuns() > 2)
        throw std::overflow_error("Too many functions! Exhaust calculation of prediction cannot be done for such cases.");


    // Generating the sequence preorder sequence -------------------------------
    std::vector< size_t > preorder(this->sequence);
    std::reverse(preorder.begin(), preorder.end());

    std::vector< std::vector< double > > res = predict_exhaust_backend(
        par, preorder
        );

    // Looping to do LOO
    std::vector< size_t > annotated_ids = this->get_annotated_nodes();
    std::vector< size_t > missing_vec(nfuns(), 9u);
    for (auto & i : annotated_ids) {

        Node & n = nodes[i];

        auto old_ann = n.annotations;
        update_annotations(i, missing_vec);

        res[n.ord] = predict_exhaust_backend(par, preorder)[n.ord];

        update_annotations(i, old_ann);

    }

    return res;

}

inline std::vector< std::vector<double> > Geese::predict_exhaust_backend(

    const std::vector< double > & par,
    const std::vector< size_t > & preorder
) {

    // Processing the probabilities --------------------------------------------
    std::vector< double > par_terms(par.begin(), par.end() - nfuns());
    std::vector< double > par_root(par.end() - nfuns(), par.end());

    // Scaling root
    for (auto& p : par_root)
        p = std::exp(p)/(std::exp(p) + 1);

    double baseline_likelihood = this->likelihood(par);

    // Computing all combinations ----------------------------------------------
    // The base PhyloArray will store the original set of annotations.
    PhyloArray base(nfuns(), nnodes());
    for (auto& n : nodes)
    {

        for (size_t f = 0u; f < nfuns(); ++f)
            base(f, n.second.ord) = n.second.annotations[f];

    }

    PhyloPowerSet pset(base);//this->nfuns(), this->nnodes());
    pset.add_rule(
            rule_empty_free<PhyloArray,PhyloRuleData>,
            PhyloRuleData()
            );
    pset.calc();
    
    // Making space for the expected values
    std::vector< double > expected(nnodes() * nfuns(), 0.0);
    
    // This vector says whether the probability has to be included in 
    // the final likelihood or not.
    for (size_t p = 0u; p < pset.size(); ++p)
    {
        
        // ith state
        const PhyloArray * s = &pset[p];
        
        // Computing the likelihood of the state s        
        double current_prob = 1.0;
        for (auto & o: preorder)
        {
            // Getting the corresponding node
            Node & n = nodes[o];

            // Nothing to do at the leaf level (leafs are calculated from parents)
            if (n.is_leaf())
                continue;

            // Extracting the parent column (without checking boundaries)
            auto par_state = s->get_col_vec(n.ord, false);

            // Need to compute the root probability (if we havent before)
            if (n.parent == nullptr)
            {

                for (size_t f = 0u; f < nfuns(); ++f)
                    current_prob *= par_state[f] ? par_root[f] : (1.0 - par_root[f]);

            }
        
            // Generating a copy of the observed array
            // (data is copied so that we can chage the state of the parent)
            PhyloArray tmparray(n.array, true);

            // Updating the state of the parent
            for (size_t f = 0u; f < nfuns(); ++f)
                tmparray.D_ptr()->states[f] = par_state[f] == 1u;

            // Updating offspring annotations
            int loc = 0;
            for (auto & off : n.offspring) {
                
                for (size_t f = 0u; f < nfuns(); ++f)
                {

                    if (s->operator()(f, off->ord) == 1u)
                        tmparray(f, loc) = 1u;
                    else
                        tmparray.rm_cell(f, loc);

                }

                // Next offspring start in the next column of the array, Duh.
                ++loc;

            }
            
            // Computing the likelihood
            current_prob *= model->likelihood(par_terms, tmparray, -1, false);

        }
            // this->update_annotations(n.second.id, s->get_col_vec(n.second.ord));
        
        // Adding to the overall probability
        for (auto & n: nodes)
            for (size_t j = 0u; j < nfuns(); ++j)
                expected[n.second.ord +  j * nnodes()] += s->operator()(j, n.second.ord) * current_prob/
                    baseline_likelihood;
        
    }

    // Coercing expected to a list vector
    std::vector< std::vector< double > > res(nnodes());
    std::vector< double > zerovec(nfuns(), 0.0);
    for (auto & n: nodes)
    {
        res[n.second.ord] = zerovec;
        for (size_t i = 0u; i < nfuns(); ++i)
            res[n.second.ord][i] = expected[n.second.ord +  i * nnodes()];
    }

    return res;

}
#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-meat-predict_exhaust.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/geese-meat-predict_sim.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include "geese-bones.hpp"

#ifndef GEESE_MEAT_PREDICT_SIM_HPP
#define GEESE_MEAT_PREDICT_SIM_HPP 1

inline std::vector< std::vector<double> > Geese::predict_sim(
    const std::vector< double > & par,
    bool use_reduced_sequence,
    size_t nsims
)
{

    INITIALIZED()

    // Preparing
    std::vector< std::vector< size_t > > tmp;

    std::vector< double > zerovec(nfuns(), 0.0);
    std::vector< std::vector< double > > res_vec(nnodes(), zerovec);
    std::vector< int > counts(nnodes(), 0);

    // We will iterate through this list everytime we need to check
    // whether we have all the annotations for the conditional prob.
    auto annotated = this->get_annotated_nodes();

    for (size_t i = 0u; i < nsims; ++i)
    {

        // Generating a sample
        tmp = this->simulate(par);

        for (auto j = nodes.begin(); j != nodes.end(); ++j)
        {
            // Retrieving node
            const Node & n = j->second;

            // Checking we have all matching
            bool includeit = true;
            for (auto & id : annotated)
            {

                // Same node need not to match (since we are not conditionin
                // each node on itself!)
                if (n.id == id)
                    continue;

                const auto & ord     = nodes[id].ord; 
                const auto & n_w_ann = nodes[id].annotations;
                for (size_t f = 0u; f < nfuns(); ++f)
                {
                    // No checking missings
                    if (n_w_ann[f] == 9u)
                        continue;

                    // If this is not matching, then we cannot use it!
                    if (n_w_ann[f] != tmp[ord][f])
                    {
                        includeit = false;
                        break;
                    }

                }

                if (!includeit)
                    break;
            }

            // If it passed the test, then we can use it for counting stuff
            if (!includeit)
                continue;

            for (size_t f = 0u; f < nfuns(); ++f)
                if (tmp[n.ord][f] == 1)
                    ++res_vec[n.ord][f];

            ++counts[n.ord];

        }

    }

    // Once the simulations have finalized, we can then approximate
    // probabilities
    for (size_t i = 0u; i < nnodes(); ++i)
    {
        // printf_barry("We used %i counts for node %i.\n", counts[i], i);
        for (size_t f = 0u; f < nfuns(); ++f)
            res_vec[i][f] /= (static_cast< double >(counts[i]) + 1e-10);
    }
    
    return res_vec;

}


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/geese-meat-predict_sim.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/flock-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


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
    size_t         nfunctions  = 0u;
    bool                 initialized = false;
    
    // Common components
    std::mt19937              rengine;
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
        std::vector< int > &                       parent,
        std::vector< bool > &                      duplication
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
        bool use_reduced_sequence = true
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
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/flock-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/flock-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef GEESE_FLOCK_MEET_HPP 
#define GEESE_FLOCK_MEET_HPP 1

// #include "flock-bones.hpp"

inline size_t Flock::add_data(
    std::vector< std::vector<size_t> > & annotations,
    std::vector< size_t > &              geneid,
    std::vector< int > &                       parent,
    std::vector< bool > &                      duplication
) {

    // Setting up the model
    if (dat.size() == 0u)
    {

        model.set_rengine(&this->rengine, false);

        model.add_hasher(keygen_full);
        
        model.store_psets();

    }
    else
    {

        if (annotations[0u].size() != nfuns())
            throw std::length_error("The number of functions in the new set of annotations does not match that of the first Geese.");

    }

    // Generating the Geese object
    dat.push_back(Geese(annotations, geneid, parent, duplication));

    if (dat.size() == 1u)
        this->nfunctions = dat[0].nfuns();
       
    return dat.size() - 1u;

}

inline void Flock::set_seed(const size_t & s)
{

    this->rengine.seed(s);

}

inline void Flock::init(size_t bar_width)
{

    // For some strange reason, pointing to model during
    // the add_data function changes addresses once its out.
    for (auto& a : dat)
    {

        if (a.delete_support)
            delete a.model;

        a.model          = &model;
        a.delete_support = false;

        if (a.delete_rengine)
            delete a.rengine;

        a.rengine         = &rengine;
        a.delete_rengine  = false;
        
    }

    // Initializing the models.
    if (bar_width > 0u)
    {

        printf_barry("Initializing nodes in Flock (this could take a while)...\n");
        barry::Progress prog_bar(this->ntrees(), bar_width);
        for (auto& d : dat)
        {

            d.init(0u);
            prog_bar.next();

        }

        prog_bar.end();

    }
    else
    {

        for (auto& d : dat)
            d.init(0u);

    }

    this->initialized = true;
    
}

inline PhyloCounters * Flock::get_counters()
{

    if (dat.size() == 0u)
        throw std::logic_error("The flock has no data yet.");

    return this->model.get_counters();

}

inline PhyloSupport *  Flock::get_support_fun()
{

    return this->model.get_support_fun();

}

inline std::vector< std::vector< double > > *  Flock::get_stats_support()
{

    return this->model.get_stats_support();

}

inline std::vector< std::vector< double > > *  Flock::get_stats_target()
{

    return this->model.get_stats_target();

}

inline PhyloModel *  Flock::get_model()
{

    return &this->model;

}

inline double Flock::likelihood_joint(
    const std::vector< double > & par,
    bool as_log,
    bool use_reduced_sequence
)
{

    INITIALIZED()

    double ans = as_log ? 0.0: 1.0;

    if (as_log) {

        for (auto& d : this->dat) 
            ans += d.likelihood(par, as_log, use_reduced_sequence);

    }
    else
    {

        for (auto& d : this->dat) 
            ans *= d.likelihood(par, as_log, use_reduced_sequence);
            
    }
    
    return ans;

}

inline size_t Flock::nfuns() const noexcept
{

    return this->nfunctions;

}

inline size_t Flock::ntrees() const noexcept
{

    return this->dat.size();

}

inline std::vector< size_t > Flock::nnodes() const noexcept
{

    std::vector< size_t > res;

    res.reserve(this->ntrees());

    for (const auto& d : dat)
        res.push_back(d.nnodes());

    return res;

}

inline std::vector< size_t > Flock::nleafs() const noexcept
{

    std::vector< size_t > res;

    res.reserve(this->ntrees());

    for (const auto& d : dat)
        res.push_back(d.nleafs());

    return res;

}

inline size_t Flock::nterms() const
{

    INITIALIZED()
    return model.nterms() + this->nfuns();

}

inline size_t Flock::support_size() const noexcept
{

    return this->model.support_size();

}

inline std::vector< std::string > Flock::colnames() const
{

    return this->model.colnames();

}

inline size_t Flock::parse_polytomies(
    bool verb,
    std::vector< size_t > * dist
) const noexcept
{

    size_t ans = 0;

    int i = 0;

    for (const auto & d : dat)
    {

        if (verb)
            printf_barry("Checking tree %i\n", i);

        size_t tmp = d.parse_polytomies(verb, dist);

        if (tmp > ans)
            ans = tmp;

    }

    return ans;

}

inline void Flock::print() const 
{

    // Information relevant to print:
    // - Number of phylogenies
    // - Number of functions
    // - Total number of annotations.

    // Computing total number of annotations and events
    size_t nzeros = 0u;

    size_t nones  = 0u;

    size_t ndpl   = 0u;

    size_t nspe   = 0u;

    for (const auto & tree : this->dat)
    {
        nzeros += tree.n_zeros;
        nones  += tree.n_ones;
        ndpl   += tree.n_dupl_events;
        nspe   += tree.n_spec_events;
        
    }

    printf_barry("FLOCK (GROUP OF GEESE)\nINFO ABOUT THE PHYLOGENIES\n");
    
    printf_barry("# of phylogenies         : %li\n", ntrees());
    
    printf_barry("# of functions           : %li\n", nfuns());
    
    printf_barry("# of ann. [zeros; ones]  : [%li; %li]\n", nzeros, nones);
    
    printf_barry("# of events [dupl; spec] : [%li; %li]\n", ndpl, nspe);
    
    printf_barry("Largest polytomy         : %li\n", parse_polytomies(false));
    
    printf_barry("\nINFO ABOUT THE SUPPORT\n");
    
    return this->model.print();

}

inline Geese* Flock::operator()(size_t i, bool check_bounds)
{

    if (check_bounds && i >= ntrees())
        throw std::logic_error("Geese not found in the flock (out of range).");

    return &dat[i];

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/flock-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



namespace phylocounters {
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/models/geese/counters.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H 1



/**
 * @ingroup counting
 * @details Details about the available counters for `PhyloArray`
 * objects can be found in the \ref counters-phylo section.
 */
///@{
#define MAKE_DUPL_VARS() \
    bool DPL = Array.D_ptr()->duplication; \
    size_t DATA_AT = data[0u];

#define IS_EITHER()      (DATA_AT == Geese::etype_either)
#define IS_DUPLICATION() ((DATA_AT == Geese::etype_duplication) & (DPL))
#define IS_SPECIATION()  ((DATA_AT == Geese::etype_speciation) & (!DPL))

#define IF_MATCHES() MAKE_DUPL_VARS() \
    if (IS_EITHER() | IS_DUPLICATION() | IS_SPECIATION())
#define IF_NOTMATCHES() MAKE_DUPL_VARS() \
    if (!IS_EITHER() & !IS_DUPLICATION() & !IS_SPECIATION())


/**
 * @brief Extension of a simple counter.
 * 
 * It allows specifying extra arguments, in particular, the corresponding
 * sets of rows to which this statistic may be relevant. This could be important
 * in the case of, for example, counting correlation type statistics between
 * function 1 and 2, and between function 1 and 3.
 * 
 * 
 */
#define PHYLO_COUNTER_LAMBDA(a) barry::Counter_fun_type<PhyloArray, PhyloCounterData> a = \
    [](const PhyloArray & Array, size_t i, size_t j, PhyloCounterData & data)

#define PHYLO_RULE_DYN_LAMBDA(a) barry::Rule_fun_type<PhyloArray, PhyloRuleDynData> a = \
    [](const PhyloArray & Array, size_t i, size_t j, PhyloRuleDynData & data)

#define PHYLO_CHECK_MISSING() if (Array.D_ptr() == nullptr) \
    throw std::logic_error("The array data is nullptr."); \
    
inline std::string get_last_name(size_t d) {return ((d == 1u)? " at duplication" : ((d == 0u)? " at speciation" : ""));}

/**
 * @weakgroup counters-phylo Phylo counters
 * @brief Counters for phylogenetic modeling
 * @param counters A pointer to a `PhyloCounters` object (`Counters`<`PhyloArray`, `PhyloCounterData`>).
 */
///@{
// -----------------------------------------------------------------------------
/**
 * @brief Overall functional gains
 * @details Total number of gains (irrespective of the function).
 */
inline void counter_overall_gains(
    PhyloCounters * counters,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        
        return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {
        IF_NOTMATCHES()
            return 0.0;
      
        return Array.D_ptr()->states[i] ? 0.0 : 1.0;
      
    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Overall gains" + get_last_name(duplication)
    );

    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Functional gains for a specific function (`nfun`).
 */
inline void counter_gains(
    PhyloCounters * counters,
    std::vector<size_t> nfun,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        IF_NOTMATCHES()
            return 0.0;

        double ngains = 0.0;
        auto   k = data[1u];
        auto   s = Array.D_ptr()->states[k];

        if (s)
            return 0.0;

        for (auto o = 0u; o < Array.ncol(); ++o)
        {
            if (!s && (Array(k,o) == 1u))
                ngains += 1.0;
        }
        
        return ngains;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is there any gain?
        if (Array.D_ptr()->states[i])
            return 0.0;

        IF_MATCHES()
            return (i == data[1u]) ? 1.0 : 0.0;
        
        return 0.0;

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            PhyloCounterData({duplication, i}),
            "Gains " + std::to_string(i) + get_last_name(duplication)
        );
    
    return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief k genes gain function nfun
 */
inline void counter_gains_k_offspring(
    PhyloCounters * counters,
    std::vector<size_t> nfun,
    size_t k = 1u,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

      PHYLO_CHECK_MISSING();
      return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this relevant?
        if (i != data[1u])
            return 0.0;

        IF_NOTMATCHES()
            return 0.0;

        // Is there any gain?
        if (Array.D_ptr()->states[i])
            return 0.0;

        // Making the counts
        int counts = 0;
        for (size_t k = 0u; k < Array.ncol(); ++k)
            if (k != j)
            {
                if (Array(i, k, false) == 1u)
                    ++counts;
            }

        // Three cases: base on the diff
        int diff = static_cast<int>(data[2u]) - counts + 1;
        // (a) counts were 1 below k, then +1
        if (diff == 1)
            return -1.0;
            // (b) counts were equal to k, then -1
        else if (diff == 0)
        {
            return 1.0;
        } else 
            // (c) Otherwise, nothing happens
            return 0.0;
      

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            PhyloCounterData({duplication, i, k}),
            std::to_string(k) + " genes gain " + std::to_string(i) +
                get_last_name(duplication)
        );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Keeps track of how many genes are changing (either 0, 1, or 2 if dealing
 * with regular trees.)
 */
inline void counter_genes_changing(
    PhyloCounters * counters,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        for (auto s : Array.D_ptr()->states)
        {

            if (s) 
                // Yup, we are loosing a function, so break
                return static_cast<double>(Array.ncol());
            
        }

        return 0.0;
      

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;

        // Need to check the other functions
        for (size_t k = 0u; k < Array.nrow(); ++k)
        {

            // Nah, this gene was already different.
            if ((k != i) && (Array.D_ptr()->states[k] != (Array(k, j, false) == 1u)))
                return 0.0;
            

        }

        // Nope, this gene is now matching its parent, so we need to 
        // take it out from the count of genes that have changed.
        return Array.D_ptr()->states[i] ? -1.0 : 1.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Num. of genes changing" + get_last_name(duplication)
    );

    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Keeps track of how many pairs of genes preserve pseudostate.
 */
inline void counter_preserve_pseudogene(
    PhyloCounters * counters,
    size_t nfunA,
    size_t nfunB,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        if (Array.D_ptr()->states[data[1u]] || Array.D_ptr()->states[data[2u]])
            return 0.0;

        double n = static_cast<double>(Array.ncol());
        return n * (n - 1.0) / 2.0;
      

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;

        auto nfunA = data[1u];
        auto nfunB = data[2u];

        if ((i != nfunA) & (i != nfunB))
            return 0.0;

        if (Array.D_ptr()->states[data[1u]] || Array.D_ptr()->states[data[2u]])
            return 0.0;

        size_t k = (i == nfunA) ? nfunB : nfunA;

        if (Array(k, j) == 1u)
            return 0.0;

        double res = 0.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {
            if (off == j)
                continue;

            if ((Array(i, off) == 0u) && (Array(k, off) == 0u))
                res -= 1.0;

        }

        return res;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Preserve pseudo gene (" + 
        std::to_string(nfunA) + ", " +
        std::to_string(nfunB) + ")" + get_last_name(duplication)
    );

    
    return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief Keeps track of how many genes are changing (either 0, 1, or 2 if dealing
 * with regular trees.)
 */
inline void counter_prop_genes_changing(
    PhyloCounters * counters,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        for (auto s : Array.D_ptr()->states)
        {
            if (s)
                return 1.0;
        }
        
        return 0.0;
      
    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;
        
        // Setup
        bool j_diverges = false;
        const std::vector< bool > & par_state = Array.D_ptr()->states;

        for (size_t f = 0u; f < Array.nrow(); ++f)
        {

            // Was the gene annotation different from the parent?
            if (par_state[f] != (Array(f,j) == 1u))
            {
                j_diverges = true;
                break;
            }

        }


        bool j_used_to_diverge = false;
        for (size_t f = 0u; f < Array.nrow(); ++f)
        {

            if (f == i)
            {
                if (par_state[f])
                {
                    j_used_to_diverge = true;
                    break;
                }
            }
            else
            {

                if (par_state[f] != (Array(f,j) == 1u))
                {
                    j_used_to_diverge = true;
                    break;
                }

            }

        }

        // Case 1: j hasn't changed
        if ((!j_used_to_diverge & !j_diverges) | (j_used_to_diverge & j_diverges))
            return 0.0;
        // Case 2: j NOW diverges
        else if (j_diverges)
            return 1.0/Array.ncol();
        // Case 3: j USED to diverge, so no more
        else
            return -1.0/Array.ncol();

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Proportion of genes changing" + get_last_name(duplication)
    );

    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Overall functional loss
 */
inline void counter_overall_loss(
    PhyloCounters * counters,
    size_t duplication = Geese::etype_default
    )
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        if (!Array.D_ptr()->states[i])
            return 0.0;

        IF_MATCHES()
            return -1.0;
        else 
            return 0.0;
        
    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        IF_NOTMATCHES()
            return 0.0;
        
        double res = 0.0;
        for (auto s : Array.D_ptr()->states)
            if (s)
                res += 1.0;

        return res * static_cast<double>(Array.ncol());

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Overall loses" + get_last_name(duplication)
    );
    
    return;

}

// -----------------------------------------------------------------------------
/**
 * @brief Cap the number of functions per gene
 */
inline void counter_maxfuns(
    PhyloCounters * counters,
    size_t            lb,
    size_t            ub,
    size_t duplication = Geese::etype_default
 )
 {

    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();    

        IF_NOTMATCHES()
            return 0.0;

        // At first, all are zero, so we need to check if the lower
        // bound is zero
        if (data[1u] == 0)
            return static_cast<double>(Array.ncol());
        
        return 0.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        int count = Array.colsum(j);
        int ub    = data[2u];
        
        // It now matches
        if (count == static_cast<int>(data[1u]))
            return 1.0;

        // Was within, but now outside
        if (count > ub && ((count - ub) == 1))
            return -1.0;

        // Otherwise nothing happens.
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, lb, ub}),
        "Genes with [" + std::to_string(lb) + ", " + std::to_string(ub) +
            "] funs" + get_last_name(duplication)
    );
    
    return;
  
}
  
// -----------------------------------------------------------------------------
/**
 * @brief Total count of losses for an specific function.
 */
inline void counter_loss(
    PhyloCounters * counters,
    std::vector<size_t> nfun,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        if (!Array.D_ptr()->states[i])
            return 0.0;
        
        return (i == data[1u]) ? -1.0 : 0.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;
        
        auto f = data[1u];

        if (!Array.D_ptr()->states[f])
            return 0.0;
        
        return static_cast<double>(Array.ncol());

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            PhyloCounterData({duplication, i}),
            "Loss " + std::to_string(i) + get_last_name(duplication)
        );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Total number of changes. Use this statistic to account for "preservation"
 */
inline void counter_overall_changes(
    PhyloCounters * counters,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        if (Array.D_ptr()->states[i])
            return -1.0;
        else 
            return 1.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        IF_NOTMATCHES()
            return 0.0;

        PHYLO_CHECK_MISSING();


        // Since we start with all the array at zero,
        // As many chances to change as offspring
        double noff   = static_cast<double> (Array.ncol());
        double counts = 0.0;
        for (size_t k = 0u; k < Array.nrow(); ++k)
            if (Array.D_ptr()->states[k])
                counts += noff;

        return counts;

        

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Overall changes" + get_last_name(duplication)
    );
    
    
    return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief Total count of Sub-functionalization events.
 * @details It requires to specify data = {funA, funB}
 */
inline void counter_subfun(
    PhyloCounters * counters,
    size_t nfunA,
    size_t nfunB,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;

        auto funA = data[1u];
        auto funB = data[2u];
        
        // Are we looking at either of the relevant functions?
        if ((funA != i) && (funB != i))
            return 0.0;
        
        // Are A and B existant? if not, no change
        if (!Array.D_ptr()->states[funA] || !Array.D_ptr()->states[funB])
            return 0.0;
        
        // Figuring out which is the first (reference) function
        size_t other = (i == funA)? funB : funA;
        double res = 0.0;
        // There are 4 cases: (first x second) x (had the second function)
        if (Array(other, j, false) == 1u)
        { 
          
            for (size_t off = 0u; off < Array.ncol(); ++off)
            {
                
                // Not on self
                if (off == j)
                    continue;
                
                if ((Array(i, off, false) == 1u) && (Array(other, off, false) == 0u))
                    res -= 1.0;
                
            }
          
        } else {
          
            for (size_t off = 0u; off < Array.ncol(); ++off)
            {
              
                // Not on self
                if (off == j)
                    continue;
                
                if ((Array(i, off, false) == 0u) && (Array(other, off, false) == 1u))
                    res += 1.0;
              
            }
          
        }
        
        return res;

    };

    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Subfun between " + std::to_string(nfunA) + " and " +
            std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Co-evolution (joint gain or loss)
 * @details Needs to specify pairs of functions (`nfunA`, `nfunB`).
 */
inline void counter_cogain(
    PhyloCounters * counters,
    size_t nfunA,
    size_t nfunB,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        auto d1 = data[1u];
        auto d2 = data[2u];
      
        // Is the function in scope relevant?
        if ((i != d1) && (i != d2))
            return 0.0;
        
        // None should have it
        if (!Array.D_ptr()->states[d1] && !Array.D_ptr()->states[d2])
        {

            size_t other = (i == d1)? d2 : d1;

            if (Array(other, j, false) == 1u)
                return 1.0;
            else
                return 0.0;

        } else 
            return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Co-gains " + std::to_string(nfunA) + " & " + std::to_string(nfunB) +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/** @brief Longest branch mutates (either by gain or by loss) */
inline void counter_longest(
    PhyloCounters * counters,
    size_t duplication = Geese::etype_default
    )
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // Figuring out which match
        std::vector< bool> is_longest(Array.ncol(), false);
        bool j_mutates = false;
        int nmutate = 0;
        int nmutate_longest = 0;

        auto states  = Array.D_ptr()->states;
        
        for (auto off = 0u; off < Array.ncol(); ++off)
        {

            // On the fly, figuring out if it is longest
            for (auto & l : data)
                if (l == off)
                    is_longest[off] = true;

            for (auto f = 0u; f < Array.nrow(); ++f)
            {
                if ((Array(f, off) == 1u) != states[f])
                {
                    
                    // If it happens that j != off and is not longest
                    // then return 0 (a not longest was mutating prev)
                    if (is_longest[off] && (off != j))
                        return 0.0;

                    if (off == j)
                        j_mutates = true;

                    if (is_longest[j])
                        nmutate_longest++;
                    else
                        nmutate++;

                    break;
                }

            }
        }

        // There was already more than one in difference
        // so nothing to change
        if (std::fabs(nmutate - nmutate_longest) > 1)
            return 0.0;

        // Figuring out previously
        bool j_mutates_prev = false;
        for (auto f = 0u; f < Array.nrow(); ++f)
        {
            // Checking the previous function... was it
            // different before?
            if ((f == i) && states[i])
            {
                j_mutates_prev = true;
                break;
            }
            else if ((Array(f, j) == 1u) != states[f])
            {
                j_mutates_prev = true;
                break;
            }

        }

        // Adjusting the previous count
        auto nmutate_prev         = nmutate;
        auto nmutate_longest_prev = nmutate_longest;
        if (j_mutates & !j_mutates_prev)
        {
            if (is_longest[j])
                nmutate_longest_prev--;
            else
                nmutate_prev--;
        }
        else if (!j_mutates & j_mutates)
        {
            if (is_longest[j])
                nmutate_longest_prev++;
            else
                nmutate_prev++;

        }
        
        // Just compute the change statistic directly
        return
            ( ((nmutate == 0) & (nmutate_longest > 0)) ? 1.0 : 0.0 ) +
            ( ((nmutate_prev == 0) & (nmutate_longest_prev > 0)) ? 1.0 : 0.0 );

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        
        if (Array.D_ptr()->blengths.size() != Array.ncol())
            throw std::logic_error(
                "longest should be initialized with a vec of size Array.ncol()."
            );
          
        // Finding the longest branch (or branches) --
        size_t longest_idx = 0u;
        double diff      = 0.0;
        data.reserve(Array.ncol()); 
        data.push_back(0u);
        for (size_t ii = 1u; ii < Array.ncol(); ++ii)
        {
            
            diff = Array.D_ptr()->blengths[longest_idx] - Array.D_ptr()->blengths[ii];
            if (diff > 0.0)
                continue;
            else if (diff < 0.0)
            {

                data.empty();
                data.push_back(ii);
                longest_idx = ii;

            }
            else if (diff == 0.0)
                data.push_back(ii);
            
        }

        data.shrink_to_fit();
        
        if (data.size() == 0u)
            throw std::logic_error("The data on the longest branch has size 0.");
        
        // Starting the counter, since all in zero, then this will be equal to
        // the number of functions in 1 x number of longest branches
        for (size_t ii = 0u; ii < Array.nrow(); ++ii)
        {
            
            if (Array.D_ptr()->states[ii])
                return (1.0 * static_cast<double>(data.size()));

        }
        
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Longest branch mutates" + get_last_name(duplication)
    );
    
    return;
  
}

//------------------------------------------------------------------------------
/**
 * @brief Total number of neofunctionalization events 
 * @details Needs to specify pairs of function.
 */
inline void counter_neofun(
    PhyloCounters * counters,
    size_t nfunA,
    size_t nfunB,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;
        
        auto funA = data[1u];
        auto funB = data[2u];

        // Is the function in scope relevant?
        if ((i != funA) && (i != funB))
            return 0.0;
        
        // Checking if the parent has both functions
        size_t other = (i == funA)? funB : funA;
        bool parent_i     = Array.D_ptr()->states[i];
        bool parent_other = Array.D_ptr()->states[other];
        
        if (!parent_i & !parent_other) 
            return 0.0;
        else if (parent_i & parent_other) 
            return 0.0;
        
        // Figuring out which is the first (reference) function
        double res = 0.0;
        
        if (Array(other, j) == 0u)
        {


            for (auto off = 0u; off < Array.ncol(); ++off)
                if ((Array(i,off) == 0) && (Array(other,off) == 1))
                    res += 1.0;

        }
        else
        {

            for (auto off = 0u; off < Array.ncol(); ++off)
                if ((Array(i,off) == 1) && (Array(other,off) == 0))
                    res -= 1.0;
                
        }
             
        return res;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Neofun between " + std::to_string(nfunA) + " and " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
  
}

//------------------------------------------------------------------------------
/**
 * @brief Total number of neofunctionalization events 
 * sum_u sum_{w < u} [x(u,a)*(1 - x(w,a)) + (1 - x(u,a)) * x(w,a)]
 * change stat: delta{x(u,a): 0->1} = 1 - 2 * x(w,a)
 */
inline void counter_pairwise_neofun_singlefun(
    PhyloCounters * counters,
    size_t nfunA,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;
        
        // Is the function in scope relevant?
        if (i != data[1u])
            return 0.0;
        
        // Checking if the parent has the function
        if (Array.D_ptr()->states[i])
            return 0.0;
        
        // Figuring out which is the first (reference) function
        double res = 0.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {

            if (off == j)
                continue;

            if ((Array(i, off) == 0))
                res += 1.0;
            else 
                res -= 1.0;

        }

        return res;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA}),
        "Pairwise neofun function " + std::to_string(nfunA) +
        get_last_name(duplication)
    );
    
    return;
  
}

//------------------------------------------------------------------------------
/**
 * @brief Total number of neofunctionalization events 
 * @details Needs to specify pairs of function.
 */
inline void counter_neofun_a2b(
    PhyloCounters * counters,
    size_t nfunA,
    size_t nfunB,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;
        
        const size_t & funA = data[1u];
        const size_t & funB = data[2u];

        // Checking scope
        if ((i != funA) && (i != funB))
            return 0.0;
        
        // Checking the parent doesn't have funA or has funB
        if (!Array.D_ptr()->states[funA] || Array.D_ptr()->states[funB]) 
            return 0.0;

        double res = 0.0;

        if (i == funA)
        {

            if (Array(funB, j) == 0u)
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {

                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 0u) && (Array(funB, off) == 1u))
                        res += 1.0;

                }

            }
            else
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {
                    
                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 1u) && (Array(funB, off) == 0u))
                        res -= 1.0;

                }

            }

        }
        else
        {

            if (Array(funA, j) == 0u)
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {

                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 1u) && (Array(funB, off) == 0u))
                        res += 1.0;

                }

            }
            else
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {
                    
                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 0u) && (Array(funB, off) == 1u))
                        res -= 1.0;

                }

            }

        }

        return res;
        
    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Neofun from " + std::to_string(nfunA) + " to " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
    
}

// -----------------------------------------------------------------------------
/**
 * @brief Function co-opting
 * @details Function co-opting of functions A and B happens when, for example,
 * function B is gained as a new featured leveraging what function A already does;
 * without losing function A. The sufficient statistic is defined as follows:
 * \f[
 * x_{pa}(1 - x_{pb})\sum_{i<j}\left[x_{ia}^p(1 - x_{ib}^p)x_{ja}^px_{jb}^p + x_{ja}^p(1 - x_{jb}^p)x_{ia}^px_{ib}^p\right]
 * \f]
 * This algorithm implements the change statistic.
 */
inline void counter_co_opt(
    PhyloCounters * counters,
    size_t nfunA,
    size_t nfunB, 
    size_t duplication = Geese::etype_default
) {
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    { 

        // Checking whether this is for duplication or not
        IF_NOTMATCHES()
            return 0.0;
        
        const size_t funA = data[1u];
        const size_t funB = data[2u];

        // If the change is out of scope, then nothing to do
        if ((i != funA) & (i != funB))
            return 0.0;

        // If the parent does not have the initial state, then it makes no sense
        if ((!Array.D_ptr()->states[funA]) || Array.D_ptr()->states[funB])
            return 0.0;

        // Checking whether function A or function B changed
        if (i == funA) {

            // What was the state of the other function? If B is present, then
            // nothing changes.
            if (Array(funB, j, false) == 1u) 
                return 0.0;

            // Iterating through the sibs
            double res = 0.0;
            for (auto c = 0u; c < Array.ncol(); ++c)
                if ((c != j) && (Array(funA, c, false) == 1u) && (Array(funB, c, false) == 1u))
                    res += 1.0;

            return res;

        } else {

            // What was the state of the other function? If A is not present, then
            // nothing changes.
            if (Array(funA, j, false) == 0u) 
                return 0.0;

            // Iterating through the sibs
            double res = 0.0;
            for (auto c = 0u; c < Array.ncol(); ++c)
                if ((c != j) && (Array(funA, c, false) == 1u))
                    res += (Array(funB, c, false) == 0u) ? 1.0 : -1.0;

            return res;

        }

        

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        if (data.size() != 3u)
            throw std::length_error("The counter data should be of length 2.");

        if (data[1u] == data[2u])
            throw std::logic_error("Functions A and B should be different from each other.");

        if (data[1u] >= Array.nrow())
            throw std::length_error("Function A in counter out of range.");

        if (data[2u] >= Array.nrow())
            throw std::length_error("Function B in counter out of range.");

        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Coopt of " + std::to_string(nfunA) + " by " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Indicator function. Equals to one if \f$k\f$ genes changed and zero
 * otherwise.
 */
inline void counter_k_genes_changing(
    PhyloCounters * counters,
    size_t k,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        for (auto s : Array.D_ptr()->states)
            if (s)
                return Array.ncol() == data[1u] ? 1.0 : 0.0;

        return data[1u] == 0 ? 1.0 : 0.0;
      
    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;
        
        // How many genes diverge the parent
        int              count = 0; 
        bool        j_diverges = false;
        const auto & par_state = Array.D_ptr()->states;

        int k = static_cast<int>(data[1u]);

        for (auto o = 0u; o < Array.ncol(); ++o)
        {

            for (auto f = 0u; f < Array.nrow(); ++f)
            {

                // Was the gene annotation different from the parent?
                if ((Array(f, o) == 1u) != par_state[f])
                {

                    if (o == j)
                        j_diverges = true;

                    count++;
                    break;

                }

            }

        }

        // Counts will only be relevant if (count - k) > 1. Otherwise,
        // having the j gene changed is not relevant
        if (std::abs(count - k) > 1)
            return 0.0;

        // Did it used to diverge?
        bool j_used_to_diverge = false;
        for (auto f = 0u; f < Array.nrow(); ++f)
        {

            if (f == i)
            {
                if (par_state[f]) // Since it is now true, it used to diverge
                {
                    j_used_to_diverge = true;
                    break;
                }
            }
            else
            {

                if (par_state[f] != (Array(f,j) == 1u))
                {
                    j_used_to_diverge = true;
                    break;
                }

            }

        }

        auto count_prev = count;
        // Case 1: j hasn't changed
        if ((!j_used_to_diverge & !j_diverges) | (j_used_to_diverge & j_diverges))
            return 0.0;
        // Case 2: j NOW diverges
        else if (j_diverges)
            count_prev--;
        // Case 3: j USED to diverge
        else
            count_prev++;

        return (count == k ? 1.0 : 0.0) - (count_prev == k ? 1.0 : 0.0);

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, k}),
        std::to_string(k) + " genes changing" + get_last_name(duplication)
    );
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Indicator function. Equals to one if \f$k\f$ genes changed and zero
 * otherwise.
 */
inline void counter_less_than_p_prop_genes_changing(
    PhyloCounters * counters,
    double p,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        for (auto s : Array.D_ptr()->states)
            if (s)
                return data[1u] == 100 ? 1.0 : 0.0;

        // Only one if it was specified it was zero
        return 1.0;
      
    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;
        
        // Setup
        double count = 0.0; ///< How many genes diverge the parent

        bool j_diverges = false;
        const std::vector< bool > & par_state = Array.D_ptr()->states;

        for (size_t o = 0u; o < Array.ncol(); ++o)
        {

            for (size_t f = 0u; f < Array.nrow(); ++f)
            {

                // Was the gene annotation different from the parent?
                if ((Array(f, o) == 1u) != par_state[f])
                {

                    if (o == j)
                        j_diverges = true;

                    count += 1.0;
                    break;

                }

            }

        }


        bool j_used_to_diverge = false;
        for (size_t f = 0u; f < Array.nrow(); ++f)
        {

            if (f == i)
            {
                if (par_state[f])
                {
                    j_used_to_diverge = true;
                    break;
                }
            }
            else
            {

                if (par_state[f] != (Array(f,j) == 1u))
                {
                    j_used_to_diverge = true;
                    break;
                }

            }

        }

        auto count_prev = count;
        // Case 1: j hasn't changed
        if ((!j_used_to_diverge & !j_diverges) | (j_used_to_diverge & j_diverges))
            return 0.0;
        // Case 2: j NOW diverges
        else if (j_diverges)
            count_prev -= 1.0;
        // Case 3: j USED to diverge
        else
            count_prev += 1.0;

        double ncol = static_cast<double>(Array.ncol());
        double p    = static_cast<double>(data[1u]) / 100.0;

        return ((count/ncol) <= p ? 1.0 : 0.0) - ((count_prev/ncol) <= p ? 1.0 : 0.0);

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, static_cast<size_t>(p * 100)}),
        std::to_string(p) + " prop genes changing" + get_last_name(duplication)
    );
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 */
inline void counter_gains_from_0(
    PhyloCounters * counters,
    std::vector< size_t > nfun,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // All must be false
        for (auto s : Array.D_ptr()->states)
        {

            if (s)
                return 0.0;

        }

        // Is this the function?
        if (i != data[1u])
            return 0.0;

        // Now computing the change stats
        double res = static_cast<double>(Array.ncol()) - 1.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {
            if (off  == j)
                continue;

            if (Array(i, off) == 1u)
                res -= 2.0;
        }


        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            PhyloCounterData({duplication, i}),
            "First gain " + std::to_string(i) +
                get_last_name(duplication)
        );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 */
inline void counter_overall_gains_from_0(
    PhyloCounters * counters,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // All must be false
        for (auto s : Array.D_ptr()->states)
        {

            if (s)
                return 0.0;

        }

        return 1.0;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Overall first gains" +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 */
inline void counter_pairwise_overall_change(
    PhyloCounters * counters,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        size_t funpar = Array.D_ptr()->states[i] == 1u;

        // All must be false
        double res = 0.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {
            if (off == j)
                continue;

            if (funpar > Array(i, off))
                res -= 1.0;
            else if (funpar < Array(i, off))
                res += 1.0;
        }
        
        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        double res = 0.0;
        double n   = static_cast<double>(Array.ncol());
        for (auto s : Array.D_ptr()->states)
            if (s)
                res += n * (n - 1.0) / 2.0;

        return res;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Pairs of genes changing" +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 * sum x(a)^3(1-x(b))^3 + x(b)^3(1-x(a))^3 + x(a)^3 * x(b)^3 + (1 - x(a))^3 * (1-x(b))^3
 */
inline void counter_pairwise_preserving(
    PhyloCounters * counters,
    size_t nfunA,
    size_t nfunB,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // Not in the scope
        auto funA = data[1u];
        auto funB = data[2u];
        if ((funA != i) && (funB != i))
            return 0.0;

        size_t k = (funA == i) ? funB : funA;

        bool parent_i = Array.D_ptr()->states[i];
        bool parent_k = Array.D_ptr()->states[k];

        // if (!parent_i & !parent_k)
        //     return 0.0;

        double res = 0.0;
        // Case 1: (0,0)
        if (!parent_i & !parent_k)
        {

            if (Array(k, j) == 1u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 0u) && (Array(k, off) == 0u))
                    res -= 1.0;

            }

        }
        else if (parent_i & !parent_k)
        {

            if (Array(k, j) == 1u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 1u) && (Array(k, off) == 0u))
                    res += 1.0;

            }

        }
        else if (!parent_i & parent_k)
        {

            if (Array(k, j) == 0u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 0u) && (Array(k, off) == 1u))
                    res += 1.0;

            }

        }
        else
        {

            if (Array(k, j) == 0u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 1u) && (Array(k, off) == 1u))
                    res += 1.0;

            }
        }

        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {


        IF_NOTMATCHES()
            return 0.0;
        
        PHYLO_CHECK_MISSING();
        
        double n = static_cast< double >(Array.ncol());
        if (!Array.D_ptr()->states[data[1u]] && !Array.D_ptr()->states[data[2u]])
            return n * (n - 1.0) / 2.0;

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Pariwise preserve (" + std::to_string(nfunA) + ", " +
            std::to_string(nfunB) + ")" +get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 * sum x(a)^3(1-x(b))^3 + x(b)^3(1-x(a))^3 + x(a)^3 * x(b)^3 + (1 - x(a))^3 * (1-x(b))^3
 */
inline void counter_pairwise_first_gain(
    PhyloCounters * counters,
    size_t nfunA,
    size_t nfunB,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // Not in the scope
        auto funA = data[1u];
        auto funB = data[2u];
        if ((funA != i) && (funB != i))
            return 0.0;

        size_t k = (funA == i) ? funB : funA;

        double res = 0.0;
        if (Array(k, j) == 1)
        {

            for (auto off = 0u; off < Array.ncol(); ++off)
            {
                if (off == j)
                    continue;

                if ((Array(i,off) == 0u) && (Array(k,off) == 0u))
                    res -= 1.0;
            }

        }
        else
        {

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 1u))
                {

                    // j: (0,0)\(1,0) -> (1,0)\(1,0), so less 1
                    if (Array(k, off) == 0u)
                        res -= 1.0;

                }
                else
                {

                    if (Array(k, off) == 1u) 
                    // j: (0,0)\(0,1) -> (1,0)\(0,1), so less 1
                        res -= 1.0;
                    else
                    // j: (0,0)\(0,0) -> (1,0)\(0,0), so plus 1
                        res += 1.0;

                }

            }

        }
        

        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "First gain (either " + std::to_string(nfunA) + " or " +
            std::to_string(nfunB) + ")" +get_last_name(duplication)
    );
    
    return;
  
}

///@}

/**
 * @weakgroup rules-phylo Phylo rules
 * @brief Rules for phylogenetic modeling
 * @param rules A pointer to a `PhyloRules` object (`Rules`<`PhyloArray`, `PhyloRuleData`>).
 */
///@{
/**
 * @brief Overall functional gains
 * @param support Support of a model.
 * @param pos Position of the focal statistic.
 * @param lb Lower bound
 * @param ub Upper bound
 * @details 
 * @return (void) adds a rule limiting the support of the model.
 */
inline void rule_dyn_limit_changes(
    PhyloSupport * support,
    size_t pos,
    size_t lb,
    size_t ub,
    size_t duplication = Geese::etype_default
)
{
  
    PHYLO_RULE_DYN_LAMBDA(tmp_rule)
    {

        size_t rule_type = data.duplication;
        if (rule_type != Geese::etype_either)
        {

            if (Array.D_ptr()->duplication & (rule_type != Geese::etype_duplication))
                return true;
            else if (!Array.D_ptr()->duplication & (rule_type != Geese::etype_speciation))
                return true;
                
        }

        if (data() < data.lb)
            return false;
        else if (data() > data.ub)
            return false;
        else
            return true;
      
    };

    // Checking whether the rule makes sense (dupl)
    if (duplication != Geese::etype_either)
    {
        if (support->get_counters()->operator[](pos).data[0] != duplication)
            throw std::logic_error(
                "The rule is not compatible with the duplication type of the model." +
                std::string("The rule is for ") + get_last_name(duplication) +
                std::string(" but the term is for ") + get_last_name(
                    support->get_counters()->operator[](pos).data[0]
                    )
        );
    }
    
    support->get_rules_dyn()->add_rule(
        tmp_rule,
        PhyloRuleDynData(
            support->get_current_stats(),
            pos, lb, ub, duplication
            ),
        std::string("Limiting changes in '") + 
            support->get_counters()->get_names()[pos] +
            "' to [" + std::to_string(lb) + ", " +
            std::to_string(ub) + std::string("]") +
            get_last_name(duplication),        
        std::string("When the support is ennumerated, the number of changes in '") + 
            support->get_counters()->get_names()[pos] +
            std::to_string(pos) + "' is limited to [" + std::to_string(lb) + ", " +
            std::to_string(ub) + "]" +
            get_last_name(duplication)
    );
    
    return;
  
}

///@}

#undef MAKE_DUPL_VARS
#undef IS_EITHER
#undef IS_DUPLICATION
#undef IS_SPECIATION
#undef IF_MATCHES
#undef IF_NOTMATCHES


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/models/geese/counters.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


}


#endif
