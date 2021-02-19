#include "../../include/barry/barry.hpp"

#ifndef LEAF_HPP
#define LEAF_HPP 1

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

    // // Free cells
    // for (auto i = 0u; i < array.nrow(); ++i)
    //     for (auto j = 0u; j < array.ncol(); ++j)
    //         dat.push_back((double) array.get_cell(i, j, false));
    
    return dat;
}

class Node {
public:
    unsigned int id;
    PhyloArray array;
    
    std::vector< unsigned int > annotations;         ///< Observed annotations (only defined for leafs)
    std::vector< PhyloArray >   arrays    = {};      ///< Arrays given all possible states
    Node *                      parent    = nullptr; ///< Parent node
    std::vector< Node* >        offspring = {};      ///< Offspring nodes
    std::vector< unsigned int > idx_cons  = {};      ///< Id of the constrained support.
    std::vector< unsigned int > idx_full  = {};
    bool                        visited   = false;          

    /**
     * @brief The probability of observing each state 
     */
    std::vector< double > probabilities; 

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


class Leafs {
public:

    PhyloModel                         model_const;
    PhyloModel                         model_full;
    unsigned int                       nfuns;
    barry::Map< unsigned int, Node >   nodes;
    std::vector< unsigned int >        sequence;
    std::vector< bool >                visited;
    std::vector< std::vector< bool > > states;
    barry::MapVec_type< unsigned int > map_to_nodes;

    Leafs();

    Leafs(
        std::vector< std::vector<unsigned int> > & annotations,
        std::vector< unsigned int > & geneid,
        std::vector< unsigned int> & parent,
        PhyloCounters & counters
        );

    ~Leafs() {};

    void init(PhyloCounters & counters);

    double operator()(std::vector< double > & par, unsigned int & i);
    void tip_prob(const std::vector< double > & par);
    void print();
    void calc_sequence(Node * n = nullptr);
    double likelihood(const std::vector< double > & par);
};

Leafs::Leafs() : model_const(), model_full(), nodes() {
    return;
}

Leafs::Leafs(
    std::vector< std::vector<unsigned int> > & annotations,
    std::vector< unsigned int > & geneid,
    std::vector< unsigned int> & parent,
    PhyloCounters & counters
) : model_const(), model_full(), nodes() {

    // Check the lengths
    if (annotations.size() == 0u)
        throw std::logic_error("Annotations is empty");

    nfuns = annotations.size();

    unsigned int n = annotations.at(0u).size();
    for (auto& iter : annotations) {
        if (iter.size() != n)
            throw std::length_error("Not all the annotations have the same length");
    }

    // Grouping up the data by parents -----------------------------------------
    for (unsigned int i = 0u; i < geneid.size(); ++i) {
        
        // Temp vector with the annotations
        std::vector< unsigned int > funs(nfuns);
        for (unsigned int j = 0u; j < nfuns; ++j)
            funs.at(j) = annotations.at(j).at(i);

        // Does the parent already exists?
        auto iter = nodes.find(parent.at(i));

        if (iter == nodes.end()) {

            // Adding parent
            auto key_par = nodes.insert({
                parent.at(i),
                Node({parent.at(i)}) 
            });

            // Adding offspring
            auto key_off = nodes.insert({
                geneid.at(i),
                Node({geneid.at(i), funs})
                });

            // Adding the offspring to the parent
            key_par.first->second.offspring.push_back(
                &key_off.first->second
            );

            // Adding the parent to the offspring
            key_off.first->second.parent = &key_par.first->second;

        } else { 
            // In this case, the parent exists, so we only need to assing the
            // offspring
            // Adding offspring
            auto key_off = nodes.insert({
                geneid.at(i),
                Node({geneid.at(i), funs})
                });

            // Adding the offspring to the parent
            iter->second.offspring.push_back(
                &key_off.first->second
            );

            // Adding the parent to the offspring
            key_off.first->second.parent = &iter->second;
        }

    }

    init(counters);

    return;

}

void Leafs::init(PhyloCounters & counters) {

    // Generating the model data -----------------------------------------------
    model_const.set_keygen(keygen_const);
    model_full.set_keygen(keygen_full);

    model_const.add_rule(rule_blocked<PhyloArray,PhyloRuleData>);

    model_const.set_counters(&counters);
    model_full.set_counters(&counters);

    model_full.store_psets();

    // All combinations of the function
    PhyloPowerSet pset(nfuns, 1u);
    pset.calc();

    states.reserve(pset.data.size());
    unsigned int i = 0u;
    for (auto& iter : pset.data) {

        states.push_back(std::vector< bool >(nfuns, false));

        for (auto iter2 = iter.get_col(0u)->begin(); iter2 != iter.get_col(0u)->end(); ++iter2)
            states.at(i).at(iter2->first) = true;
        
        // Adding to map so we can look at it later on
        map_to_nodes.insert({iter.get_col_vec(0u), i});

        i++;
    }

    unsigned int narrays = 0u;

    // Iterating throught the nodes
    for (auto& iter : nodes) {

        // Only parents get a node
        if (!iter.second.is_leaf()) {

            // Creating the phyloarray, nfuns x noffspring
            iter.second.array = PhyloArray(nfuns, iter.second.offspring.size());
            iter.second.probabilities.resize(pset.size(), 0.0);

            // Adding the data, first through functions
            for (unsigned int k = 0u; k < nfuns; ++k) {

                // Then through the offspring
                unsigned int j = 0;
                for (auto& o : iter.second.offspring) {
                    
                    // If leaf, then it may have an annotation
                    if (o->is_leaf()) {
                        if (o->annotations.at(k) != 0) {
                            iter.second.array.insert_cell(
                                k, j, o->annotations.at(k), false, false
                                );
                        }
                    } else {
                        // Otherwise, we fill it with a 9.
                        iter.second.array.insert_cell(
                            k, j, 9u, false, false
                            );

                    }

                    ++j;

                }

            }

            // We then need to set the powerset           
            unsigned int i = 0u;
            std::vector< double > blen(iter.second.offspring.size(), 1.0);
            for (auto& s : states) {
                
                iter.second.arrays.push_back(iter.second.array);
                iter.second.arrays.at(i).set_data(
                    new NodeData(blen, s),
                    true
                );

                // Once the array is ready, we can add it to the model
                iter.second.idx_cons.push_back(
                    model_const.add_array(iter.second.arrays.at(i))
                    );
                iter.second.idx_full.push_back(
                    model_full.add_array(iter.second.arrays.at(i++))
                    );

            }

            narrays++;

        }
    }

    // Finally, setting this variable for later, we will need this for generating
    // the pruning sequence.
    visited.resize(nodes.size(), false);

    // Computing the pruning sequence.
    calc_sequence();

    // Resetting the sequence
    for (auto& n: this->nodes) {
        n.second.visited = false;
    }

    return;
}

void Leafs::tip_prob(const std::vector< double > & par) {

    std::vector< double > probs;
    std::vector< unsigned int > ids;
    unsigned int i = 0;

    for (auto& iter : nodes) {
        
        // Counting the probability assuming equal likelihood of each 
        // state
        if (!iter.second.is_leaf()) {

            double tmp = 0.0;
            ids.push_back(iter.first);

            // Iterating through the states
            for (unsigned int j = 0u; j < iter.second.idx_cons.size(); ++j) {

                model_const.likelihood(par, iter.second.idx_cons.at(j));
                model_full.likelihood(par, iter.second.idx_full.at(j));
                
                tmp += (( 
                    model_const.get_norm_const(par, iter.second.idx_cons.at(j)) /
                    model_full.get_norm_const(par, iter.second.idx_full.at(j))
                    ) / iter.second.idx_full.size());
            }
            iter.second.array.print();
            model_full.print_stats(i);

            probs.push_back(tmp);

            ++i;

        }
    }

    // Printing the probabilities
    i = 0u;
    for (auto& iter : probs) {
        printf("Pr(%4i) = %.4f\n", ids[i++], iter);
    }

    return;
}

void Leafs::print() {

    // Overall information
    std::cout << "\n--- General information ----------" << std::endl;
    std::cout << " Number of arrays (full)               : " << this->model_full.size() << std::endl;
    std::cout << " Number of unique models (full)        : " << this->model_full.size_unique() << std::endl;
    std::cout << " Number of arrays (constrained)        : " << this->model_const.size() << std::endl;
    std::cout << " Number of unique models (constrained) : " << this->model_const.size_unique() << std::endl;

    // Node information
    std::cout << "\n--- Information about the genes ---" << std::endl;
    for (auto& iter : this->nodes) {
        std::cout << " Gene id: " << iter.first;
        std::cout << " Parent id: " << iter.second.get_parent();
        std::cout << " offsprings: {";
        for (auto& o : iter.second.offspring) {
            std::cout << o->id;
            
            if (o->id != (iter.second.offspring.at(iter.second.offspring.size() - 1u)->id))
                std::cout << ", ";

        }
        std::cout << "}" << std::endl;
    }

    std::cout << "\n--- Raw probabilities -------------" << std::endl;
    std::vector< double > p(model_const.counters.size(), 1);
    this->tip_prob(p);

    return;
}

void Leafs::calc_sequence(Node * n) {

    if (sequence.size() == nodes.size())
        return;

    // First iteration
    if (n == nullptr) {
    
        // pointing to something
        n = &(nodes.begin()->second);

    }

    // Here before?
    if (n->visited) 
        return;

    n->visited = true;

    if (!n->is_leaf()) {

        // iterating over its offspring, only if not there before
        for (auto& it : n->offspring) {
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


double Leafs::likelihood(const std::vector< double > & par) {

    // Following the prunning sequence
    for (auto& i : this->sequence) {

        // We cannot compute probability at the leaf, we need to continue
        if (this->nodes[i].is_leaf())
            continue;

        // If it is located right before a leaf, then the computations
        // are a bit different
        if (this->nodes[i].offspring.at(0u)->is_leaf()) {

            // Iterating through the different parent states
            for (unsigned int s = 0u; s < states.size(); ++s) {

                // Update the normalizing constants
                double numer = model_const.get_norm_const(par, nodes[i].idx_cons[s]);
                double denom = model_full.get_norm_const(par, nodes[i].idx_full[s]);

                // Computing the probability at "leaf" level
                nodes[i].probabilities[s] = numer/denom;

            }

        } else {
            // Iterating through states
            
            for (unsigned int s = 0u; s < states.size(); ++s) {

                // Starting the prob
                double totprob = 0.0;
                
                // Retrieving the sets of arrays
                const std::vector< PhyloArray > * psets = model_full.get_pset(
                    nodes[i].idx_full[s]
                    );

                // Summation over all possible values of X
                for (auto x = psets->begin(); x != psets->end(); ++x) {

                    // Extracting the possible values of each offspring
                    double off_mult = 1.0;
                    for (auto o = 0u; o < x->ncol(); ++o) {

                        // First, getting what is the corresponding state
                        unsigned int loc = nodes[i].offspring[o]->idx_full[map_to_nodes[x->get_col_vec(o)]];
                        off_mult *= nodes[i].offspring[o]->probabilities[loc];

                    }

                    // Multiplying by P(x|x_n)
                    off_mult *= model_full.likelihood(
                        par,
                        nodes[i].idx_full[s]
                        );

                    // Adding to the total probabilities
                    totprob += off_mult;

                }

                // Setting the probability at the node
                nodes[i].probabilities[s] = totprob;

            }

            

        }


    }

    return 0.0;

}

#endif