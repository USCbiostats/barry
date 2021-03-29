#include "geese-bones.hpp"

#ifndef GEESE_MEAT_HPP
#define GEESE_MEAT_HPP 1

inline Geese::Geese() : model(nullptr), nodes() {
    return;
}

inline Geese::Geese(
    std::vector< std::vector<unsigned int> > & annotations,
    std::vector< unsigned int > & geneid,
    std::vector< int > &          parent,
    std::vector< bool > &         duplication
) : model(nullptr), nodes() {

    // Check the lengths
    if (annotations.size() == 0u)
        throw std::logic_error("Annotations is empty");

    nfunctions = annotations.at(0u).size();

    // unsigned int n = annotations.size();
    for (auto& iter : annotations) {
        if (iter.size() != nfunctions)
            throw std::length_error("Not all the annotations have the same length");
    }

    // Grouping up the data by parents -----------------------------------------
    for (unsigned int i = 0u; i < geneid.size(); ++i) {

        // Temp vector with the annotations
        std::vector< unsigned int > funs(annotations.at(i));

        if ((parent.at(i) >= 0) && (nodes.find(parent.at(i)) == nodes.end())) {

            // Adding parent
            auto key_par = nodes.insert({
                parent.at(i),
                Node(parent.at(i), true)
            });

            // Adding offspring
            if (nodes.find(geneid.at(i)) == nodes.end()) {

                auto key_off = nodes.insert({
                    geneid.at(i),
                    Node(geneid.at(i), funs, duplication.at(i))
                    });

                // Adding the offspring to the parent
                key_par.first->second.offspring.push_back(
                    &key_off.first->second
                );

                // Adding the parent to the offspring
                key_off.first->second.parent = &key_par.first->second;
                key_off.first->second.annotations = funs;
                key_off.first->second.duplication = duplication.at(i);

            } else {

                // We just need to make sure that we update it!
                nodes[geneid.at(i)].duplication = duplication.at(i);
                nodes[geneid.at(i)].annotations = funs;
                nodes[geneid.at(i)].parent      = &nodes[parent.at(i)];

                nodes[parent.at(i)].offspring.push_back(
                    &nodes[geneid.at(i)]
                );

            }

        } else {
            // In this case, the parent exists, so we only need to assing the
            // offspring
            // Adding offspring

            // Does the offspring exist?
            if (nodes.find(geneid.at(i)) == nodes.end()) {

                auto key_off = nodes.insert({
                    geneid.at(i),
                    Node(geneid.at(i), funs, duplication.at(i))
                    });

                // Adding the offspring to the parent
                if (parent.at(i) >= 0) {
                    nodes[parent.at(i)].offspring.push_back(
                        &key_off.first->second
                    );

                    // Adding the parent to the offspring
                    key_off.first->second.parent = &nodes[parent.at(i)];
                }

                key_off.first->second.annotations = funs;
                key_off.first->second.duplication = duplication.at(i);

            } else {

                // We just need to make sure that we update it!
                nodes[geneid.at(i)].duplication = duplication.at(i);
                nodes[geneid.at(i)].annotations = funs;

                if (parent.at(i) >= 0) {
                    nodes[geneid.at(i)].parent = &nodes[parent.at(i)];
                    nodes[parent.at(i)].offspring.push_back(
                        &nodes[geneid.at(i)]
                    );
                }

            }
        }

    }

    // init(counters);

    return;

}

inline void Geese::init_node(Node & n) {

    // Creating the phyloarray, nfunctions x noffspring
    n.array = phylocounters::PhyloArray(nfunctions, n.offspring.size());
    std::vector< bool > tmp_state = vector_caster<bool,uint>(n.annotations);
    std::vector< double > blen(n.offspring.size(), 1.0);
    n.array.set_data(
        new phylocounters::NodeData(blen, tmp_state, n.duplication),
        true
    );

    n.subtree_prob.resize(states.size(), 0.0);

    // Adding the data, first through functions
    for (unsigned int k = 0u; k < nfunctions; ++k) {

        // Then through the offspring
        unsigned int j = 0;
        for (auto& o : n.offspring) {

            // If leaf, then it may have an annotation
            if (o->is_leaf()) {
                if (o->annotations.at(k) != 0) {
                    n.array.insert_cell(
                        k, j, o->annotations.at(k), false, false
                        );
                }
            } else {
                // Otherwise, we fill it with a 9.
                n.array.insert_cell(
                    k, j, 9u, false, false
                    );

            }

            ++j;

        }

    }

    // We then need to set the powerset
    if (n.arrays.size() != states.size()) {
        n.arrays.resize(states.size());
        n.narray.resize(states.size());
    }
    
    for (unsigned int s = 0u; s < states.size(); ++s) {

        n.arrays[s] = phylocounters::PhyloArray(n.array, true);
        n.arrays[s].set_data(
            new phylocounters::NodeData(blen, states[s], n.duplication),
            true
        );

        // Once the array is ready, we can add it to the model
        n.narray[s] = support->add_array(n.arrays[s]);

    }

    return;
}

inline Geese::~Geese() {
    if (delete_model)
        delete model;
    return;
}

inline void Geese::init() {

    // Initializing the model, if it is null
    if (this->model == nullptr) {
        model = new phylocounters::PhyloModel();
        this->delete_model = true;
    }

    // Generating the model data -----------------------------------------------
    support->set_keygen(keygen_full);
    support->set_counters(&counters);
    support->store_psets();
    support->set_rengine(&this->rengine, false);

    // All combinations of the function
    phylocounters::PhyloPowerSet pset(nfunctions, 1u);
    pset.calc();

    states.reserve(pset.data.size());
    unsigned int i = 0u;
    for (auto& iter : pset.data) {

        states.push_back(std::vector< bool >(nfunctions, false));

        for (auto iter2 = iter.get_col(0u, false)->begin(); iter2 != iter.get_col(0u, false)->end(); ++iter2)
            states.at(i).at(iter2->first) = true;

        // Adding to map so we can look at it later on
        map_to_nodes.insert({iter.get_col_vec(0u, false), i});

        i++;
    }

    // Iterating throught the nodes
    for (auto& iter : nodes) {

        // Only parents get a node
        if (!iter.second.is_leaf()) 
            this->init_node(iter.second);
            
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

    // So that others now know it was initialized
    initialized = true;

    return;
}

inline void Geese::inherit_support(Geese & model_, bool delete_model_) {
    
    if (this->model != nullptr)
        throw std::logic_error("There is already a model in this Geese. Cannot set a model after one is present.");

    this->model = model_.model;
    this->delete_model = delete_model_;
    return;

}

inline void Geese::set_support(phylocounters::PhyloModel * model_, bool delete_model_) {

    if (this->model != nullptr)
        throw std::logic_error("There is already a model in this Geese. Cannot set a model after one is present.");

    this->model = model_;
    this->delete_model = delete_model_;
    return;

}

inline void Geese::update_annotations(
    unsigned int nodeid,
    std::vector< unsigned int > newann
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
    init_node(*nodes[nodeid].parent);

    return;
}

inline void Geese::calc_sequence(Node * n) {

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

inline std::vector< double > Geese::get_probabilities() const {

    std::vector< double > res;
    res.reserve(
        this->states.size() * nodes.size()
        );
    
    for (auto& i : sequence) {
        for (auto& p : this->nodes.at(i).subtree_prob)
            res.push_back(p);
    }

    return res;
    
}

inline unsigned int Geese::nfuns() const {
    return this->nfunctions;
}

inline unsigned int Geese::nnodes() const {
    return this->nodes.size();
}

inline unsigned int Geese::nleafs() const {

    unsigned int n = 0u;
    for (auto& iter : this->nodes)
        if (iter.second.is_leaf())
            n++;

    return n;
}

inline unsigned int Geese::nterms() const {

    INITIALIZED()

    return support->nterms() + this->nfuns();
}

inline std::vector< std::vector<double> > Geese::observed_counts() {

    // Making room for the output
    std::vector<std::vector<double>> ans;
    ans.reserve(nnodes());

    // Creating counter
    phylocounters::PhyloStatsCounter tmpcount;
    tmpcount.set_counters(&support->counters);

    // Iterating through the nodes
    for (auto& n : nodes) {

        if (n.second.is_leaf()) {
            ans.push_back({});
            continue;
        }

        phylocounters::PhyloArray tmparray(nfuns(), n.second.offspring.size());

        uint j = 0u;
        for (auto& o : n.second.offspring) {
            for (uint k = 0u; k < nfuns(); ++k) {
                if (o->annotations.at(k) != 0) {
                    tmparray.insert_cell(
                        k, j, o->annotations.at(k), false, false
                        );
                }
            }
            ++j;
        }

        std::vector< bool > tmp_state =vector_caster<bool,uint>(n.second.annotations);
        std::vector< double > blen(n.second.offspring.size(), 1.0);
        tmparray.set_data(
            new phylocounters::NodeData(blen, tmp_state, n.second.duplication),
            true
        );

        tmpcount.reset_array(&tmparray);
        ans.push_back(tmpcount.count_all());

    }

    return ans;

}

inline void Geese::print_observed_counts() {

    // Making room for the output
    std::vector<std::vector<double>> ans;
    ans.reserve(nnodes());

    // Creating counter
    phylocounters::PhyloStatsCounter tmpcount;
    tmpcount.set_counters(&support->counters);

    // Iterating through the nodes
    for (auto& n : nodes) {

        if (n.second.is_leaf()) {
            ans.push_back({});
            continue;
        }

        phylocounters::PhyloArray tmparray(nfuns(), n.second.offspring.size());

        uint j = 0u;
        for (auto& o : n.second.offspring) {
            for (uint k = 0u; k < nfuns(); ++k) {
                if (o->annotations.at(k) != 0) {
                    tmparray.insert_cell(
                        k, j, o->annotations.at(k), false, false
                        );
                }
            }
            ++j;
        }

        std::vector< bool > tmp_state =vector_caster<bool,uint>(n.second.annotations);
        std::vector< double > blen(n.second.offspring.size(), 1.0);
        tmparray.set_data(
            new phylocounters::NodeData(blen, tmp_state, n.second.duplication),
            true
        );

        tmpcount.reset_array(&tmparray);
        std::vector< double > counts = tmpcount.count_all();

        // Printing
        auto dpl = n.second.duplication ? "duplication" : "speciation";
        std::cout << "----------\n" <<
            "nodeid: " << n.second.id << " (" << dpl <<
            ")\nstate : [";
        for (uint f = 0u; f < nfuns(); ++f)
            std::cout << tmparray.data->states[f] << ", ";
        std::cout << "]; Array:" << std::endl;
        tmparray.print();
        std::cout << "Counts: ";
        for (auto& c : counts)
            std::cout << c << ", ";
        std::cout << std::endl;

    }

    return;

}

#endif