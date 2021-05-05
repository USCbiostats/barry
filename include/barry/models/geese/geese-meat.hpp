// #include "geese-bones.hpp"

#ifndef GEESE_MEAT_HPP
#define GEESE_MEAT_HPP 1

inline void Geese::init_node(Node & n) {

    // Creating the phyloarray, nfunctions x noffspring
    n.array = phylocounters::PhyloArray(nfunctions, n.offspring.size());
    std::vector< bool > tmp_state = vector_caster<bool,uint>(n.annotations);
    std::vector< double > blen(n.offspring.size(), 1.0);
    n.array.set_data(
        new phylocounters::NodeData(blen, tmp_state, n.duplication),
        true
    );

    // We initialize all with a zero since, if excluded from the pruning process,
    // We need to set it to one (as the result of the full integration).
    n.subtree_prob.resize(states.size(), 1.0);

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

    if (delete_support)
        delete support;

    if (delete_rengine)
        delete rengine;

    return;
}

inline void Geese::init() {

    // Initializing the model, if it is null
    if (this->support == nullptr) {

        this->support = new phylocounters::PhyloModel();
        this->delete_support = true;
        this->support->set_keygen(keygen_full);
        this->support->store_psets();

    }

    // Checking rseed, this is relevant when dealing with a flock. In the case of
    // flock, both support and rengine are shared.
    if (this->support->get_rengine() == nullptr) 
        this->support->set_rengine(this->rengine, false);

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

    // Resetting the sequence
    for (auto& n: this->nodes) {
        n.second.visited = false;
    }

    // So that others now know it was initialized
    initialized = true;

    return;
}

inline void Geese::inherit_support(const Geese & model_, bool delete_support_) {
    
    if (this->support != nullptr)
        throw std::logic_error("There is already a -support- in this Geese. Cannot set a -support- after one is present.");

    this->support = model_.support;
    this->delete_support = delete_support_;

    // And random number generation
    if (this->delete_rengine) {
        delete this->rengine;
        this->delete_rengine = false;
    }
    
    this->rengine = model_.rengine;
    
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

inline void Geese::calc_reduced_sequence() {

    // The criteria, if none of its decendants is annotated, then we can remove
    // the node from the model
    std::vector< bool > includeit(nodes.size(), false);
    for (auto& i : sequence) {

        Node & n = nodes[i];

        // We will count this at the end
        if (n.is_leaf())
        {

            for (unsigned int k = 0u; k < nfuns(); ++k)
                if (n.annotations[k] != 9u)
                {

                    includeit[i] = true;
                    reduced_sequence.push_back(i);
                    break;

                }

        } else {

            // Checking, am I including any of my offspring?
            for (auto& o : n.offspring) 

                if (includeit[o->id])
                {
                    
                    includeit[i] = true;
                    reduced_sequence.push_back(i);
                    break;

                }

        }

    }

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

inline unsigned int Geese::nfuns() const noexcept {
    return this->nfunctions;
}

inline unsigned int Geese::nnodes() const noexcept {
    return this->nodes.size();
}

inline unsigned int Geese::nleafs() const noexcept {

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
    tmpcount.set_counters(this->support->get_counters());

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
    tmpcount.set_counters(this->support->get_counters());

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
            std::cout << tmparray.D()->states[f] << ", ";
        std::cout << "]; Array:" << std::endl;
        tmparray.print();
        std::cout << "Counts: ";
        for (auto& c : counts)
            std::cout << c << ", ";
        std::cout << std::endl;

    }

    return;

}

inline std::mt19937 * Geese::get_rengine() {
    return this->rengine;
}

inline phylocounters::PhyloCounters * Geese::get_counters() {
    return this->support->get_counters();
}

inline phylocounters::PhyloSupport * Geese::get_support() {
    return this->support->get_support();
}

inline std::vector< std::vector< bool > > Geese::get_states() {
    return this->states;
}


#endif