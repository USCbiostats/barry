#include "geese-bones.hpp"

#ifndef GEESE_MEAT_CONSTRUCTORS_HPP
#define GEESE_MEAT_CONSTRUCTORS_HPP 1

inline Geese::Geese() {

    // In order to start...
    this->counters        = new phylocounters::PhyloCounters();
    this->delete_counters = true;
    this->rengine         = new std::mt19937;
    this->delete_rengine  = true;

    return;
}

inline Geese::Geese(
    std::vector< std::vector<unsigned int> > & annotations,
    std::vector< unsigned int > & geneid,
    std::vector< int > &          parent,
    std::vector< bool > &         duplication
) {

    // In order to start...
    this->counters        = new phylocounters::PhyloCounters();
    this->delete_counters = true;
    this->rengine         = new std::mt19937;
    this->delete_rengine  = true;

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

    // Computing the pruning sequence.
    calc_sequence();
    calc_likelihood_sequence();

    return;

}

inline Geese::Geese(const Geese & model_, bool copy_data) : 
    states(model_.states),
    nfunctions(model_.nfunctions),
    nodes(model_.nodes),
    map_to_nodes(model_.map_to_nodes),
    sequence(model_.sequence),
    likelihood_sequence(model_.likelihood_sequence),
    initialized(model_.initialized) {

    
    // Replicating -------------------------------------------------------------
    if (copy_data) {

        if (model_.rengine != nullptr) {
            rengine = new std::mt19937(*(model_.rengine));
            delete_rengine = true;
        }

        if (model_.counters != nullptr) {
            counters = new phylocounters::PhyloCounters(*(model_.counters));
            delete_counters = true;
        }

        if (model_.support != nullptr) {
            support = new phylocounters::PhyloModel(*(model_.support));
            delete_support = true;
        }

    } else {
        
        if (model_.rengine != nullptr) {
            rengine = model_.rengine;
            delete_rengine = false;
        }

        if (model_.counters != nullptr) {
            counters = model_.counters;
            delete_counters = false;
        }

        if (model_.support != nullptr) {
            support = model_.support;
            delete_support = false;
        }

    }

    // Dealing with the nodes is a bit different -------------------------------
    auto revseq = this->sequence;
    std::reverse(revseq.begin(), revseq.end());

    for (auto& i : revseq) {

        // Leaf do not have offspring
        if (this->nodes[i].is_leaf())
            continue;

        // Clearing offspring
        this->nodes[i].offspring.clear();

        // I cannot directly access the node since, if non existent, it will 
        // create an entry with it (alegedly).
        auto n = model_.nodes.find(i);

        for (const auto& off : n->second.offspring) {
            this->nodes[i].offspring.push_back(&this->nodes[off->id]);
        }

    }

    return;
  
}

// Constructor move
inline Geese::Geese(Geese && x) noexcept :
    rengine(nullptr),
    counters(nullptr),
    support(nullptr),
    states(std::move(x.states)),
    nfunctions(x.nfunctions),
    nodes(std::move(x.nodes)),
    map_to_nodes(std::move(x.map_to_nodes)),
    sequence(std::move(x.sequence)),
    likelihood_sequence(std::move(x.likelihood_sequence)),
    initialized(x.initialized)
{

    // To point or not to point
    if (x.delete_counters) {

        counters = new phylocounters::PhyloCounters(*x.counters);
        delete_counters = true;

    } else {

        counters = x.counters;
        delete_counters = false;

    }

    if (x.delete_rengine) {

        rengine = new std::mt19937(*x.rengine);
        delete_rengine = true;

    } else {

        rengine = x.rengine;
        delete_rengine = false;

    }

    if (x.delete_support) {

        support = new phylocounters::PhyloModel(*x.support);
        delete_support = true;

    } else {

        support = x.support;
        delete_support = false;
        
    }

    // Figuring out if support needs to be updated
    if ((support != nullptr) && (x.delete_support | x.delete_rengine))
        support->set_rengine(rengine, false);

    if ((support != nullptr) && (x.delete_support | x.delete_counters))
        support->set_counters(counters);

    return;

}

// // Copy assignment
// inline Geese & Geese::operator=(const Geese & model_) {

// }

// // Move assignment
// inline Geese & Geese::operator=(Geese && model_) noexcept {

// }

#endif