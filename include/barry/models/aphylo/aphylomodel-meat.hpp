#include "aphylomodel-bones.hpp"

#ifndef APHYLOMODEL_MEAT_HPP
#define APHYLOMODEL_MEAT_HPP 1

APhyloModel::APhyloModel() : model_const(), model_full(), nodes() {
    return;
}

APhyloModel::APhyloModel(
    std::vector< std::vector<unsigned int> > & annotations,
    std::vector< unsigned int > & geneid,
    std::vector< int > &          parent,
    std::vector< bool > &         duplication
) : model_const(), model_full(), nodes() {

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

        if ((parent.at(i) > 0) && (nodes.find(parent.at(i)) == nodes.end())) {

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

            } else {

                // We just need to make sure that we update it!
                nodes[geneid.at(i)].duplication = duplication.at(i);
                nodes[geneid.at(i)].annotations = funs;
                nodes[geneid.at(i)].parent = &nodes[parent.at(i)];

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
                nodes[parent.at(i)].offspring.push_back(
                    &key_off.first->second
                );

                // Adding the parent to the offspring
                key_off.first->second.parent = &nodes[parent.at(i)];

            } else {

                // We just need to make sure that we update it!
                nodes[geneid.at(i)].duplication = duplication.at(i);
                nodes[geneid.at(i)].annotations = funs;

                if (parent.at(i) > 0) {
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

void APhyloModel::init() {

    // Generating the model data -----------------------------------------------
    model_const.set_keygen(keygen_const);
    model_full.set_keygen(keygen_full);

    model_const.add_rule(
            rule_empty_free<phylocounters::PhyloArray,phylocounters::PhyloRuleData>
            );

    model_const.set_counters(&counters);
    model_full.set_counters(&counters);

    model_full.store_psets();

    model_const.set_rengine(&this->rengine, false);
    model_full.set_rengine(&this->rengine, false);

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
        if (!iter.second.is_leaf()) {

            // Creating the phyloarray, nfunctions x noffspring
            iter.second.array = phylocounters::PhyloArray(nfunctions, iter.second.offspring.size());
            iter.second.probabilities.resize(pset.size(), 0.0);

            // Adding the data, first through functions
            for (unsigned int k = 0u; k < nfunctions; ++k) {

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

                iter.second.arrays.push_back(
                    phylocounters::PhyloArray(iter.second.array, true));
                iter.second.arrays.at(i).set_data(
                    new phylocounters::NodeData(blen, s),
                    true
                );

                // Once the array is ready, we can add it to the model
                iter.second.idx_cons.push_back(
                    model_const.add_array(iter.second.arrays.at(i))
                    );
                iter.second.idx_full.push_back(
                    model_full.add_array(iter.second.arrays.at(i++))
                    );

                // model_full.print_stats(0u);

            }
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

    // So that others now know it was initialized
    initialized = true;

    return;
}

void APhyloModel::calc_sequence(Node * n) {

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

std::vector< double > APhyloModel::get_probabilities() const {

    std::vector< double > res;
    res.reserve(
        this->states.size() * nodes.size()
        );
    
    for (auto& i : sequence) {
        for (auto& p : this->nodes.at(i).probabilities)
            res.push_back(p);
    }

    return res;
    
}

unsigned int APhyloModel::nfuns() const {
    return this->nfunctions;
}

unsigned int APhyloModel::nnodes() const {
    return this->nodes.size();
}

unsigned int APhyloModel::nleafs() const {

    unsigned int n = 0u;
    for (auto& iter : this->nodes)
        if (iter.second.is_leaf())
            n++;

    return n;
}

unsigned int APhyloModel::nterms() const {

    INITIALIZED()

    return model_const.nterms() + this->nfuns();
}

#endif