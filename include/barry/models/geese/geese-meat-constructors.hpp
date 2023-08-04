// #include "geese-bones.hpp"

#ifndef GEESE_MEAT_CONSTRUCTORS_HPP
#define GEESE_MEAT_CONSTRUCTORS_HPP 1

inline Geese::Geese() {

    // In order to start...
    this->rengine         = new std::mt19937;
    this->delete_rengine  = true;
    this->model           = new phylocounters::PhyloModel();
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
    this->model           = new phylocounters::PhyloModel();
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
            model = new phylocounters::PhyloModel(*(model_.model));
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

        model = new phylocounters::PhyloModel(*x.model);
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