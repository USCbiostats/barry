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
    for (unsigned int k = 0u; k < nfunctions; ++k)
    {

        // Then through the offspring
        unsigned int j = 0;
        for (auto& o : n.offspring) {

            // If leaf, then it may have an annotation
            if (o->is_leaf())
            {

                if (o->annotations[k] != 0)
                    n.array.insert_cell(k, j, o->annotations[k], false, false);

            } else {
                // Otherwise, we fill it with a 9.
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

    }
    
    for (unsigned int s = 0u; s < states.size(); ++s)
    {

        n.arrays[s] = phylocounters::PhyloArray(n.array, true);
        n.arrays[s].set_data(
            new phylocounters::NodeData(blen, states[s], n.duplication),
            true
        );

        // Once the array is ready, we can add it to the model
        n.narray[s] = model->add_array(n.arrays[s]);

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

inline void Geese::init(bool verb) {

    // Initializing the model, if it is null
    if (this->model == nullptr)
    {

        this->model = new phylocounters::PhyloModel();
        this->delete_support = true;
        this->model->set_keygen(keygen_full);
        this->model->store_psets();

    }

    // Checking rseed, this is relevant when dealing with a flock. In the case of
    // flock, both model and rengine are shared.
    if (this->model->get_rengine() == nullptr) 
        this->model->set_rengine(this->rengine, false);

    // All combinations of the function
    phylocounters::PhyloPowerSet pset(nfunctions, 1u);
    pset.calc();

    states.reserve(pset.data.size());
    unsigned int i = 0u;
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

    if (verb)
        printf_barry("Initializing nodes in Geese (this could take a while)...\n");

    double width   = 73.0;
    int n_internal = this->nnodes() - this->nleafs();
    int n_steps    = (n_internal > width) ?
        static_cast<int>(floor(static_cast<double>(n_internal) / width)) : n_internal;
    int step_size  = static_cast<int>(n_internal / n_steps);
    int n_bars     = static_cast<int>(std::max(1.0, floor(width / static_cast<double>(n_steps))));

    // Iterating throught the nodes
    int node_i = 0;
    for (auto& iter : nodes) {

        // Only parents get a node
        if (!iter.second.is_leaf())
        {
            this->init_node(iter.second); 

            if (verb & !(node_i++ % step_size))
            {
                for (int j = 0; j < n_bars; ++j)
                    printf_barry("|");
            }

        }
        
    }

    // Adding the extra bars, if needed
    if (verb)
    {

        int reminder = static_cast<int>(width) - n_bars * n_steps;
        for (int j = 0; j < reminder; ++j)
            printf_barry("|");
        
        printf_barry(" done.\n");

    }

    // Resetting the sequence
    for (auto& n: this->nodes)
        n.second.visited = false;

    // So that others now know it was initialized
    initialized = true;

    return;

}

inline void Geese::inherit_support(const Geese & model_, bool delete_support_) {
    
    if (this->model != nullptr)
        throw std::logic_error("There is already a -model- in this Geese. Cannot set a -model- after one is present.");

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

    // This only makes sense (for now) if it is a tip 
    if (!nodes[nodeid].is_leaf())
        return;

    init_node(*nodes[nodeid].parent);

    return;

}

inline void Geese::calc_sequence(Node * n) {

    if (sequence.size() == nodes.size())
        return;

    // First iteration
    if (n == nullptr)
        n = &(nodes.begin()->second);

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

                    includeit[n.ord] = true;
                    reduced_sequence.push_back(i);
                    break;

                }

        } else {

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

inline unsigned int Geese::nterms() const
{

    INITIALIZED()
    return model->nterms() + this->nfuns();

}

inline unsigned int Geese::support_size() const noexcept
{

    if (model == nullptr)
        return 0u;

    return model->support_size();
    
}

inline std::vector< unsigned int > Geese::nannotations() const noexcept
{
    std::vector< unsigned int > ans = {this->n_zeros, this->n_ones};
    return ans;
}

inline std::vector< std::string > Geese::colnames() const
{

    return this->model->colnames();

}

inline unsigned int Geese::parse_polytomies(bool verb) const noexcept {

    unsigned int largest = 0u;
    for (const auto& n : this->nodes)
    {

        unsigned int noff = n.second.noffspring();
        if (noff > 2u)
        {

            if (verb)
                printf_barry("Node id: %i has polytomy size %i\n", n.second.id, noff);
                
        }

        if (noff > largest)
            largest = noff;

    }

    return largest;

}

inline std::vector< std::vector<double> > Geese::observed_counts() {

    // Making room for the output
    std::vector<std::vector<double>> ans;
    ans.reserve(nnodes());

    // Creating counter
    phylocounters::PhyloStatsCounter tmpcount;
    tmpcount.set_counters(this->model->get_counters());

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
    tmpcount.set_counters(this->model->get_counters());

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
        printf_barry("----------\n");
        printf_barry("nodeid: % 3i (%s)\nstate: [", n.second.id, dpl);
        for (uint f = 0u; f < nfuns(); ++f)
            printf_barry("%i, ", (tmparray.D()->states[f] ? 1 : 0));

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
    printf_barry("# of functions           : %i\n", this->nfuns());
    printf_barry("# of nodes [int; leaf]   : [%i; %i]\n", this->nnodes(), this->nleafs());
    printf_barry("# of ann. [zeros; ones]  : [%i; %i]\n", this->n_zeros, this->n_ones);
    printf_barry("# of events [dupl; spec] : [%i; %i]\n", this->n_dupl_events, this->n_spec_events);
    printf_barry("Largest polytomy         : %i\n", parse_polytomies(false));
    printf_barry("\nINFO ABOUT THE SUPPORT\n");
    this->model->print();

}

inline std::mt19937 * Geese::get_rengine()
{
    return this->rengine;
}

inline phylocounters::PhyloCounters * Geese::get_counters()
{
    return this->model->get_counters();
}

inline phylocounters::PhyloModel * Geese::get_model() {
    return this->model;
}

inline phylocounters::PhyloSupport * Geese::get_support() {
    return this->model->get_support();
}

inline std::vector< std::vector< bool > > Geese::get_states() const {
    return this->states;
}

inline std::vector< unsigned int > Geese::get_annotated_nodes() const {

    std::vector< unsigned int > ids(0u);
    for (auto & n : nodes)
    {

        // Counting non-9 annotations
        for (unsigned int f = 0u; f < nfuns(); ++f)
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