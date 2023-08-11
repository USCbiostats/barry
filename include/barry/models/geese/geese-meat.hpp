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

    }
    
    // Here we have an issue: Some transitions may not be right
    // under the dynamic rules. So not all states can be valid.
    // The arrays and narrays need to be updated once the model
    // is initialized.
    //
    // The later is especially true for leaf nodes, where the
    // limitations are not known until the model is initialized.
    // PhyloStatsCounter stats_counter;
    // stats_counter.set_counters(model->get_counters());
    for (size_t s = 0u; s < states.size(); ++s)
    {

        n.arrays[s] = PhyloArray(n.array, false);

        n.arrays[s].set_data(
            new NodeData(blen, states[s], n.duplication),
            true
        );

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
        map_to_state_id.insert({iter.get_col_vec(0u, false), i});

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
                pset_loc[s][a].push_back(map_to_state_id[tmpstate]);
                
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

inline void Geese::print_nodes() const
{

    printf_barry("GEESE\nINFO ABOUT NODES\n");

    for (const auto & n: nodes)
    {            
        printf_barry("% 4li - Id: %li -- ", n.second.ord, n.second.id);

        // Node type
        printf_barry(
            "node type: %s -- ",
            n.second.is_leaf() ? 
                std::string("leaf").c_str() :
                std::string("internal").c_str()
            );
        
        // Event type
        printf_barry(
            "event type: %s -- ",
            n.second.duplication ?
                std::string("duplication").c_str() :
                std::string("speciation").c_str()
            );

        // Annotations
        printf_barry("ann: [");
        for (const auto & a: n.second.annotations)
        {
            // Print with ']' if last element
            if (&a == &n.second.annotations.back())
            {
                printf_barry("%i] -- ", a);
            }
            else
            {
                printf_barry("%i, ", a);
            }
        }

        // Parent information
        if (n.second.parent == nullptr)
        {
            printf_barry("parent id: (none) -- ");
        } else {
            printf_barry("parent id: %li -- ", n.second.parent->id);
        }

        // Offspring information
        if (n.second.offspring.size() > 0u)
        {
            printf_barry("off ids: [");
            for (const auto & o: n.second.offspring)
            {
                // Same as in previous loop
                if (&o == &n.second.offspring.back())
                {
                    printf_barry("%li].", o->id);
                }
                else
                {
                    printf_barry("%li, ", o->id);
                }
            }
        }

        printf_barry("\n");

    }


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
