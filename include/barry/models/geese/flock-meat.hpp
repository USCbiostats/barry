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