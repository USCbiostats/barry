#ifndef BARRY_MODEL_BONES_HPP 
#define BARRY_MODEL_BONES_HPP 1

/**
 * @ingroup stat-models
 * @brief Methods available for approximate (MCMC) sampling
 * @details See `Model::sample_approx()`.
 * - `tnt`: Metropolis-Hastings with the tie/no-tie (TNT) proposal.
 */
enum class SampleApproxMethod {
    tnt
};

/**
 * @ingroup stat-models
 * @brief General framework for discrete exponential models.
 * This class allows generating discrete exponential models in the form of a linear
 * exponential model:
 * \f[
 * \mathrm{Pr}_{\mathcal{A}}(A = a; \theta) = \frac{
 *    \exp{\left(\theta^{\mathbf{t}}c(A)\right)}
 *  }{
 *    \sum_{A'\in \mathcal{A}}\exp{\left(\theta^{\mathbf{t}}c(A')\right)}
 *  }
 * \f]
 * 
 * This implementation aims to reduce the number of times that the support
 * needs to be computed. Models included here use more than a single array, and
 * thus allow the function to recycle support sets as needed. For example,
 * if we are looking at directed graphs all of the same size and without
 * vertex level features, i.e. a model that only counts edges, triangles, etc.
 * then the support needs to be fully computed only once.
 * 
 * @tparam Array_Type Class of `BArray` object.
 * @tparam Data_Counter_Type Any type.
 * @tparam Data_Rule_Type Any type.
 */
template<
    typename Array_Type = BArray<>,
    typename Data_Counter_Type = bool,
    typename Data_Rule_Type  = bool,
    typename Data_Rule_Dyn_Type = bool
    >
class Model {

protected:
    /**
     * @name Random number generation
     * @brief Random number generation
     */
    ///@{
    std::mt19937 * rengine = nullptr;
    bool delete_rengine    = false;

    /**
     * @name Information about the arrays used in the model 
     * @details `stats_target` holds the observed sufficient statistics for each
     * array in the dataset. `array_frequency` contains the frequency with which
     * each of the target stats_target (arrays) shows in the support. `array2support` 
     * maps array indices (0, 1, ...) to the corresponding support.
     * 
     * Each vector of `stats_support` has the data stored in a row-wise order,
     * with each row starting with the weights, e.g., in a model with `k` terms
     * the first k + 1 elements of `stats_support` would be:
     * - weights
     * - term 1
     * - term 2
     * - ...
     * - term k
     */
    ///@{
    std::vector< double >                stats_support;           ///< Sufficient statistics of the model (support)
    std::vector< size_t >                stats_support_sizes;     ///< Number of vectors included in the support.
    std::vector< size_t >                stats_support_sizes_acc; ///< Accumulated number of vectors included in the support.
    std::vector< size_t >                stats_support_n_arrays;  ///< Number of arrays included per support.
    std::vector< std::vector< double > > stats_target;            ///< Target statistics of the model
    std::vector< double >                stats_likelihood;
    std::vector< size_t >                arrays2support;
    ///@}

    /**
      * @brief Map of types of arrays to support sets
      * @details This is of the same length as the vector `stats_target`.
      */
    MapVec_type< double, size_t > keys2support;

    /**
     * @name Container space for the powerset (and its sufficient stats_target)
     * @details This is useful in the case of using simulations or evaluating
     * functions that need to account for the full set of states.
     */
    ///@{
    bool with_pset = false;
    std::vector< std::vector< Array_Type > > pset_arrays; ///< Arrays of the support(s)
    std::vector< double > pset_stats;     ///< Statistics of the support(s)
    std::vector< double > pset_probs;     ///< Probabilities of the support(s)
    std::vector< size_t > pset_sizes;     ///< Number of vectors included in the support.
    std::vector< size_t > pset_locations; ///< Accumulated number of vectors included in the support.
    ///@}
    
    /**
      * @name Functions to compute statistics
      * @details Arguments are recycled to save memory and computation.
      */
    ///@{
    Counters<Array_Type,Data_Counter_Type> *                                counters;
    Rules<Array_Type,Data_Rule_Type> *                                      rules;
    Rules<Array_Type,Data_Rule_Dyn_Type> *                                  rules_dyn;
    Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> support_fun;
    StatsCounter<Array_Type,Data_Counter_Type>                              counter_fun;
    ///@}
    
    /**@brief Vector of the previously used parameters */
    std::vector< std::vector<double> > params_last;
    std::vector< double > normalizing_constants;
    std::vector< bool > first_calc_done;

    bool delete_counters  = false;
    bool delete_rules     = false;
    bool delete_rules_dyn = false;

    /**
     * @brief Transformation of the model
     * 
     * @details When specified, this function will update the model by modifying
     * the linear equation. For example, if the user wanted to add interaction
     * terms, rescale, or apply other operations of the sorts, the user can do such
     * through this function.
     * 
     * The function should return `void` and receive the following arguments:
     * - `data` Pointer to the first element of the set of sufficient statistics
     * - `k` size_t indicating the number of sufficient statistics
     * 
     * @returns
     * Nothing, but it will modify the model data.
     */
    std::function<std::vector<double>(double *, size_t k)>
        transform_model_fun = nullptr;

    std::vector< std::string > transform_model_term_names;
    
    /**
     * @brief Draws a uniform random number in [0, 1), creating the random
     * engine if it does not exist yet.
     */
    double runif_();
public:

    /**
     * @brief Computes the normalizing constant for a given set of parameters
     * @details This function will compute the normalizing constant for a given
     * set of parameters. It will also update the `normalizing_constants` member
     * variable.
    */
    void update_normalizing_constants(
        const std::vector< double > & params,
        BARRY_NCORES_ARG(=1),
        int i = -1
        );

    void update_likelihoods(
        const std::vector< double > & params,
        BARRY_NCORES_ARG(=1)
        );

    void update_pset_probs(
        const std::vector< double > & params,
        BARRY_NCORES_ARG(=1),
        int i = -1
        );
    
    void set_rengine(std::mt19937 * rengine_, bool delete_ = false) {

        if (delete_rengine)
            delete rengine;

        rengine        = rengine_;
        delete_rengine = delete_;
        
    };

    void set_seed(size_t s) {

        if (rengine == nullptr)
        {
            rengine = new std::mt19937;
            delete_rengine = true;
        }

        rengine->seed(s);

    };
    ///@}
        
    Model();
    Model(size_t size_);
    Model(const Model<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> & Model_);
    Model<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> & operator=(
        const Model<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> & Model_
    );

    virtual ~Model() {
        if (delete_counters)
            delete counters;

        if (delete_rules)
            delete rules;

        if (delete_rules_dyn)
            delete rules_dyn;

        if (delete_rengine)
            delete rengine;
    };
    
    void store_psets() noexcept;
    std::vector< double > gen_key(const Array_Type & Array_);
    
    /**
     * @name Wrappers for the `Counters` member. 
     * @details These will add counters to the model, which are shared by the
     * support and the actual counter function.
     */
    ///@{
    void add_counter(Counter<Array_Type, Data_Counter_Type> & counter);
    void add_counter(
        Counter_fun_type<Array_Type,Data_Counter_Type> count_fun_,
        Counter_fun_type<Array_Type,Data_Counter_Type> init_fun_    = nullptr,
        Data_Counter_Type                              data_        = nullptr
    );
    void set_counters(Counters<Array_Type,Data_Counter_Type> * counters_);
    void add_hasher(Hasher_fun_type<Array_Type,Data_Counter_Type> fun_);
    ///@}
    
    /**
     * @name Wrappers for the `Rules` member. 
     * @details These will add rules to the model, which are shared by the
     * support and the actual counter function.
     */
    ///@{
    void add_rule(Rule<Array_Type, Data_Rule_Type> & rule);
    void add_rule(
        Rule_fun_type<Array_Type, Data_Rule_Type> count_fun_,
        Data_Rule_Type                            data_
    );
    
    void set_rules(Rules<Array_Type,Data_Rule_Type> * rules_);

    void add_rule_dyn(Rule<Array_Type, Data_Rule_Dyn_Type> & rule);
    void add_rule_dyn(
        Rule_fun_type<Array_Type, Data_Rule_Dyn_Type> count_fun_,
        Data_Rule_Dyn_Type                            data_
    );
    
    void set_rules_dyn(Rules<Array_Type,Data_Rule_Dyn_Type> * rules_);
    ///@}
    

    /**
     * @brief Adds an array to the support of not already included.
     * @param Array_ array to be added
     * @param force_new If `false`, it will use `keygen` to obtain a double vector
     * and create a hash of it. If the hash has been computed earlier, the support
     * is recycled.
     * 
     * @return The number of the array.
     */
    size_t add_array(const Array_Type & Array_, bool force_new = false);
    
    
    /**
     * @name Likelihood functions.
     * @details Calculation of likelihood functions is done reusing normalizing
     * constants. Before recalculating the normalizing constant, the function 
     * checks whether `params` matches the last set vector of parameters used
     * to compute it.
     * 
     * 
     * @param params Vector of parameters
     * @param as_log When `true`, the function returns the log-likelihood.
     */
    ///@{
    double likelihood(
        const std::vector<double> & params,
        const size_t & i,
        bool as_log = false,
        bool no_update_normconst = false
    );
    
    double likelihood(
        const std::vector<double> & params,
        const Array_Type & Array_,
        int i = -1,
        bool as_log = false,
        bool no_update_normconst = false
    );
    
    double likelihood(
        const std::vector<double> & params,
        const std::vector<double> & target_,
        const size_t & i,
        bool as_log = false,
        bool no_update_normconst = false
    );

    double likelihood(
        const std::vector<double> & params,
        const double * target_,
        const size_t & i,
        bool as_log = false,
        bool no_update_normconst = false
    );
    
    double likelihood_total(
        const std::vector<double> & params,
        bool as_log = false,
        BARRY_NCORES_ARG(=2),
        bool no_update_normconst = false
    );
    ///@}

    /**
     * @name Derivatives of the total log-likelihood
     * @details Exact gradient and Hessian of the total log-likelihood
     * \f[
     * \ell(\theta) = \sum_a \theta^{\mathbf{t}}s(y_a) - \sum_u m_u\log \eta_u(\theta)
     * \f]
     * where \f$s(y_a)\f$ are the observed statistics of array \f$y_a\f$,
     * \f$\eta_u\f$ is the normalizing constant of support \f$u\f$, and
     * \f$m_u\f$ is the number of arrays sharing that support. Then
     * \f[
     * \nabla\ell(\theta) = \sum_a s(y_a) - \sum_u m_u \mathbb{E}_{\mathcal{Y}_u}[s(Y); \theta],
     * \quad
     * \nabla^2\ell(\theta) = -\sum_u m_u \mbox{Cov}_{\mathcal{Y}_u}[s(Y); \theta].
     * \f]
     * These functions do not modify the model (nor its cached normalizing
     * constants), and are computed with a stable (max-shifted) exponent.
     * The derivation and implementation details are in
     * `doc/likelihood-derivatives.md`.
     *
     * @param params Vector of parameters.
     * @param grad Output vector, resized to `nterms()` and filled with the
     * gradient.
     * @return `likelihood_total_and_grad()` returns the log-likelihood,
     * `likelihood_total_grad()` the gradient, and `likelihood_total_hess()`
     * the Hessian as a row-major `nterms() x nterms()` vector.
     */
    ///@{
    double likelihood_total_and_grad(
        const std::vector<double> & params,
        std::vector<double> & grad,
        BARRY_NCORES_ARG(=1)
    ) const;

    std::vector< double > likelihood_total_grad(
        const std::vector<double> & params,
        BARRY_NCORES_ARG(=1)
    ) const;

    std::vector< double > likelihood_total_hess(
        const std::vector<double> & params,
        BARRY_NCORES_ARG(=1)
    ) const;
    ///@}

    /**
     * @name Extract elements by index 
     * @param i Index relative to the array in the model.
     * @param params A new vector of model parameters to compute the normalizing
     * constant.
     * @param as_log When `true` returns the logged version of the normalizing
     * constant.
     */
    ///@{
    const std::vector< double > & get_normalizing_constants() const;
    const std::vector< double > & get_likelihoods() const;

    const std::vector< Array_Type > * get_pset(
        const size_t & i
    );

    const double * get_pset_stats(
        const size_t & i
    );
    ///@}
    
    void print_stats(size_t i) const;

    /**
     * @brief Prints information about the model
     */
    virtual void print() const;
    
    Array_Type sample(const Array_Type & Array_, const std::vector<double> & params = {});
    Array_Type sample(const size_t & i, const std::vector<double> & params);
    
    /**
     * @brief Conditional probability ("Gibbs sampler")
     * 
     * @details Computes the conditional probability of observing
     * \mathrm{Pr}_{\mathcal{Y}}\left(Y(i,j) = 1 \mid Y^C; \theta\right), i.e., the probability of observing the entry Y(i,j) equal
     * to one given the rest of the array.
     * 
     * @param Array_ Array to check
     * @param params Vector of parameters
     * @param i Row entry
     * @param j Column entry
     * @return double The conditional probability
     */
    double conditional_prob(
        const Array_Type & Array_,
        const std::vector< double > & params,
        size_t i,
        size_t j
    );

    /**
     * @name Approximate (MCMC) sampling
     * @brief Sampling arrays without computing the support
     *
     * @details These functions draw arrays from
     * \f[
     * \mathrm{Pr}_{\mathcal{Y}}(Y = y; \theta) =
     *   \frac{\exp\left(\theta^{\mathbf{t}}s(y)\right)}{\eta(\theta)}
     * \f]
     * using Markov chain Monte Carlo, so neither the support nor
     * `store_psets()` is needed. The default (and currently only) method is
     * Metropolis-Hastings with the tie/no-tie (TNT) proposal: with
     * probability 1/2 a random non-zero free cell is proposed for removal;
     * otherwise a random free cell is toggled. A proposal moving from
     * \f$y\f$ to \f$y'\f$ is accepted with probability
     * \f[
     * \min\left\{1,
     *   \exp\left(\theta^{\mathbf{t}}\left[s(y') - s(y)\right]\right)
     *   \frac{q(y' \to y)}{q(y \to y')}\right\},
     * \f]
     * where \f$q\f$ is the TNT proposal probability. Statistics are updated
     * with change statistics only (never recounted), and the chain state
     * persists across calls through `SampleApproxState`. This allows running
     * a burn-in, changing \f$\theta\f$ between calls (as in stochastic
     * approximation), and streaming sufficient statistics (as needed by
     * MCMC-MLE) without storing arrays.
     *
     * Free cells are those allowed by the model rules (`get_rules()`);
     * cells locked by the rules keep their values from the starting array,
     * exactly as in the support used by `sample()` (e.g., the Markov rows
     * of a DEFM). Proposals leading to states that fail the dynamic rules
     * (`get_rules_dyn()`) are rejected. Dynamic rules are evaluated on the
     * statistics of the proposed state, so rules that inspect the array
     * itself are not supported.
     *
     * See `doc/sample-approx.md` for details.
     */
    ///@{
    /**
     * @brief Persistent state of an approximate sampler
     */
    struct SampleApproxState {
        Array_Type            array;          ///< Current state.
        std::vector< double > stats;          ///< Raw sufficient statistics of `array` (before `transform_model_fun`).
        std::vector< size_t > free;           ///< Free cells as (i, j) pairs.
        std::vector< size_t > ties;           ///< Indices (into `free`) of the non-zero free cells.
        std::vector< int >    pos;            ///< Position of each free cell in `ties` (-1 if empty).
        size_t                nsteps    = 0u; ///< Total number of proposals so far.
        size_t                naccepted = 0u; ///< Total number of accepted proposals so far.
    };

    /**
     * @brief Creates the state of an approximate sampler
     * @param Array_ Starting array. Must satisfy the dynamic rules. Cells
     * locked by the rules are kept at their values in `Array_` (the chain
     * conditions on them), so they should hold the intended values, e.g.,
     * an empty diagonal when using `rules_zerodiag`.
     */
    SampleApproxState sample_approx_init(const Array_Type & Array_);

    /**
     * @brief Runs an approximate sampler, updating its state in place
     * @param state State created by `sample_approx_init()`.
     * @param params Vector of parameters (can change between calls).
     * @param n Number of statistics to record.
     * @param thin Number of proposals between recorded statistics. In total,
     * `n * thin` proposals are made.
     * @param stats_out When not `nullptr`, the statistics of the current
     * state (after `transform_model_fun`, i.e., on the scale of
     * `get_stats_target()`) are appended every `thin` proposals, giving `n`
     * rows of `nterms()` values (row-major). Use `nullptr` for burn-in.
     * @param method Sampling method.
     */
    void sample_approx(
        SampleApproxState & state,
        const std::vector< double > & params,
        size_t n                          = 1u,
        size_t thin                       = 1u,
        std::vector< double > * stats_out = nullptr,
        SampleApproxMethod method         = SampleApproxMethod::tnt
    );

    /**
     * @brief Draws a single array with an approximate sampler
     * @param Array_ Starting array. Must satisfy the dynamic rules.
     * @param params Vector of parameters.
     * @param nsteps Number of proposals.
     * @param stats_out When not `nullptr`, the statistics are appended after
     * every proposal (`nsteps` rows of `nterms()` values).
     * @param method Sampling method.
     * @return The array after `nsteps` proposals.
     */
    Array_Type sample_approx(
        const Array_Type & Array_,
        const std::vector< double > & params,
        size_t nsteps                     = 16384u,
        std::vector< double > * stats_out = nullptr,
        SampleApproxMethod method         = SampleApproxMethod::tnt
    );
    ///@}

protected:

    /**
     * @brief One Metropolis-Hastings step with the tie/no-tie proposal
     */
    void sample_approx_tnt_step_(
        SampleApproxState & state,
        const std::vector< double > & params
    );

    /**
     * @brief Applies `transform_model_fun` (if any) to raw statistics
     */
    std::vector< double > transform_stats_(
        const std::vector< double > & stats
    );

public:
    
    /**
     * @name Size of the model
     * 
     * @brief Number of different supports included in the model
     * 
     * This will return the size of `stats_target`.
     * 
     * @return `size()` returns the number of arrays in the model.
     * @return `size_unique()` returns the number of unique arrays (according to
     * the hasher) in the model.
     * @return `nterms()` returns the number of terms in the model.
     */
    ///@{
    size_t size() const noexcept;
    size_t size_unique() const noexcept;
    size_t nterms() const noexcept;
    size_t nrules() const noexcept;
    size_t nrules_dyn() const noexcept;
    size_t support_size() const noexcept;
    std::vector< std::string > colnames() const;
    ///@}

    const std::mt19937 * get_rengine() const;

    Counters<Array_Type,Data_Counter_Type> * get_counters();
    Rules<Array_Type,Data_Rule_Type>       * get_rules();
    Rules<Array_Type,Data_Rule_Dyn_Type>   * get_rules_dyn();
    Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> * get_support_fun();

    /**
     * @brief Raw pointers to the support and target statistics
     * @details 
     * The support of the model is stored as a vector of vector<double>. Each
     * element of it contains the support for an specific type of array included.
     * It represents an array of size `(k + 1) x n unique elements`, with the data
     * stored by-row. The last element of each entry corresponds to the weights,
     * i.e., the frequency with which such sufficient statistics are observed in
     * the support.
     */
    ///@{
    std::vector< std::vector< double > > * get_stats_target();
    std::vector< double > * get_stats_support(); ///< Sufficient statistics of the support(s)
    std::vector< size_t > * get_stats_support_sizes(); ///< Number of vectors included in the support.
    std::vector< size_t > * get_stats_support_sizes_acc(); ///< Accumulated number of vectors included in the support.
    std::vector< size_t > * get_arrays2support();
    std::vector< std::vector< Array_Type > > * get_pset_arrays();
    std::vector< double > * get_pset_stats();  ///< Statistics of the support(s)
    std::vector< double > * get_pset_probs(); 
    std::vector< size_t > * get_pset_sizes();
    std::vector< size_t > * get_pset_locations();
    ///@}

    /**
     * @brief Set the transform_model_fun object
     * @details The transform_model function is used to transform the data
     * 
     * @param data 
     * @param target 
     * @param n_arrays 
     * @param arrays2support 
     */
    ///@{
    void set_transform_model(
        std::function<std::vector<double>(double*,size_t)> fun,
        std::vector< std::string > names
        );
    std::vector<double> transform_model(
        double * data,
        size_t k
    );
    ///@}

};


#endif
