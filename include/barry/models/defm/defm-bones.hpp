#ifndef DEFM_BONES_HPP
#define DEFM_BONES_HPP 1

// #include <vector>
// #include <algorithm>
// #include <random>
// #include <stdexcept>
// #include <memory>

#define DEFM_RANGES(a) \
    size_t __CONCAT(start_,a) = start_end[a * 2u];\
    size_t __CONCAT(end_,a)   = start_end[a * 2u + 1u];\
    size_t __CONCAT(nobs_,a)  = __CONCAT(end_,i) - __CONCAT(start_,i) + 1u;

#define DEFM_LOOP_ARRAYS(a) \
    for (size_t a = 0u; a < (nobs_i - (M_order + 1u) + 1u); ++a)

class DEFM {
private:

    std::shared_ptr< std::mt19937 > rengine = nullptr;
    std::shared_ptr< defmcounters::DEFMModel > model = nullptr;

    /**
     * @brief Model data
     */
    ///@{
    const int * Y = nullptr;    ///< Outcome variable
    const int * ID = nullptr;   ///< Individual ids
    const double * X = nullptr; ///< Covariates
    
    size_t N;             ///< Number of agents/individuals
    size_t ID_length;     ///< Length of the vector IDs
    size_t Y_ncol;        ///< Number of columns in the response
    size_t Y_length;      ///< Length of the vector Y
    size_t X_ncol;        ///< Number of columns in the features
    size_t X_length;      ///< Length of the vector X
    size_t M_order;       ///< Markov order of the model

    std::vector< size_t > start_end;
    std::vector< size_t > model_ord;
    ///@}

public:

    DEFM(
        const int * id,
        const int * y,
        const double * x,
        size_t id_length,
        size_t y_ncol,
        size_t x_ncol,
        size_t m_order
    );

    ~DEFM() {};

    defmcounters::DEFMModel & get_model() {
        return *model;
    };

    void init();

    double likelihood(std::vector< double > & par, bool as_log = false);
    void simulate(std::vector< double > par, int * y_out);

};

inline void DEFM::simulate(
    std::vector< double > par,
    int * y_out
) {

    size_t model_num = 0u; 
    size_t n_entry = M_order * Y_ncol;
    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        DEFM_RANGES(i)
        
        DEFM_LOOP_ARRAYS(proc_n)
        {

            defmcounters::DEFMArray tmp_array = model->sample(model_num++, par);
            for (size_t y = 0u; y < Y_ncol; ++y)
                *(y_out + n_entry++) = tmp_array(0u, y, false);

        }

        n_entry += M_order * Y_ncol;

    }

}

inline DEFM::DEFM(
    const int * id,
    const int * y,
    const double * x,
    size_t id_length,
    size_t y_ncol,
    size_t x_ncol,
    size_t m_order
) {

    // Pointers
    ID = id;
    Y  = y;
    X  = x;

    // Overall dimmensions
    ID_length = id_length;

    Y_ncol    = y_ncol;
    Y_length  = y_ncol * id_length;

    X_ncol    = x_ncol;
    X_length  = x_ncol * id_length;

    M_order   = m_order;

    // Creating the model and engine
    rengine = std::make_shared< std::mt19937 >();
    model   = std::make_shared< defmcounters::DEFMModel >();

    model->set_rengine(&(*(rengine)));

    // Iterating for adding observations
    start_end.reserve(id_length);
    start_end.push_back(0);

    // Identifying the start and end of each observation
    N = 0u;
    for (size_t row = 1u; row < id_length; ++row)
    {

        // Still in the individual
        if (*(id + row) != *(id + row - 1u))
        {

            // End of the previous observation
            start_end.push_back(row - 1u);

            // In the case that the start and end do not fit
            // within the markov process order, then it should fail
            size_t n_rows_i = (row - 1u) - start_end[N++ * 2u] + 1;
            if (n_rows_i < (M_order + 1u))
                throw std::length_error(
                    "Obs. id: " + std::to_string(*(id + row - 1u)) + " (row " +
                    std::to_string(row) + ") has fewer rows (" +
                    std::to_string(n_rows_i) + ") than those needed (" +
                    std::to_string(M_order + 1) + ") for the Markov Model."
                );

            // Beginning of the current
            start_end.push_back(row);

        }
        
    }

    start_end.push_back(id_length - 1u);

    N++;

    return;    

}


inline void DEFM::init() 
{
    // Creating the arrays
    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        size_t start_i = start_end[i * 2u];
        size_t end_i   = start_end[i * 2u + 1u];
        size_t nobs_i  = end_i - start_i + 1u;

        // Creating the observations
        for (size_t n_proc = 0u; n_proc < (nobs_i - (M_order + 1u) + 1u); ++n_proc)
        {

            // Creating the array for process n_proc and setting the data
            defmcounters::DEFMArray array(M_order + 1u, Y_ncol);
            array.set_data(
                new defmcounters::DEFMData(X, (start_i + n_proc), ID_length, X_ncol),
                true // Delete the data
            );

            // Filling-out the array
            for (size_t k = 0u; k < Y_ncol; ++k)
                for (size_t o = 0u; o < (M_order + 1u); ++o)
                    array(o, k) = *(Y + k * ID_length + start_i + n_proc);

            // Adding the rule
            defmcounters::rules_markov_fixed(model->get_rules(), M_order);

            // Adding to the model
            model_ord.push_back( model->add_array(array) );

        }

    }
}

#undef DEFM_RANGES
#undef DEFM_LOOP_ARRAYS

#endif

