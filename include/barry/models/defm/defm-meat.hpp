#ifndef DEFM_MEAT_HPP
#define DEFM_MEAT_HPP 1

inline std::vector< double > keygen_defm(
    const defmcounters::DEFMArray & Array_,
    defmcounters::DEFMCounterData * data
    ) {
    
    size_t nrow = Array_.nrow();
    size_t ncol = Array_.ncol();

    std::vector< double > res(
        2u +                // Rows + cols
        ncol * (nrow - 1u) // Markov cells
        );

    res[0u] = static_cast<double>(nrow);
    res[1u] = static_cast<double>(ncol);

    size_t iter = 2u;
    // Adding the cells
    for (size_t i = 0u; i < (nrow - 1); ++i)
        for (size_t j = 0u; j < ncol; ++j)
            res[iter++] = Array_(i, j);

    return res;

}

#define DEFM_RANGES(a) \
    size_t __CONCAT(start_,a) = start_end[a * 2u];\
    size_t __CONCAT(end_,a)   = start_end[a * 2u + 1u];\
    size_t __CONCAT(nobs_,a)  = __CONCAT(end_,i) - __CONCAT(start_,i) + 1u;

#define DEFM_LOOP_ARRAYS(a) \
    for (size_t a = 0u; a < (nobs_i - M_order); ++a)

inline void DEFM::simulate(
    std::vector< double > par,
    int * y_out
) {

    size_t model_num = 0u; 
    size_t n_entry = M_order * Y_ncol;
    auto idx = model->get_arrays2support();
    defmcounters::DEFMArray last_array;
    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        DEFM_RANGES(i)
        
        DEFM_LOOP_ARRAYS(proc_n)
        {

            // In the first process, we take the data as is
            if (proc_n == 0u)
            {
                last_array = model->sample(idx->at(model_num++), par);
                for (size_t y = 0u; y < Y_ncol; ++y)
                    *(y_out + n_entry++) = last_array(M_order, y, false);

                // last_array.print("i: %li, proc_n: %li\n", i, proc_n);

            }
            else
            // Otherwise, we need to continue using the previous data!
            {
                // Removing the previous row
                defmcounters::DEFMArray tmp_array(M_order + 1u, Y_ncol);
                for (size_t t_i = 1u; t_i < (M_order + 1u); ++t_i)
                    for (size_t t_j = 0u; t_j < Y_ncol; ++t_j)
                        tmp_array(t_i - 1u, t_j) = last_array(t_i, t_j);

                // Setting the data
                tmp_array.set_data(
                    new defmcounters::DEFMData(&tmp_array, X, (start_i + proc_n), X_ncol, ID_length),
                    true // Delete the data
                );

                // Baseline
                // tmp_array.print("baseline i: %li, proc_n: %li\n", i, proc_n);
                // tmp_array.D().print();

                model_num++;
                last_array = model->sample(tmp_array, par);
                for (size_t y = 0u; y < Y_ncol; ++y)
                    *(y_out + n_entry++) = last_array(M_order, y, false);

                // last_array.print("generated i: %li, proc_n: %li\n", i, proc_n);

            }


            
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
    model->store_psets();
    auto kgen = keygen_defm;
    model->add_hasher(kgen);

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

    // Creating the names
    for (auto i = 0u; i < Y_ncol; ++i)
        Y_names.push_back(std::string("y") + std::to_string(i));

    for (auto i = 0u; i < X_ncol; ++i)
        X_names.push_back(std::string("X") + std::to_string(i));

    return;    

}


inline void DEFM::init() 
{

    // Adding the rule
    defmcounters::rules_markov_fixed(model->get_rules(), M_order);

    // Creating the arrays
    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        size_t start_i = start_end[i * 2u];
        size_t end_i   = start_end[i * 2u + 1u];
        size_t nobs_i  = end_i - start_i + 1u;

        // Creating the observations.
        // Number of processes : (N rows) - (Process size)
        for (size_t n_proc = 0u; n_proc < (nobs_i - M_order); ++n_proc)
        {

            // Creating the array for process n_proc and setting the data
            defmcounters::DEFMArray array(M_order + 1u, Y_ncol);
            array.set_data(
                new defmcounters::DEFMData(&array, X, (start_i + n_proc), X_ncol, ID_length),
                true // Delete the data
            );

            // Filling-out the array
            for (size_t k = 0u; k < Y_ncol; ++k)
                for (size_t o = 0u; o < (M_order + 1u); ++o)
                    array(o, k) = *(Y + k * ID_length + start_i + n_proc + o);

            // Adding to the model
            model_ord.push_back( model->add_array(array, true) );

        }

    }

}

inline size_t DEFM::get_n_y() const
{
    return Y_ncol;
}

inline size_t DEFM::get_n_obs() const
{
    return N;
}

inline size_t DEFM::get_n_covars() const
{
    return X_ncol;
}

inline size_t DEFM::get_m_order() const
{
    return M_order;
}

inline size_t DEFM::get_n_rows() const
{
    return ID_length;
}

inline const int * DEFM::get_Y() const
{
    return Y;
}

inline const int * DEFM::get_ID() const
{
    return ID;
}

inline const double * DEFM::get_X() const
{
    return X;
}


inline barry::FreqTable<int> DEFM::motif_census(
        std::vector< size_t > idx
) {

    // Checking all sizes
    for (const auto & i : idx)
        if (i >= Y_ncol)
            throw std::range_error("The -idx- for motif accounting is out of range.");

    barry::FreqTable<int> ans;
    std::vector<int> array(idx.size() * (M_order + 1));

    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        DEFM_RANGES(i)
        
        DEFM_LOOP_ARRAYS(proc_n)
        {

            // Generating an integer array between the parts
            size_t nele = 0u;

            for (size_t o = 0u; o < (M_order + 1u); ++o)
                for (auto & k : idx)
                    array[nele++] = *(Y + k * ID_length + start_i + proc_n + o);

            ans.add(array, nullptr);

        }

    }

    return ans;

}

inline std::vector< double > DEFM::logodds(
    const std::vector< double > & par,
    size_t i_,
    size_t j_
) {
    

    std::vector< double > res(ID_length, std::nan(""));

    for (size_t i = 0u; i < N; ++i)
    {

        // Figuring out how many processes can we observe
        DEFM_RANGES(i)
        
        DEFM_LOOP_ARRAYS(n_proc)
        {

            // Creating the array for process n_proc and setting the data
            defmcounters::DEFMArray array(M_order + 1u, Y_ncol);
            array.set_data(
                new defmcounters::DEFMData(&array, X, (start_i + n_proc), X_ncol, ID_length),
                true // Delete the data
            );

            // Filling-out the array
            for (size_t k = 0u; k < Y_ncol; ++k)
                for (size_t o = 0u; o < (M_order + 1u); ++o)
                    array(o, k) = *(Y + k * ID_length + start_i + n_proc + o);

            double p_1 = model->conditional_prob(array, par, i_, j_);
            res[M_order + start_i + n_proc] = std::log(p_1/(1.0 - p_1));

        }

    }

    return res;


}

inline void DEFM::set_names(
    std::vector< std::string > Y_names_,
    std::vector< std::string > X_names_
) {

    // Checking the length
    if (Y_names_.size() != Y_ncol)
        throw std::length_error("The length of Y_names_ doesn't match the number of dependent variables.");

    if (X_names_.size() != X_ncol)
        throw std::length_error("The length of X_names_ doesn't match the number of dependent variables.");

    Y_names = Y_names_;
    X_names = X_names_;

}

inline const std::vector<std::string > & DEFM::get_Y_names() {
    return Y_names;
}

inline const std::vector<std::string > & DEFM::get_X_names() {
    return X_names;
}

#undef DEFM_RANGES
#undef DEFM_LOOP_ARRAYS

#endif