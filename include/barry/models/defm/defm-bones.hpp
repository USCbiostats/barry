#ifndef DEFM_BONES_HPP
#define DEFM_BONES_HPP 1

class DEFM : public DEFMModel {
private:

    /**
     * @brief Model data
     */
    ///@{
    int * Y           = nullptr; ///< Outcome variable
    int * ID          = nullptr; ///< Individual ids
    double * X        = nullptr; ///< Covariates
    bool column_major = true;    ///< Whether the data is column major or not

    // In case we need a copy of the data
    std::shared_ptr<std::vector< int >> Y_shared;    ///< Outcome variable
    std::shared_ptr<std::vector< int >> ID_shared;   ///< Individual ids
    std::shared_ptr<std::vector< double >> X_shared; ///< Covariates
    
    size_t N;         ///< Number of agents/individuals
    size_t ID_length; ///< Length of the vector IDs
    size_t Y_ncol;    ///< Number of columns in the response
    size_t Y_length;  ///< Length of the vector Y
    size_t X_ncol;    ///< Number of columns in the features
    size_t X_length;  ///< Length of the vector X
    size_t M_order;   ///< Markov order of the model

    std::vector< std::string > Y_names; ///< Names of the response variables
    std::vector< std::string > X_names; ///< Names of the covariates
    std::vector< size_t > start_end;    ///< Start and end of each observation
    std::vector< size_t > model_ord;    ///< Order of the model
    ///@}

public:

    DEFM(
        int * id,
        int * y,
        double * x,
        size_t id_length,
        size_t y_ncol,
        size_t x_ncol,
        size_t m_order,
        bool copy_data = true,
        bool column_major = true
    );

    DEFMModel & get_model() {
        return *this;
    };

    void init();

    double likelihood(std::vector< double > & par, bool as_log = false);
    void simulate(std::vector< double > par, int * y_out);

    size_t get_n_y() const;
    size_t get_n_obs() const;
    size_t get_n_covars() const;
    size_t get_m_order() const;
    size_t get_n_rows() const;

    const int * get_Y() const;
    const int * get_ID() const;
    const double * get_X() const;

    barry::FreqTable<int> motif_census(
        std::vector< size_t > idx
    );

    std::vector< double > logodds(
        const std::vector< double > & par,
        size_t i,
        size_t j
    );

    void set_names(
        std::vector< std::string > Y_names_,
        std::vector< std::string > X_names_
    );

    const std::vector< std::string > & get_Y_names() const;
    const std::vector< std::string > & get_X_names() const;

    void print() const;

    std::vector< bool > is_motif();

};

#endif

