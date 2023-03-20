#ifndef DEFM_BONES_HPP
#define DEFM_BONES_HPP 1

class DEFM : public defmcounters::DEFMModel {
private:

    // std::shared_ptr< std::mt19937 > rengine = nullptr;
    // std::shared_ptr< defmcounters::DEFMModel > model = nullptr;

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

    std::vector< std::string > Y_names;
    std::vector< std::string > X_names;
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

    // ~DEFM() {
    //     defmcounters::DEFMModel::~Model();
    // };

    defmcounters::DEFMModel & get_model() {
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

