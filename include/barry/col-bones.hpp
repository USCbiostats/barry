// #include "typedefs.hpp"
// #include "barray-bones.hpp"

// #ifndef BARRY_COL_BONES_HPP
// #define BARRY_COL_BONES_HPP 1

// template<typename Cell_Type, typename Data_Type>
// class BCol {
// protected:
//     friend class BArray<Cell_Type,Data_Type>;
//     Col_type<Cell_Type> * dat;
//     bool deleted = false;
// public:
//     BCol() : dat(new Col_type<Cell_Type>()) {};
//     BCol(Col_type<Cell_Type> & dat_);
//     BCol(BArray<Cell_Type,Data_Type> & array_, size_t col);
//     ~BCol();

//     std::vector< Cell_Type > as_vector() const;
// };

// template<typename Cell_Type, typename Data_Type>
// inline BCol<Cell_Type,Data_Type>::~BCol() {
//     if (!deleted)
//         delete dat;
// };

// template<typename Cell_Type, typename Data_Type>
// inline BCol<Cell_Type,Data_Type>::BCol(Col_type<Cell_Type> & dat_) {
//     delete = true;
//     dat = &dat_;
// };

// template<typename Cell_Type, typename Data_Type>
// inline BCol<Cell_Type,Data_Type>::BCol(
//     BArray<Cell_Type,Data_Type> & array_, size_t col
//     ) {

//     delete = true;
//     dat = &(array_.get_col(col));

// }

// template<typename Cell_Type, typename Data_Type>
// inline std::vector<Cell_Type> BCol<Cell_Type,Data_Type>::as_vector() const {
//     std::vector<>
// }

// template<typename Cell_Type, typename Data_Type>
// class BCols {
// protected:
//     friend class BArray<Cell_Type,Data_Type>;
//     std::vector< BCol > dat;
// public:
// };

// #endif