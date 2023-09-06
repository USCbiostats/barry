#ifndef BARRY_TYPEDEFS_HPP
#define BARRY_TYPEDEFS_HPP 1

// Configuration ---------------------------------------------------------------
#include "barry-configuration.hpp"

// Debug
#include "barry-debug.hpp"

// Progress bar
#include "progress.hpp"

// -----------------------------------------------------------------------------

// Basic types
// See this thread
// https://stackoverflow.com/questions/35055042/difference-between-size_t8-t-size_t-fast8-t-and-size_t-least8-t

// Mostly relevant for the BArray definition -----------------------------------

// Constants
/**
  * @brief Integer constants used to specify which cell
  * should be check.
  */
namespace CHECK {
    const int BOTH = -1;
    const int NONE = 0;
    const int ONE  = 1;
    const int TWO  = 2;
}

/**
  * @brief Integer constants used to specify which cell
  * should be check to exist or not.
  */
namespace EXISTS {
    const int BOTH = -1;
    const int NONE = 0;
    const int ONE  = 1;
    const int TWO  = 1;
    
    const int UKNOWN  = -1;
    const int AS_ZERO = 0;
    const int AS_ONE  = 1;
}

/***
  * A single count
  */
typedef std::vector< std::pair< std::vector<double>, size_t > > Counts_type;

// class Counts_type
// {
// private:
//     std::vector< std::size_t_fast32_t > stats_counts;
//     std::vector< double > stats_values;
//     size_t n_stats;
//     size_t n_obs;
// public:
//     std::vector< double > operator()
// }

template <class Type_A > class Cell;

template<typename Cell_Type>
using Row_type = Map< size_t, Cell<Cell_Type> >;

template<typename Cell_Type>
using Col_type = Map< size_t, Cell<Cell_Type>* >;

/**
  * @brief A wrapper class to store `source`, `target`, `val` from a `BArray` object.
  * 
  * @tparam Cell_Type Any type
  */
template<typename Cell_Type>
class Entries {
public:
    std::vector< size_t > source;
    std::vector< size_t > target;
    std::vector< Cell_Type > val;
    
    Entries() : source(0u), target(0u), val(0u) {};
    Entries(size_t n) {
        source.reserve(n);
        target.reserve(n);
        val.reserve(n);
        return;
    };
    
    ~Entries() {};
    
    void resize(size_t n) {
        source.resize(n);
        target.resize(n);
        val.resize(n);
        return;
    }
    
};

// Relevant for anything using vecHasher function ------------------------------
template <typename T>
struct vecHasher
{

    std::size_t operator()(std::vector< T > const&  dat) const noexcept
    {
        
        std::hash< T > hasher;
        std::size_t hash = hasher(dat[0u]);
        
        // ^ makes bitwise XOR
        // 0x9e3779b9 is a 32 bit constant (comes from the golden ratio)
        // << is a shift operator, something like lhs * 2^(rhs)
        if (dat.size() > 1u)
            for (size_t i = 1u; i < dat.size(); ++i)
                hash ^= hasher(dat[i]) + 0x9e3779b9 + (hash<<6) + (hash>>2);
        
        return hash;
        
    }

};

template<typename Ta = double, typename Tb = size_t> 
using MapVec_type = std::unordered_map< std::vector< Ta >, Tb, vecHasher<Ta>>;

/**
 * @brief Ascending sorting an array
 * 
 * It will sort an array solving ties using the next column. Data is
 * stored column-wise.
 * 
 * @tparam T 
 * @param v 
 * @param nrows 
 * @return std::vector<size_t> The sorting index.
 */
inline std::vector< size_t > sort_array(
    const double * v,
    size_t start,
    size_t ncols,
    size_t nrows
    ) {

    // initialize original index locations
    std::vector<size_t> idx(nrows);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),
       [&v,nrows,ncols,start](size_t i1, size_t i2) {

            for (size_t j = 0u; j < ncols; ++j)
            {
                if (*(v + (nrows * j + i1+start)) == *(v + (nrows * j + i2 + start)))
                    continue;   
                else 
                    return *(v + (nrows * j + i1+start)) < *(v + (nrows * j + i2 + start));
            }

            return false;
        });

    return idx;

}   


// Mostly relevant in the case of the stats count functions -------------------
template <typename Cell_Type, typename Data_Type> class BArray;
template <typename Array_Type, typename Counter_Type> class Counter;
template <typename Cell_Type, typename Data_Type> class BArrayDense;

/**
 * @brief Counter and rule functions
 * @param Array_Type a BArray
 * @param unit, size_t Focal cell
 * @param Data_Type Data associated with the function, for example, id of the attribute
 *  in the Array.
 * @return `Counter_fun_type` a double (the change statistic)
 * @return `Rule_fun_type` a bool. True if the cell is blocked.
 */
///@{
template <typename Array_Type, typename Data_Type>
using Counter_fun_type = std::function<double(const Array_Type &, size_t, size_t, Data_Type &)>;

template <typename Array_Type, typename Data_Type>
using Rule_fun_type = std::function<bool(const Array_Type &, size_t, size_t, Data_Type &)>;
///@}

/**
 * @brief Hasher function used by the counter
 * @details Used to characterize the support of the array.
 * 
 * @tparam Array_Type 
 */
template <typename Array_Type, typename Data_Type>
using Hasher_fun_type = std::function<std::vector<double>(const Array_Type &, Data_Type *)>;

// Misc ------------------------------------------------------------------------
/**
 * @brief Compares if -a- and -b- are equal
 * @param a,b Two vectors of the same length
 * @return `true` if all elements are equal.
 */
///@{
template <typename T>
inline bool vec_equal(
    const std::vector< T > & a,
    const std::vector< T > & b
) {
    
    if (a.size() != b.size())
    {
        
        std::string err = "-a- and -b- should have the same length. length(a) = " +
            std::to_string(a.size()) + " and length(b) = " + std::to_string(b.size()) +
            std::string(".");
        throw std::length_error(err);

    }
    
    size_t i = 0;
    while (a[i] == b[i]) {
        if (++i == a.size())
            return true;
    }
    
    return false;
}

template <typename T>
inline bool vec_equal_approx(
    const std::vector< T > & a,
    const std::vector< T > & b,
    double eps = 1e-100
) {
    
    if (a.size() != b.size())
    {
        std::string err = "-a- and -b- should have the same length. length(a) = " +
            std::to_string(a.size()) + " and length(b) = " + std::to_string(b.size()) +
            std::string(".");
        throw std::length_error(err);
    }
    
    size_t i = 0;
    while (static_cast<double>(std::fabs(a[i] - b[i])) < eps) {
        if (++i == a.size())
            return true;
    }
    
    return false;
}
///@}

#ifdef __OPENM
#pragma omp declare simd
#endif
template <typename T>
inline T vec_inner_prod(
    const T * a,
    const T * b,
    size_t n
) {
    
    double res = 0.0;
    #ifdef __OPENM 
    #pragma omp simd reduction(+:res)
    #else
        #ifdef __GNUC__
            #ifndef __clang__
            #pragma GCC ivdep
            #endif
        #endif
    #endif
    for (size_t i = 0u; i < n; ++i)
        res += (*(a + i) * *(b + i));
    
    return res;

}

#ifdef __OPENM
#pragma omp declare simd
#endif
template <>
inline double vec_inner_prod(
    const double * a,
    const double * b,
    size_t n
) {
    
    double res = 0.0;
    #ifdef __OPENMP
    #pragma omp simd reduction(+:res)
    #else
        #ifdef __GNUC__
            #ifndef __clang__
            #pragma GCC ivdep
            #endif
        #endif
    #endif
    for (size_t i = 0u; i < n; ++i)
        res += (*(a + i) * *(b + i));
    
    return res;

}

#endif

