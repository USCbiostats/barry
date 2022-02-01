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
// https://stackoverflow.com/questions/35055042/difference-between-uint8-t-uint-fast8-t-and-uint-least8-t
typedef unsigned int uint;

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
typedef std::vector< std::pair< std::vector<double>, uint > > Counts_type;

// class Counts_type
// {
// private:
//     std::vector< std::uint_fast32_t > stats_counts;
//     std::vector< double > stats_values;
//     size_t n_stats;
//     unsigned int n_obs;
// public:
//     std::vector< double > operator()
// }

template <class Type_A > class Cell;

template<typename Cell_Type>
using Row_type = Map< uint, Cell<Cell_Type> >;

template<typename Cell_Type>
using Col_type = Map< uint, Cell<Cell_Type>* >;

/**
  * @brief A wrapper class to store `source`, `target`, `val` from a `BArray` object.
  * 
  * @tparam Cell_Type Any type
  */
template<typename Cell_Type>
class Entries {
public:
    std::vector< uint > source;
    std::vector< uint > target;
    std::vector< Cell_Type > val;
    
    Entries() : source(0u), target(0u), val(0u) {};
    Entries(uint n) {
        source.reserve(n);
        target.reserve(n);
        val.reserve(n);
        return;
    };
    
    ~Entries() {};
    
    void resize(uint n) {
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
            for (unsigned int i = 1u; i < dat.size(); ++i)
                hash ^= hasher(dat[i]) + 0x9e3779b9 + (hash<<6) + (hash>>2);
        
        return hash;
        
    }

};

template<typename Ta = double, typename Tb = uint> 
using MapVec_type = std::unordered_map< std::vector< Ta >, Tb, vecHasher<Ta>>;
  

// Mostly relevant in the case of the stats count functions -------------------
template <typename Cell_Type, typename Data_Type> class BArray;
template <typename Array_Type, typename Counter_Type> class Counter;
template <typename Cell_Type, typename Data_Type> class BArrayDense;

/**
 * @brief Counter and rule functions
 * @param Array_Type a BArray
 * @param unit, uint Focal cell
 * @param Data_Type Data associated with the function, for example, id of the attribute
 *  in the Array.
 * @return `Counter_fun_type` a double (the change statistic)
 * @return `Rule_fun_type` a bool. True if the cell is blocked.
 */
///@{
template <typename Array_Type, typename Data_Type>
using Counter_fun_type = std::function<double(const Array_Type &, uint, uint, Data_Type *)>;

template <typename Array_Type, typename Data_Type>
using Rule_fun_type = std::function<bool(const Array_Type &, uint, uint, Data_Type *)>;
///@}

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
        throw std::length_error("-a- and -b- should have the same length.");
    
    unsigned int i = 0;
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
        throw std::length_error("-a- and -b- should have the same length.");
    
    unsigned int i = 0;
    while (static_cast<double>(std::fabs(a[i] - b[i])) < eps) {
        if (++i == a.size())
            return true;
    }
    
    return false;
}
///@}

template <typename T>
inline T vec_inner_prod(
const std::vector< T > & a,
const std::vector< T > & b
) {

    
    if (a.size() != b.size())
        throw std::length_error("-a- and -b- should have the same length.");
    
    double res = 0.0;
    #pragma GCC ivdep
    for (unsigned int i = 0u; i < a.size(); ++i)
        res += (a[i] * b[i]);
    
    return res;

}

template <>
inline double vec_inner_prod(
const std::vector< double > & a,
const std::vector< double > & b
) {
    
    
    if (a.size() != b.size())
        throw std::length_error("-a- and -b- should have the same length.");
    
    double res = 0.0;
    #pragma GCC ivdep
    for (unsigned int i = 0u; i < a.size(); ++i)
        res += (a[i] * b[i]);
    
    return res;

}

#endif
