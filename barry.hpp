#include <iostream>
#include <cstdarg>
#include <vector>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <map>
#include <algorithm>
#include <utility>
#include <random>
#include <climits>
#include <cfloat>
#include <string>
#include <cstdint>
#include <memory>
#include <regex>
#include <iterator>

#if defined(__OPENMP) || defined(_OPENMP)
#include <omp.h>

// Set the number of threads to match the number of cores
// in the machine

#endif

#ifndef BARRY_HPP
#define BARRY_HPP 

/* Versioning */
#define BARRY_VERSION_MAYOR 0
#define BARRY_VERSION_MINOR 2
#define BARRY_VERSION_PATCH 1
#define BARRY_VERSION BARRY_VERSION_MAYOR ## . ## BARRY_VERSION_MINOR ## . ## BARRY_VERSION_PATCH

static const int barry_version_major = BARRY_VERSION_MAYOR;
static const int barry_version_minor = BARRY_VERSION_MINOR;
static const int barry_version_patch = BARRY_VERSION_PATCH;

/**
  * @brief barry: Your go-to motif accountant
  */
namespace barry {
    
    //! Tree class and TreeIterator class
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/typedefs.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_TYPEDEFS_HPP
#define BARRY_TYPEDEFS_HPP 1

// Configuration ---------------------------------------------------------------
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry//barry-configuration.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_CONFIGURATION_HPP
#define BARRY_CONFIGURATION_HPP

/**
  * @name Configuration MACROS
  * @details These are mostly related to performance. The definitions follow:
  * 
  * - `BARRY_USE_UNORDERED_MAP` If specified, then barry is compiled using
  *   `std::unordered_map`. Otherwise it will use `std::map` for the arrays.
  * 
  * - `BARRY_USE_SAFE_EXP` When specified, it will multiply all likelihoods
  *   in `Model` by (1/-100)/(1/-100) so that numerical overflows are avoided.
  * 
  * - `BARRY_USE_ISFINITE` When specified, it will introduce a macro that 
  *   checks whether the likelihood is finite or not.
  * 
  * - `printf_barry` If not specified, will be defined as `printf`.
  * 
  * - `BARRY_DEBUG_LEVEL`, when defined, will make things verbose.
  */
///@{
#ifdef BARRY_USE_UNORDERED_MAP
    template<typename Ta,typename Tb>
    using Map = std::unordered_map<Ta,Tb>;
#else
    template<typename Ta,typename Tb>
    using Map = std::map<Ta,Tb>;
#endif

#ifdef BARRY_USE_SAFE_EXP
    #define BARRY_SAFE_EXP 
#else
    #define BARRY_SAFE_EXP -100.0
#endif

#ifdef BARRY_USE_ISFINITE
    #define BARRY_ISFINITE(a) if (!std::isfinite( (a) )) \
        throw std::overflow_error("The likelihood function has overflowed.");
#else
    #define BARRY_ISFINITE(a) 
#endif

#ifdef BARRAY_USE_CHECK_SUPPORT
    #define BARRY_CHECK_SUPPORT(x, maxs) if ((x).size() > (maxs)) \
        throw std::length_error("The support has exceeded its maximum size.");
#else
    #define BARRY_CHECK_SUPPORT(x, maxs)
#endif

#ifndef printf_barry
    #define printf_barry printf
#endif

#ifndef BARRY_MAX_NUM_ELEMENTS
    #define BARRY_MAX_NUM_ELEMENTS static_cast< size_t >(std::numeric_limits< size_t >::max() /2u)
#endif

#if defined(__OPENMP) || defined(_OPENMP)
    #define BARRY_WITH_OMP
    #include <omp.h>
#endif


#ifdef BARRY_USE_LATEX
    #define BARRY_WITH_LATEX
#else
    #undef BARRY_WITH_LATEX
#endif

// BARRY_DEBUG_LEVEL: See barry-debug.hpp

// BARRY_PROGRESS_BAR_WIDTH: See progress.hpp

///@}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry//barry-configuration.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



// Debug
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry//barry-debug.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_DEBUG_HPP
#define BARRY_DEBUG_HPP

#ifndef BARRY_DEBUG_LEVEL
    #define BARRY_DEBUG_LEVEL 0
#else
    // The start of the line in every debug print
    #define BARRY_DEBUG_HEADER "[barry]"
    #define BARRY_DEBUG_MSG(a) \
        printf_barry("%s %s\n", BARRY_DEBUG_HEADER, (a));

    // Generic printer (default)
    template <typename T>
    void BARRY_DEBUG_VEC_PRINT(const std::vector<T> & a) {
        printf_barry("%s  [", BARRY_DEBUG_HEADER);
        for(const auto & iter : (a)) {
            printf_barry("%.4f ", static_cast<double>(iter));
        }
        printf_barry("]\n");
        return;
    }

    // Specialization for the printer
    template<>
    inline void BARRY_DEBUG_VEC_PRINT(const std::vector< int > & a) {
        printf_barry("%s  [", BARRY_DEBUG_HEADER);
        for(const auto & iter : (a)) 
        {
            printf_barry("%i ", iter);
        }
        printf_barry("]\n");
        return;
    }

    template<>
    inline void BARRY_DEBUG_VEC_PRINT(const std::vector< std::string > & a) {
        printf_barry("%s \n", BARRY_DEBUG_HEADER);
        for(const auto & iter : (a)) 
        {
            printf_barry("%s %s\n", BARRY_DEBUG_HEADER, iter.c_str());
        }
        printf_barry("%s \n", BARRY_DEBUG_HEADER);
        return;
    }
#endif

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry//barry-debug.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



// Progress bar
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry//progress.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_PROGRESS_HPP
#define BARRY_PROGRESS_HPP

#ifndef BARRY_PROGRESS_BAR_WIDTH
#define BARRY_PROGRESS_BAR_WIDTH 80
#endif

/**
 * @brief A simple progress bar
  */
class Progress {
private:
    int    width;     ///< Total width size (number of bars)
    int    n;         ///< Total number of iterations
    double step_size; ///< Size of the step
    int last_loc;     ///< Last location of the bar
    int cur_loc;      ///< Last location of the bar
    int i;            ///< Current iteration step
    
public:

    Progress(int n_, int width_);
    ~Progress() {};

    void next();
    void end();

};

inline Progress::Progress(int n_, int width_) {


    width     = std::max(7, width_ - 7);
    n         = n_;
    step_size = static_cast<double>(width)/static_cast<double>(n);
    last_loc  = 0;
    i         = 0;

}

inline void Progress::next() {

    cur_loc = std::floor((++i) * step_size);

    for (int j = 0; j < (cur_loc - last_loc); ++j)
    {
        printf_barry("|");
    }

    last_loc = cur_loc;

}

inline void Progress::end() {

    printf_barry(" done.\n");

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry//progress.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



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

#if defined(__OPENMP) || defined(_OPENMP)
#pragma omp declare simd
#endif
template <typename T>
inline T vec_inner_prod(
    const T * a,
    const T * b,
    size_t n
) {
    
    double res = 0.0;
    #if defined(__OPENMP) || defined(_OPENMP) 
    #pragma omp simd reduction(+:res)
    #elif defined(__GNUC__) && !defined(__clang__)
        #pragma GCC ivdep
    #endif
    for (size_t i = 0u; i < n; ++i)
        res += (*(a + i) * *(b + i));
    
    return res;

}

#if defined(__OPENMP) || defined(_OPENMP)
#pragma omp declare simd
#endif
template <>
inline double vec_inner_prod(
    const double * a,
    const double * b,
    size_t n
) {
    
    double res = 0.0;
    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp simd reduction(+:res)
    #elif defined(__GNUC__) && !defined(__clang__)
        #pragma GCC ivdep
    #endif
    for (size_t i = 0u; i < n; ++i)
        res += (*(a + i) * *(b + i));
    
    return res;

}

#endif

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/typedefs.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barry-macros.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_BARRY_MACROS_HPP
#define BARRY_BARRY_MACROS_HPP

#define BARRY_ZERO       Cell<Cell_Type>(0.0)
#define BARRY_ZERO_DENSE static_cast<Cell_Type>(0.0)

#define BARRY_ONE       Cell<Cell_Type>(1.0)
#define BARRY_ONE_DENSE static_cast<Cell_Type>(1.0)

#define BARRY_UNUSED(expr) do { (void)(expr); } while (0);

#if defined(_OPENMP) || defined(__OPENMP)
#define BARRY_NCORES_ARG(default) size_t ncores default
#else 
#define BARRY_NCORES_ARG(default) size_t ncores default
#endif


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barry-macros.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/freqtable.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_STATSDB_HPP 
#define BARRY_STATSDB_HPP 1
  
/**
 * @brief Frequency table of vectors
 * 
 * This is mostly used in `Support`. The main data is contained in the
 * `data` double vector. The matrix is stored in a row-wise fashion, where
 * the first element is the frequency with which the vector is observed.
 * 
 * For example, in a model with `k` terms the first k + 1 elements of
 * `data` would be:
 * 
 * - weights
 * - term 1
 * - term 2
 * - ...
 * - term k
 * 
 */
template<typename T = double> 
class FreqTable {
private:

    std::unordered_map<size_t, size_t> index;
    std::vector< double > data;
    size_t k = 0u;
    size_t n = 0u;

    typename std::unordered_map<size_t, size_t>::iterator iter;
        
public:
    // size_t ncols;
    FreqTable() {};
    ~FreqTable() {};
    
    size_t add(const std::vector< T > & x, size_t * h_precomp);
    
    Counts_type                 as_vector() const;
    const std::vector< double > & get_data() const {return data;};
    const std::unordered_map<size_t,size_t> & get_index() const {return index;};
    
    void clear();
    void reserve(size_t n, size_t k);
    void print() const;

    /**
     * @brief Number of unique elements in the table.
     * (
     * @return size_t 
     */
    size_t size() const noexcept;

    size_t make_hash(const std::vector< T > & x) const;
    
};

template<typename T>  
inline size_t FreqTable<T>::add(
    const std::vector< T > & x,
    size_t * h_precomp
    ) { 
    
    // The term exists, then we add it to the list and we initialize it
    // with a single count
    size_t h;
    if (h_precomp == nullptr)
        h = make_hash(x);
    else
        h = *h_precomp;

    if (k == 0u)
    {

        index.insert({h, 0u});

        data.push_back(1.0);
        data.insert(data.end(), x.begin(), x.end());

        k = x.size();
        n++;

        return h;

    }
    else
    {


        if (x.size() != k)
            throw std::length_error(
                "The value you are trying to add doesn't have the same lenght used in the database."
                );
        
        #if __cplusplus > 201700L
        auto iter2 = index.try_emplace(h, data.size());
        
        if (!iter2.second)
        {
            
            data[(iter2.first)->second] += 1.0;

        }
        else
        {
            data.push_back(1.0);
            data.insert(data.end(), x.begin(), x.end());
            n++;
        }
        #else
        iter = index.find(h);

        if (iter == index.end())
        {
        

            index.insert({h, data.size()});
            data.push_back(1.0);
            data.insert(data.end(), x.begin(), x.end());

            n++;
            
            return h;

        }

        data[(*iter).second] += 1.0;
        
        #endif
        

    }
    
    return h;

}

template<typename T>
inline Counts_type FreqTable<T>::as_vector() const
{ 
    
    Counts_type ans;

    ans.reserve(index.size());

    for (size_t i = 0u; i < n; ++i)
    {
        
        std::vector< double > tmp(k, 0.0);

        for (size_t j = 1u; j < (k + 1u); ++j)
            tmp[j - 1u] = data[i * (k + 1) + j];
        
        ans.push_back(
            std::make_pair<std::vector<double>,size_t>(
                std::move(tmp),
                static_cast<size_t>(data[i * (k + 1u)])
                )
        );

    }
    
    
    return ans;
}

template<typename T>
inline void FreqTable<T>::clear()
{

    index.clear();
    data.clear();

    n = 0u;
    k = 0u;

    return;

}

template<typename T>
inline void FreqTable<T>::reserve(
    size_t n,
    size_t k
)
{

    // Figuring out the max size
    auto nk = std::min(BARRY_MAX_NUM_ELEMENTS, n * k);
    n = nk / k;
    data.reserve(nk);
    index.reserve(n);

    return;

}

// inline void StatsDB::rehash() {
//   stats.rehash();
//   return;
// }

template<typename T>
inline void FreqTable<T>::print() const
{

    size_t grand_total = 0u;

    printf_barry("%7s | %s\n", "Counts", "Stats");

    for (size_t i = 0u; i < n; ++i)
    {

        printf_barry("%7i | ", static_cast<int>(data[i * (k + 1u)]));

        for (size_t j = 1u; j < (k + 1u); ++j)
        {
            printf_barry(" %.2f", data[i * (k + 1) + j]);
        }
        printf_barry("\n");

        grand_total += static_cast<size_t>(data[i * (k + 1u)]);

    }

    printf_barry("Grand total: %i\n", static_cast<int>(grand_total));

    return;

}

template<typename T>
inline size_t FreqTable<T>::size() const noexcept
{

    return index.size();

}

template<typename T>
inline size_t FreqTable<T>::make_hash(const std::vector< T > & x) const
{

    std::hash< T > hasher;
    std::size_t hash = hasher(x[0u]);
    
    // ^ makes bitwise XOR
    // 0x9e3779b9 is a 32 bit constant (comes from the golden ratio)
    // << is a shift operator, something like lhs * 2^(rhs)
    if (x.size() > 1u)
        for (size_t i = 1u; i < x.size(); ++i)
            hash ^= hasher(x[i]) + 0x9e3779b9 + (hash<<6) + (hash>>2);
    
    return hash;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/freqtable.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/cell-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_CELL_BONES_HPP 
#define BARRY_CELL_BONES_HPP 1

/**
 * @brief Entries in BArray.
 * For now, it only has two members:
 * - value: the content
 * - visited: boolean (just a convenient)
 */
template <class Cell_Type > class Cell {
public:
    Cell_Type value;
    bool visited;
    bool active;
    Cell();
    Cell(Cell_Type value_, bool visited_ = false, bool active_ = true) :
        value(value_), visited(visited_), active(active_) {};
    ~Cell() {};
    
    // This is an explicit declaration since in other cases it seems
    // to try to use the move operator, which I do not intent to use.
    Cell(const Cell<Cell_Type>& arg) :
        value(arg.value), visited(arg.visited), active(arg.active) {};
    
    // Copy by assignment
    Cell<Cell_Type>& operator=(const Cell<Cell_Type>& other);
    
    // Move constructor
    Cell(Cell<Cell_Type>&& arg) noexcept:
        value(std::move(arg.value)),
        visited(std::move(arg.visited)),
        active(std::move(arg.active)) {} ;
    
    // Move assign operator
    Cell<Cell_Type>& operator=(Cell<Cell_Type>&& other) noexcept;
    
    void add(Cell_Type x);
    
    // Casting operator (implicit and explicit)
    // int x = Cell<int>(1); // returns 1
    operator Cell_Type() const {return this->value;};

    bool operator==(const Cell<Cell_Type>& rhs ) const;
    bool operator!=(const Cell<Cell_Type>& rhs ) const;
  
};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/cell-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/cell-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_CELL_MEAT_HPP
#define BARRY_CELL_MEAT_HPP 1

template <typename Cell_Type>
Cell<Cell_Type>& Cell<Cell_Type>::operator=(const Cell<Cell_Type>& other) {
    this->value   = other.value;
    this->visited = other.visited;
    this->active  = other.active;
    return *this;
}

template <typename Cell_Type>
Cell<Cell_Type>& Cell<Cell_Type>::operator=(Cell<Cell_Type>&& other) noexcept {
    this->value   = std::move(other.value);
    this->visited = std::move(other.visited);
    this->active  = std::move(other.active);
    return *this;
}

template<typename Cell_Type>
bool Cell<Cell_Type>::operator==(const Cell<Cell_Type>& rhs ) const {

    if (this == *rhs)
        return true;
    
    return this->value == rhs.value;

}

template<typename Cell_Type>
bool Cell<Cell_Type>::operator!=(const Cell<Cell_Type>& rhs ) const {

    return !this->operator==(rhs);
    
}


/***
 * Specializations
 */

template <> inline void Cell<double>::add(double x) {
    value += x;
    return;
}

template <> inline void Cell<size_t>::add(size_t x) {
    value += x;
    return;
}

template <> inline void Cell<int>::add(int x) {
    value += x;
    return;
}

template <> inline void Cell<bool>::add(bool x) {
    value = true;
    return;
}

template<> inline Cell< double >::Cell() : value(1.0), visited(false), active(true) {}
template<> inline Cell< size_t >::Cell() : value(1u), visited(false), active(true) {}
template<> inline Cell< int >::Cell() : value(1), visited(false), active(true) {}
template<> inline Cell< bool >::Cell() : value(true), visited(false), active(true) {}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/cell-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barray-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include <vector>
// #include <unordered_map>
// #include "typedefs.hpp"
// #include "cell-bones.hpp"
// #include "barraycell-bones.hpp"

#ifndef BARRAY_BONES_HPP 
#define BARRAY_BONES_HPP 1

template<typename Cell_Type, typename Data_Type>
class BArrayCell;

template<typename Cell_Type, typename Data_Type>
class BArrayCell_const;

/**
 * @brief Baseline class for binary arrays.
 * 
 * `BArray` class objects are arbitrary arrays
 * in which non-empty cells hold data of type `Cell_Type`. The non-empty cells
 * are stored by row and indexed using `unordered_map`s, i.e.
 * `std::vector< std::unordered_map<size_t,Cell_Type> >`.
 *
 * @tparam Cell_Type Type of cell (any type).
 * @tparam Data_Type Data type of the array (bool default).
 */
template <typename Cell_Type = bool, typename Data_Type = bool>
class BArray {
    friend class BArrayCell<Cell_Type,Data_Type>;
    friend class BArrayCell_const<Cell_Type,Data_Type>;
    // friend class Support<Cell_Type,Data_Type>;
    // friend class StatsCounter<Cell_Type,Data_Type>;
private:
    size_t N;
    size_t M;
    size_t NCells = 0u;
    std::vector< Row_type< Cell_Type > > el_ij;
    std::vector< Col_type< Cell_Type > > el_ji;
    Data_Type * data = nullptr;
    bool delete_data = false;

    static Cell< Cell_Type > Cell_default;
    static const bool dense = false;

public:
    
    /** 
     * This is as a reference, if we need to iterate through the cells and we need
     * to keep track which were visited, we use this as a reference. So that if
     * cell.visited = true and visited = true, it means that we haven't been here
     * yet. Ideally, any routine using this->visited should switch it at the
     * beginning of the routine.
     */
    bool visited = false;
    

    /**
     * @name Constructors
     * 
     * @param N_ Number of rows
     * @param M_ Number of columns
     * @param source An unsigned vector ranging from 0 to N_
     * @param target An size_t vector ranging from 0 to M_
     * @param target When `true` tries to add repeated observations.
     */
    ///@{
    
    /** @brief Zero-size array */
    BArray() : N(0u), M(0u), NCells(0u), el_ij(0u), el_ji(0u) {};
    
    /** @brief Empty array */
    BArray (size_t N_, size_t M_) : N(N_), M(M_), NCells(0u), el_ij(N_), el_ji(M_) {};
    
    /** @brief Edgelist with data */
    BArray (
        size_t N_, size_t M_,
        const std::vector< size_t > & source,
        const std::vector< size_t > & target,
        const std::vector< Cell_Type > & value,
        bool add = true
    );
    
    /** @brief Edgelist with no data (simpler) */
    BArray (
        size_t N_, size_t M_,
        const std::vector< size_t > & source,
        const std::vector< size_t > & target,
        bool add = true
    );
    
    /** @brief Copy constructor */
    BArray(const BArray<Cell_Type,Data_Type> & Array_, bool copy_data = false);
    
    /** @brief Assignment constructor */
    BArray<Cell_Type,Data_Type> & operator=(const BArray<Cell_Type,Data_Type> & Array_);

    /** @brief Move operator */
    BArray(BArray<Cell_Type,Data_Type> && x) noexcept;

    /** @brief Move assignment */
    BArray<Cell_Type,Data_Type> & operator=(BArray<Cell_Type,Data_Type> && x) noexcept;
    ///@}
    
    bool operator==(const BArray<Cell_Type,Data_Type> & Array_);

    ~BArray();
    
    // In principle, copy can be faster by using openmp on the rows
    // since those are independent.
    // BArray(BArray & A);
    
    /**
     * @brief Set the data object
     * 
     * @param data_ 
     * @param delete_data_ 
     */
    ///@{
    void set_data(Data_Type * data_, bool delete_data_ = false);
    Data_Type * D_ptr();
    const Data_Type * D_ptr() const;
    Data_Type & D();
    const Data_Type & D() const;
    void flush_data();
    ///@}
    
    // Function to access the elements
    // bool check_cell
    void out_of_range(size_t i, size_t j) const;
    Cell_Type get_cell(size_t i, size_t j, bool check_bounds = true) const; 
    std::vector< Cell_Type >      get_col_vec(size_t i, bool check_bounds = true) const;
    std::vector< Cell_Type >      get_row_vec(size_t i, bool check_bounds = true) const;
    void                          get_col_vec(std::vector< Cell_Type > * x, size_t i, bool check_bounds = true) const;
    void                          get_row_vec(std::vector< Cell_Type > * x, size_t i, bool check_bounds = true) const;
    const Row_type< Cell_Type > & row(size_t i, bool check_bounds = true) const;
    const Col_type< Cell_Type > & col(size_t i, bool check_bounds = true) const;

    /**
     * @brief Get the edgelist
     * 
     * `Entries` is a class with three objects: Two `std::vector` with the row and
     * column coordinates respectively, and one `std::vector` with the corresponding
     * value of the cell.
     * 
     * @return Entries<Cell_Type> 
     */
    Entries<Cell_Type> get_entries() const;

    /**
     * @name Queries
     * @details `is_empty` queries a single cell. `nrow`, `ncol`, and `nnozero`
     * return the number of rows, columns, and non-zero cells respectively.
     * @param i,j Coordinates
     * @param check_bounds If `false` avoids checking bounds.
     */
    ///@{
    bool is_empty(size_t i, size_t j, bool check_bounds = true) const;
    size_t nrow() const noexcept;
    size_t ncol() const noexcept;
    size_t nnozero() const noexcept;
    Cell<Cell_Type> default_val() const;
    ///@}

    /**
     * @name Cell-wise insertion/deletion
     * @param i,j Row,column
     * @param check_bounds When `true` and out of range, the function throws an
     * error.
     * @param check_exists Wither check if the cell exists (before trying to
     * delete/add), or, in the case of `swap_cells`, check if either of both
     * cells exists/don't exist.
     */
    ///@{  
    BArray<Cell_Type,Data_Type> & operator+=(const std::pair<size_t, size_t> & coords);
    BArray<Cell_Type,Data_Type> & operator-=(const std::pair<size_t, size_t> & coords);
    BArrayCell<Cell_Type,Data_Type> operator()(size_t i, size_t j, bool check_bounds = true);
    const Cell_Type operator()(size_t i, size_t j, bool check_bounds = true) const;
    
    void rm_cell(size_t i, size_t j, bool check_bounds = true, bool check_exists = true);
    
    void insert_cell(size_t i, size_t j, const Cell< Cell_Type > & v, bool check_bounds, bool check_exists);
    void insert_cell(size_t i, size_t j, Cell< Cell_Type > && v, bool check_bounds, bool check_exists);
    void insert_cell(size_t i, size_t j, Cell_Type v, bool check_bounds, bool check_exists);
    
    void swap_cells(
        size_t i0, size_t j0, size_t i1, size_t j1, bool check_bounds = true,
        int check_exists = CHECK::BOTH,
        int * report     = nullptr
        );
    
    void toggle_cell(size_t i, size_t j, bool check_bounds = true, int check_exists = EXISTS::UKNOWN);
    void toggle_lock(size_t i, size_t j, bool check_bounds = true);
    ///@}
    
    /**@name Column/row wise interchange*/
    ///@{
    void swap_rows(size_t i0, size_t i1, bool check_bounds = true);
    void swap_cols(size_t j0, size_t j1, bool check_bounds = true);
    
    void zero_row(size_t i, bool check_bounds = true);
    void zero_col(size_t j, bool check_bounds = true);
    ///@}
    
    void transpose();
    void clear(bool hard = true);
    void resize(size_t N_, size_t M_);
    void reserve();

    // Advances operators
    // void toggle_iterator
    
    // Misc
    void print(const char * fmt = nullptr, ...) const;
    void print_n(size_t nrow, size_t ncol, const char * fmt = nullptr, ...) const;

    /**
     * @name Arithmetic operators
     * 
     */
    ///@{
    BArray<Cell_Type,Data_Type>& operator+=(const BArray<Cell_Type,Data_Type>& rhs);
    BArray<Cell_Type,Data_Type>& operator+=(const Cell_Type & rhs);

    BArray<Cell_Type,Data_Type>& operator-=(const BArray<Cell_Type,Data_Type>& rhs);
    BArray<Cell_Type,Data_Type>& operator-=(const Cell_Type & rhs);
    
    BArray<Cell_Type,Data_Type>& operator/=(const Cell_Type & rhs);
    BArray<Cell_Type,Data_Type>& operator*=(const Cell_Type & rhs);
    ///@}
    
    // /**
    //  * @name Casting between types
    //  */
    // ///@{
    // operator BArray<double,bool>() const;
    // operator BArray<int,bool>() const;
    // operator BArray<size_t,bool>() const;
    // operator BArray<bool,bool>() const;
    // ///@}

    bool is_dense() const noexcept {return dense;};

};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barray-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barraycell-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include "typedefs.hpp"

#ifndef BARRY_BARRAYCELL_BONES_HPP
#define BARRY_BARRAYCELL_BONES_HPP 1

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayCell {
private:
  
    BArray<Cell_Type,Data_Type> * Array;
    size_t i;
    size_t j;
  
public:
  
    BArrayCell(BArray<Cell_Type,Data_Type> * Array_, size_t i_, size_t j_, bool check_bounds = true) : 
    Array(Array_), i(i_), j(j_) {

        if (check_bounds)
        {

            if (i >= Array->nrow())
                throw std::length_error("Row out of range.");
            if (j >= Array->ncol())
                throw std::length_error("Col out of range.");

        }

    };

    ~BArrayCell(){};
    void operator=(const Cell_Type & val);
    void operator+=(const Cell_Type & val);
    void operator-=(const Cell_Type & val);
    void operator*=(const Cell_Type & val);
    void operator/=(const Cell_Type & val);

    operator Cell_Type() const;
    bool operator==(const Cell_Type & val) const;
  
};



template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayCell_const {
private:
    
    const BArray<Cell_Type,Data_Type> * Array;
    size_t i;
    size_t j;
    
public:
  
    BArrayCell_const(const BArray<Cell_Type,Data_Type> * Array_, size_t i_, size_t j_, bool check_bounds = true) : 
    Array(Array_), i(i_), j(j_) {
        if (check_bounds) {

            if (i >= Array->nrow())
                throw std::length_error("Row out of range.");
            if (j >= Array->ncol())
                throw std::length_error("Col out of range.");

        }
    };
    
    ~BArrayCell_const(){};
    
    operator Cell_Type() const;
    bool operator==(const Cell_Type & val) const;
    bool operator!=(const Cell_Type & val) const;
    bool operator<(const Cell_Type & val) const;
    bool operator>(const Cell_Type & val) const;
    bool operator<=(const Cell_Type & val) const;
    bool operator>=(const Cell_Type & val) const;
  
};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barraycell-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barray-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include <stdexcept>
// #include "barray-bones.hpp"

template<typename Cell_Type>
class Cell;

template<typename Cell_Type>
class Cell_const;

#ifndef BARRY_BARRAY_MEAT_HPP
#define BARRY_BARRAY_MEAT_HPP 

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]


template<typename Cell_Type, typename Data_Type>
Cell<Cell_Type> BArray<Cell_Type,Data_Type>::Cell_default = Cell<Cell_Type>(static_cast<Cell_Type>(1.0)); 


// Edgelist with data
template<typename Cell_Type, typename Data_Type> inline BArray<Cell_Type, Data_Type>::BArray (
    size_t N_, size_t M_,
    const std::vector< size_t > & source,
    const std::vector< size_t > & target,
    const std::vector< Cell_Type > & value,
    bool add
) {
  
    if (source.size() != target.size())
        throw std::length_error("-source- and -target- don't match on length.");
    if (source.size() != value.size())
        throw std::length_error("-sorce- and -value- don't match on length.");
    
    // Initializing
    N = N_;
    M = M_;

    el_ij.resize(N);
    el_ji.resize(M);
    
    
    // Writing the data
    for (size_t i = 0u; i < source.size(); ++i) {
      
        // Checking range
        bool empty = this->is_empty(source[i], target[i], true);
        if (add && !empty) {
            ROW(source[i])[target[i]].add(value[i]);
            continue;
        } 
        
        if (!empty)
            throw std::logic_error("The value already exists. Use 'add = true'.");
          
        this->insert_cell(source[i], target[i], value[i], false, false);
    }
    
    return;
  
}

// Edgelist with data
template<typename Cell_Type, typename Data_Type>
inline BArray<Cell_Type, Data_Type>::BArray (
    size_t N_, size_t M_,
    const std::vector< size_t > & source,
    const std::vector< size_t > & target,
    bool add
) {
  
    std::vector< Cell_Type > value(source.size(), (Cell_Type) 1.0);

    if (source.size() != target.size())
      throw std::length_error("-source- and -target- don't match on length.");
    if (source.size() != value.size())
      throw std::length_error("-sorce- and -value- don't match on length.");
    
    // Initializing
    N = N_;
    M = M_;
    
    el_ij.resize(N);
    el_ji.resize(M);
    
    
    // Writing the data
    for (size_t i = 0u; i < source.size(); ++i) {
      
        // Checking range
        if ((source[i] >= N_) || (target[i] >= M_))
            throw std::range_error("Either source or target point to an element outside of the range by (N,M).");
        
        // Checking if it exists
        auto search = ROW(source[i]).find(target[i]);
        if (search != ROW(source[i]).end()) {
            if (!add)
                throw std::logic_error("The value already exists. Use 'add = true'.");
          
            // Increasing the value (this will automatically update the
            // other value)
            ROW(source[i])[target[i]].add(value[i]);
            continue;
        }
        
        // Adding the value and creating a pointer to it
        ROW(source[i]).emplace(
            std::pair<size_t, Cell< Cell_Type> >(
                target[i],
                Cell< Cell_Type >(value[i], visited)
            )
        );
        
        COL(target[i]).emplace(
            source[i],
            &ROW(source[i])[target[i]]
        );

        NCells++;

    }
    
    return;
  
}

template<typename Cell_Type, typename Data_Type>
inline BArray<Cell_Type, Data_Type>::BArray (
    const BArray<Cell_Type,Data_Type> & Array_,
    bool copy_data
) : N(Array_.N), M(Array_.M)
{
  
    // Dimensions
    // el_ij.resize(N);
    // el_ji.resize(M);
    
    std::copy(Array_.el_ij.begin(), Array_.el_ij.end(), std::back_inserter(el_ij));
    std::copy(Array_.el_ji.begin(), Array_.el_ji.end(), std::back_inserter(el_ji));

    // Taking care of the pointers
    for (size_t i = 0u; i < N; ++i)
    {

        for (auto& r: row(i, false))
            COL(r.first)[i] = &ROW(i)[r.first];

    }

    this->NCells  = Array_.NCells;
    this->visited = Array_.visited;
    
    // Data
    if (Array_.data != nullptr)
    {

        if (copy_data)
        {

            data = new Data_Type(* Array_.data );
            delete_data = true;

        } else {

            data = Array_.data;
            delete_data = false;

        }

    }
    
    return;
  
}

template<typename Cell_Type, typename Data_Type>
inline BArray<Cell_Type, Data_Type> &  BArray<Cell_Type, Data_Type>:: operator= (
    const BArray<Cell_Type,Data_Type> & Array_
) {
  
    // Clearing
    if (this != &Array_)
    {
      
        this->clear(true);
        this->resize(Array_.N, Array_.M);
        
        // Entries
        for (size_t i = 0u; i < N; ++i)
        {
          
            if (Array_.nnozero() == nnozero())
                break;
            
            for (auto& r : Array_.row(i, false)) 
                this->insert_cell(i, r.first, r.second.value, false, false);
          
        }
      
        // Data
        if (data != nullptr)
        {

            if (delete_data)
                delete data;
                
            data = nullptr;
            delete_data = false;

        }

        if (Array_.data != nullptr)
        {

            data = new Data_Type(*Array_.data);
            delete_data = true;

        }
      
    }
      
    return *this;
  
}

template<typename Cell_Type, typename Data_Type> inline BArray<Cell_Type, Data_Type>::BArray (
    BArray<Cell_Type, Data_Type> && x
  ) noexcept :
  N(0u), M(0u), NCells(0u),
  data(nullptr),
  delete_data(x.delete_data)
  {

    this->clear(true);
    this->resize(x.N, x.M);
    
    // Entries
    for (size_t i = 0u; i < N; ++i) {
      
        if (x.nnozero() == nnozero())
            break;
        
        for (auto& r : x.row(i, false)) 
            this->insert_cell(i, r.first, r.second.value, false, false);
      
    }

    // Managing data
    if (x.data != nullptr)
    {

        if (x.delete_data)
        {

            data = new Data_Type(*x.data);
            delete_data = true;

        } else {
            data = x.data;
            delete_data = false;
        }


    }

}

template<typename Cell_Type, typename Data_Type> inline BArray<Cell_Type, Data_Type> &  BArray<Cell_Type, Data_Type>:: operator= (
    BArray<Cell_Type, Data_Type> && x
) noexcept {
  
    // Clearing
    if (this != &x) {
      
        this->clear(true);
        this->resize(x.N, x.M);
        
        // Entries
        for (size_t i = 0u; i < N; ++i) {
          
            if (x.nnozero() == nnozero())
                break;
            
            for (auto& r : x.row(i, false)) 
                this->insert_cell(i, r.first, r.second.value, false, false);
          
        }
      
        // Data
        if (data != nullptr)
        {

            if (delete_data)
                delete data;
            data = nullptr;
            delete_data = false;

        }

        if (x.data != nullptr)
        {

            data = new Data_Type( *x.data );
            delete_data = true;

        }

        // x.data = nullptr;
        // x.delete_data = false;
      
    }
      
    return *this;
  
}

template<typename Cell_Type, typename Data_Type> inline bool  BArray<Cell_Type, Data_Type>:: operator== (
    const BArray<Cell_Type, Data_Type> & Array_
) {
    
    // Dimension and number of cells used
    if ((N != Array_.nrow()) | (M != Array_.ncol()) | (NCells != Array_.nnozero()))
        return false;
    
    // One holds, and the other doesn't.
    if ((!data & Array_.data) | (data & !Array_.data))
        return false;
    
    if (this->el_ij != Array_.el_ij)
        return false;
    
    return true;
}

template<typename Cell_Type, typename Data_Type> inline BArray<Cell_Type, Data_Type>::~BArray () {
    
    if (delete_data && (data != nullptr))
        delete data;
    
    return;
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: set_data (
    Data_Type * data_, bool delete_data_
) {  

    if ((data != nullptr) && delete_data)
        delete data;
    
    data        = data_;
    delete_data = delete_data_;
    
    return;
    
}

template<typename Cell_Type, typename Data_Type> inline Data_Type *  BArray<Cell_Type, Data_Type>:: D_ptr ()
{
    return this->data;
}

template<typename Cell_Type, typename Data_Type>
inline const Data_Type * BArray<Cell_Type,Data_Type>::D_ptr() const
{
    return this->data;
}

template<typename Cell_Type, typename Data_Type> inline Data_Type &  BArray<Cell_Type, Data_Type>:: D ()
{
    return *this->data;
}

template<typename Cell_Type, typename Data_Type>
inline const Data_Type & BArray<Cell_Type,Data_Type>::D() const
{
    return *this->data;
}

template<typename Cell_Type, typename Data_Type>
inline void BArray<Cell_Type,Data_Type>::flush_data()
{

    if (delete_data)
    {
        delete data;
        delete_data = false;
    }

    data = nullptr;

    return;

}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: out_of_range (
    size_t i,
    size_t j
) const {

    if (i >= N)
        throw std::range_error("The row is out of range.");
    else if (j >= M)
        throw std::range_error("The column is out of range.");
    return;

}
    
template<typename Cell_Type, typename Data_Type> inline Cell_Type  BArray<Cell_Type, Data_Type>:: get_cell (
    size_t i,
    size_t j,
    bool check_bounds
) const {
    
    // Checking boundaries  
    if (check_bounds)
        out_of_range(i,j);
    
    if (ROW(i).size() == 0u)
        return (Cell_Type) 0.0;
    
    // If it is not empty, then find and return
    auto search = ROW(i).find(j);
    if (search != ROW(i).end())
        return search->second.value;
    
    // This is if it is empty
    return (Cell_Type) 0.0;
    
}

template<typename Cell_Type, typename Data_Type> inline std::vector< Cell_Type >  BArray<Cell_Type, Data_Type>:: get_row_vec (
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    std::vector< Cell_Type > ans(ncol(), (Cell_Type) false);
    for (const auto & iter : row(i, false)) 
        ans[iter.first] = iter.second.value; //this->get_cell(i, iter->first, false);
    

    return ans;
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: get_row_vec (
    std::vector< Cell_Type > * x,
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    for (const auto & iter : row(i, false)) 
        x->at(iter.first) = iter.second.value; // this->get_cell(i, iter->first, false);
    
}

template<typename Cell_Type, typename Data_Type> inline std::vector< Cell_Type >  BArray<Cell_Type, Data_Type>:: get_col_vec (
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    std::vector< Cell_Type > ans(nrow(), (Cell_Type) false);
    for (const auto iter : col(i, false)) 
        ans[iter.first] = iter.second->value;//this->get_cell(iter->first, i, false);
    
    return ans;

}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: get_col_vec (
    std::vector<Cell_Type> * x,
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    for (const auto & iter : col(i, false)) 
        x->at(iter.first) = iter.second->value;//this->get_cell(iter->first, i, false);
    
}

template<typename Cell_Type, typename Data_Type> inline const Row_type< Cell_Type > &  BArray<Cell_Type, Data_Type>:: row (
    size_t i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, 0u);

    return this->el_ij[i];

}

template<typename Cell_Type, typename Data_Type> inline const Col_type< Cell_Type > &  BArray<Cell_Type, Data_Type>:: col (
    size_t i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(0u, i);

    return this->el_ji[i];
    
}

template<typename Cell_Type, typename Data_Type> inline Entries< Cell_Type >  BArray<Cell_Type, Data_Type>:: get_entries () const {
    
    Entries<Cell_Type> res(NCells);
    
    for (size_t i = 0u; i < N; ++i) {
        
        if (ROW(i).size() == 0u)
            continue;
        
        for (auto col = ROW(i).begin(); col != ROW(i).end(); ++col) {
            res.source.push_back(i),
            res.target.push_back(col->first),
            res.val.push_back(col->second.value);
        }
    }
    
    return res;
}

template<typename Cell_Type, typename Data_Type> inline bool  BArray<Cell_Type, Data_Type>:: is_empty (
    size_t i,
    size_t j,
    bool check_bounds
) const {
    
    if (check_bounds)
        out_of_range(i, j);
    
    if (ROW(i).size() == 0u)
        return true;
    else if (COL(j).size() == 0u)
        return true;
    
    if (ROW(i).find(j) == ROW(i).end())
        return true;
    
    return false;
    
}


template<typename Cell_Type, typename Data_Type> inline size_t  BArray<Cell_Type, Data_Type>:: nrow () const noexcept {
    return N;
}


template<typename Cell_Type, typename Data_Type> inline size_t  BArray<Cell_Type, Data_Type>:: ncol () const noexcept {
    return M;
}


template<typename Cell_Type, typename Data_Type> inline size_t  BArray<Cell_Type, Data_Type>:: nnozero () const noexcept {
    return NCells;
}

template<typename Cell_Type, typename Data_Type> inline Cell< Cell_Type >  BArray<Cell_Type, Data_Type>:: default_val () const {
    return this->Cell_default;
}

template<typename Cell_Type, typename Data_Type> inline BArray<Cell_Type, Data_Type> &  BArray<Cell_Type, Data_Type>:: operator+= (
    const std::pair<size_t,size_t> & coords
) {
    
    this->insert_cell(
            coords.first,
            coords.second,
            this->Cell_default,
            true, true
    );
    
    return *this;
    
}

template<typename Cell_Type, typename Data_Type> inline BArray<Cell_Type, Data_Type> &  BArray<Cell_Type, Data_Type>:: operator-= (
    const std::pair<size_t,size_t> & coords
) {
    
    this->rm_cell(
            coords.first,
            coords.second,
            true, true
    );
    
    return *this;
    
}

template<typename Cell_Type, typename Data_Type>
inline BArrayCell<Cell_Type,Data_Type> BArray<Cell_Type, Data_Type>::operator()(  
    size_t i,
    size_t j,
    bool check_bounds
) {
    
    return BArrayCell<Cell_Type,Data_Type>(this, i, j, check_bounds);
    
}

template<typename Cell_Type, typename Data_Type>
inline const Cell_Type BArray<Cell_Type, Data_Type>::operator() (  
    size_t i,
    size_t j,
    bool check_bounds
) const {
    
    return get_cell(i, j, check_bounds);
    
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: rm_cell (
    size_t i,
    size_t j,
    bool check_bounds,
    bool check_exists
) {
    
    // Checking the boundaries
    if (check_bounds)
        out_of_range(i,j);
    
    if (check_exists) {
        // Nothing to do
        if (ROW(i).size() == 0u)
            return;
        
        // Checking the counter part
        if (COL(j).size() == 0u)
            return;
        
        // Hard work, need to remove it from both, if it exist
        if (ROW(i).find(j) == ROW(i).end())
            return;
    }
    
    // Remove the pointer first (so it wont point to empty)
    COL(j).erase(i);
    ROW(i).erase(j);
    
    NCells--;
    
    return;
    
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: insert_cell (
        size_t i,
        size_t j,
        const Cell< Cell_Type> & v,
        bool check_bounds,
        bool check_exists
    ) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    if (check_exists) {
        
        // Checking if nothing here, then we move along
        if (ROW(i).size() == 0u) {
            
            ROW(i).insert(std::pair< size_t, Cell<Cell_Type>>(j, v));
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            return;
            
        }
        
        // In this case, the row exists, but we are checking that the value is empty  
        if (ROW(i).find(j) == ROW(i).end()) {
            
            ROW(i).insert(std::pair< size_t, Cell<Cell_Type>>(j, v)); 
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            
        } else {
            throw std::logic_error("The cell already exists.");
        }
        
        
    } else {
        
        ROW(i).insert(std::pair< size_t, Cell<Cell_Type>>(j, v));
        COL(j).emplace(i, &ROW(i)[j]);
        NCells++;
        
    }
    
    return;
    
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: insert_cell (
        size_t i,
        size_t j,
        Cell< Cell_Type> && v,
        bool check_bounds,
        bool check_exists
    ) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    if (check_exists) {
        
        // Checking if nothing here, then we move along
        if (ROW(i).size() == 0u) {
            
            ROW(i).insert(std::pair< size_t, Cell<Cell_Type>>(j, v));
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            return;
            
        }
        
        // In this case, the row exists, but we are checking that the value is empty  
        if (ROW(i).find(j) == ROW(i).end()) {
            
            ROW(i).insert(std::pair< size_t, Cell<Cell_Type>>(j, v)); 
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            
        } else {
            throw std::logic_error("The cell already exists.");
        }
        
        
    } else {
        
        ROW(i).insert(std::pair< size_t, Cell<Cell_Type>>(j, v));
        COL(j).emplace(i, &ROW(i)[j]);
        NCells++;
        
    }
    
    return;
    
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: insert_cell (
    size_t i,
    size_t j,
    Cell_Type v,
    bool check_bounds,
    bool check_exists
) {
        
    return insert_cell(i, j, Cell<Cell_Type>(v, visited), check_bounds, check_exists);

}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: swap_cells (
    size_t i0, size_t j0,
    size_t i1, size_t j1,
    bool check_bounds,
    int check_exists,
    int * report
) {
    
    if (check_bounds) {
        out_of_range(i0,j0);
        out_of_range(i1,j1);
    }
    
    // Simplest case, we know both exists, so we don't need to check anything
    if (check_exists == CHECK::NONE)
    {
        
        // Just in case, if this was passed
        if (report != nullptr)
            (*report) = EXISTS::BOTH;
        
        // If source and target coincide, we do nothing
        if ((i0 == i1) && (j0 == j1)) 
            return;
        
        // Using the initializing by move, after this, the cell becomes
        // invalid. We use pointers instead as this way we access the Heap memory,
        // which should be faster to access.
        Cell<Cell_Type> c0(std::move(ROW(i0)[j0]));
        rm_cell(i0, j0, false, false);
        Cell<Cell_Type> c1(std::move(ROW(i1)[j1]));
        rm_cell(i1, j1, false, false);
        
        // Inserting the cells by reference, these will be deleted afterwards
        insert_cell(i0, j0, c1, false, false);
        insert_cell(i1, j1, c0, false, false);
        
        return;
        
    }
    
    bool check0, check1;
    if (check_exists == CHECK::BOTH)
    {
        
        check0 = !is_empty(i0, j0, false);
        check1 = !is_empty(i1, j1, false);
        
    } else if (check_exists == CHECK::ONE) {
        
        check0 = !is_empty(i0, j0, false);
        check1 = true;
        
    } else if (check_exists == CHECK::TWO) {
        
        check0 = true;
        check1 = !is_empty(i1, j1, false);
        
    }
    
    if (report != nullptr) 
        (*report) = EXISTS::NONE;
    
    // If both cells exists
    if (check0 & check1)
    {
        
        if (report != nullptr) 
            (*report) = EXISTS::BOTH;
        
        // If source and target coincide, we do nothing
        if ((i0 == i1) && (j0 == j1)) 
            return;
        
        Cell<Cell_Type> c0(std::move(ROW(i0)[j0]));
        rm_cell(i0, j0, false, false);
        Cell<Cell_Type> c1(std::move(ROW(i1)[j1]));
        rm_cell(i1, j1, false, false);
        
        insert_cell(i0, j0, c1, false, false);
        insert_cell(i1, j1, c0, false, false);
        
    } else if (!check0 & check1) { // If only the second exists
        
        if (report != nullptr) 
            (*report) = EXISTS::TWO;
        
        insert_cell(i0, j0, ROW(i1)[j1], false, false);
        rm_cell(i1, j1, false, false);
        
    } else if (check0 & !check1) {
        
        if (report != nullptr) 
            (*report) = EXISTS::ONE;
        
        insert_cell(i1, j1, ROW(i0)[j0], false, false);
        rm_cell(i0, j0, false, false);
        
    }
    
    return;
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: toggle_cell (
    size_t i,
    size_t j,
    bool check_bounds,
    int check_exists
) {
    
    if (check_bounds)
        out_of_range(i, j);
    
    if (check_exists == EXISTS::UKNOWN) {
        
        if (is_empty(i, j, false)) {
            insert_cell(i, j, BArray<Cell_Type, Data_Type>::Cell_default, false, false);
            ROW(i)[j].visited = visited;
        } else
            rm_cell(i, j, false, false);
        
    } else if (check_exists == EXISTS::AS_ONE) {
        
        rm_cell(i, j, false, false);
        
    } else if (check_exists == EXISTS::AS_ZERO) {
        
        insert_cell(i, j, BArray<Cell_Type,Data_Type>::Cell_default, false, false);
        ROW(i)[j].visited = visited;
        
    }
    
    return;
    
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: swap_rows (
    size_t i0,
    size_t i1,
    bool check_bounds
) {
  
    if (check_bounds) {
        out_of_range(i0,0u);
        out_of_range(i1,0u);
    }
    
    bool move0=true, move1=true;
    if (ROW(i0).size() == 0u) move0 = false;
    if (ROW(i1).size() == 0u) move1 = false;
    
    if (!move0 && !move1)
        return;
    
    // Swapping happens naturally, need to take care of the pointers
    // though
    ROW(i0).swap(ROW(i1));
    
    // Delete the thing
    if (move0)
        for (auto& i: row(i1, false))
            COL(i.first).erase(i0);
    
    if (move1)
        for (auto& i: row(i0, false))
            COL(i.first).erase(i1);
    
    // Now, point to the thing, if it has something to point at. Recall that
    // the indices swapped.
    if (move1)
        for (auto& i: row(i0, false))
            COL(i.first)[i0] = &ROW(i0)[i.first];
    
    if (move0)
        for (auto& i: row(i1, false))
            COL(i.first)[i1] = &ROW(i1)[i.first];
    
    return;

}

// This swapping is more expensive overall
template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: swap_cols (
    size_t j0,
    size_t j1,
    bool check_bounds
) {
  
    if (check_bounds) {
        out_of_range(0u, j0);
        out_of_range(0u, j1);
    }
    
    // Which ones need to be checked
    bool check0 = true, check1 = true;
    if (COL(j0).size() == 0u) check0 = false;
    if (COL(j1).size() == 0u) check1 = false;
    
    if (check0 && check1) {
      
        // Just swapping one at a time
        int status;
        Col_type<Cell_Type> col_tmp = COL(j1);
        Col_type<Cell_Type> col1 = COL(j0);
        for (auto iter = col1.begin(); iter != col1.end(); ++iter) {
            
            // Swapping values (col-wise)
            swap_cells(iter->first, j0, iter->first, j1, false, CHECK::TWO, &status);
            
            // Need to remove it, so we don't swap that as well
            if (status == EXISTS::BOTH)
                col_tmp.erase(iter->first);
        }
        
        // If there's anything left to move, we start moving it, otherwise, we just
        // skip it
        if (col_tmp.size() != 0u) {
          
            for (auto iter = col_tmp.begin(); iter != col_tmp.end(); ++iter) {
                insert_cell(iter->first, j0, *iter->second, false, false);
                rm_cell(iter->first, j1);
            }
          
        }
      
    } else if (check0 && !check1) {
      
        // 1 is empty, so we just add new cells and remove the other ones
        for (auto iter = COL(j0).begin(); iter != COL(j0).begin(); ++iter)
            insert_cell(iter->first, j1, *iter->second, false, false);
        
        // Setting the column to be zero
        COL(j0).empty();
      
    } else if (!check0 && check1) {
      
        // 1 is empty, so we just add new cells and remove the other ones
        for (auto iter = COL(j1).begin(); iter != COL(j1).begin(); ++iter) {
          
            // Swapping values (col-wise)
            insert_cell(iter->first, j0, *iter->second, false, false);

        }
        
        // Setting the column to be zero
        COL(j1).empty();
      
    }
    
    
    return;
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: zero_row (
    size_t i,
    bool check_bounds
) {
  
    if (check_bounds)
        out_of_range(i, 0u);
    
    // Nothing to do
    if (ROW(i).size() == 0u)
        return;
    
    // Else, remove all elements
    auto row0 = ROW(i);
    for (auto row = row0.begin(); row != row0.end(); ++row) 
        rm_cell(i, row->first, false, false);
    
    return;
  
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: zero_col (
    size_t j,
    bool check_bounds
) {
  
    if (check_bounds)
        out_of_range(0u, j);
    
    // Nothing to do
    if (COL(j).size() == 0u)
        return;
    
    // Else, remove all elements
    auto col0 = COL(j);
    for (auto col = col0.begin(); col != col0.end(); ++col) 
        rm_cell(col->first, j, false, false);
    
    return;
  
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: transpose () {
  
    // Start by flipping the switch 
    visited = !visited;
    
    // Do we need to resize (increase) either?
    if      (N > M) el_ji.resize(N);
    else if (N < M) el_ij.resize(M);
    
    // size_t N0 = N, M0 = M;
    int status;
    for (size_t i = 0u; i < N; ++i)
    {
        
        // Do we need to move anything?
        if (ROW(i).size() == 0u)
            continue;
        
        // We now iterate changing rows
        Row_type<Cell_Type> row = ROW(i);
        for (auto col = row.begin(); col != row.end(); ++col)
        {
          
            // Skip if in the diagoal
            if (i == col->first)
            {
                ROW(i)[i].visited = visited;
                continue;
            }
            
            // We have not visited this yet, we need to change that
            if (ROW(i)[col->first].visited != visited)
            {
                
                // First, swap the contents
                swap_cells(i, col->first, col->first, i, false, CHECK::TWO, &status);
                
                // Changing the switch
                if (status == EXISTS::BOTH)
                    ROW(i)[col->first].visited = visited;
                
                ROW(col->first)[i].visited = visited;
              
            }
          
        }
      
    }
    
    // Shreding. Note that no information should have been lost since, hence, no
    // change in NCells.
    if (N > M) el_ij.resize(M);
    else if (N < M) el_ji.resize(N);
    
    // Swapping the values
    std::swap(N, M);
    
    return;

}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: clear (
    bool hard
) {
    
    if (hard)
    {
      
        el_ji.clear();
        el_ij.clear();
        
        el_ij.resize(N);
        el_ji.resize(M);
        NCells = 0u;
      
    } else {
        
        for (size_t i = 0u; i < N; ++i)
            zero_row(i, false);
        
    }
    
    return;
    
}

template<typename Cell_Type, typename Data_Type> inline void  BArray<Cell_Type, Data_Type>:: resize (
    size_t N_,
    size_t M_
) {
  
    // Removing rows
    if (N_ < N)
        for (size_t i = N_; i < N; ++i)
            zero_row(i, false);
    
    // Removing cols
    if (M_ < M)
        for (size_t j = M_; j < M; ++j)
            zero_col(j, false);
    
    // Resizing will invalidate pointers and values out of range
    if (M_ != M) {
        el_ji.resize(M_);
        M = M_;
    }
    
    if (N_ != N) {
        el_ij.resize(N_);
        N = N_;
    }
    
    
    return;

}

template<typename Cell_Type, typename Data_Type>
inline void  BArray<Cell_Type, Data_Type>:: reserve () {
#ifdef BARRAY_USE_UNORDERED_MAP
    for (size_t i = 0u; i < N; i++)
        ROW(i).reserve(M);
    
    for (size_t i = 0u; i < M; i++)
        COL(i).reserve(N);
#endif
    return;
  
}

template<typename Cell_Type, typename Data_Type>
inline void  BArray<Cell_Type, Data_Type>:: print (
    const char * fmt,
    ...
) const {
  

    std::va_list args;
    va_start(args, fmt);
    print_n(N, M, fmt, args);
    va_end(args);    
    
    return;

}

template<typename Cell_Type, typename Data_Type>
inline void  BArray<Cell_Type, Data_Type>:: print_n (
    size_t nrow,
    size_t ncol,
    const char * fmt,
    ...
) const {

    if (nrow > N)
        nrow = N;

    if (ncol > M)
        ncol = M;

    std::va_list args;
    va_start(args, fmt);
    if (fmt != nullptr)
        printf_barry(fmt, args);
    va_end(args);

    for (size_t i = 0u; i < nrow; ++i)
    {

        #ifdef BARRY_DEBUG_LEVEL
            #if BARRY_DEBUG_LEVEL > 1
                printf_barry("%s [%3i,]", BARRY_DEBUG_HEADER, i);
            #endif
        #else
        printf_barry("[%3i,] ", i);
        #endif
        for (size_t j = 0u; j < ncol; ++j) {
            if (this->is_empty(i, j, false))
            {
                printf_barry("    . ");
            } else {
                printf_barry(" %.2f ", static_cast<double>(this->get_cell(i, j, false)));
            }
            
        }

        printf_barry("\n");

    }

    if (nrow < N) {
        printf_barry("Skipping %i rows. ", static_cast< int >(N - nrow));
    }
    
    if (ncol < M) {
        printf_barry("Skipping %i columns. ", static_cast< int >(M - ncol));
    }

    if (nrow < N || ncol < M) {
        printf_barry("\n");
    }
    
    
    return;

}

#undef ROW
#undef COL

#endif

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barray-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barraycell-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include "barraycell-bones.hpp"

#ifndef BARRY_BARRAYCELL_MEAT_HPP
#define BARRY_BARRAYCELL_MEAT_HPP 1

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {
    
    if (Array->is_empty(i, j, false)) {
        Array->insert_cell(i, j, val, false, false);
    } else {
        Array->el_ij.at(i).at(j).value = val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
    
    if (Array->is_empty(i, j, false)) {
        Array->insert_cell(i, j, val, false, false);
    } else {
        Array->el_ij.at(i).at(j).value += val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
    
    if (Array->is_empty(i, j, false)) {
        Array->insert_cell(i, j, -val, false, false);
    } else {
        Array->el_ij.at(i).at(j).value -= val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
    
    if (!Array->is_empty(i, j, false)) {
        Array->el_ij.at(i).at(j).value *= val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayCell<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
    
    if (!Array->is_empty(i, j, false)) {
        Array->el_ij.at(i).at(j).value /= val;
    }

}

template<typename Cell_Type,typename Data_Type>
inline BArrayCell<Cell_Type,Data_Type>::operator Cell_Type() const {
        return Array->get_cell(i, j, false);
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) == static_cast<Cell_Type>(val);  
}

template<typename Cell_Type,typename Data_Type>
inline BArrayCell_const<Cell_Type,Data_Type>::operator Cell_Type() const {
        return Array->get_cell(i, j, false);
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) == static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator!=(const Cell_Type & val) const {
    return !(this->operator==(val));
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator<(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) < static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator>(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) > static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator<=(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) <= static_cast<Cell_Type>(val);    
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayCell_const<Cell_Type,Data_Type>::operator>=(const Cell_Type & val) const {
    return Array->get_cell(i, j, false) >= static_cast<Cell_Type>(val);    
}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barraycell-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barray-meat-operators.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include <stdexcept>
// #include "barray-bones.hpp"

#ifndef BARRY_BARRAY_MEAT_OPERATORS_HPP
#define BARRY_BARRAY_MEAT_OPERATORS_HPP 1

#define BARRAY_TYPE() BArray<Cell_Type, Data_Type>

#define BARRAY_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BARRAY_TEMPLATE(a,b) \
    template BARRAY_TEMPLATE_ARGS() inline a BARRAY_TYPE()::b

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]

template BARRAY_TEMPLATE_ARGS()
inline void checkdim_(
    const BARRAY_TYPE()& lhs,
    const BARRAY_TYPE()& rhs
) {

    if (lhs.ncol() != rhs.ncol())
        throw std::length_error("Number of columns do not match.");

    if (lhs.nrow() != rhs.nrow())
        throw std::length_error("Number of rows do not match.");

    return;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator+=) (
    const BArray<Cell_Type, Data_Type>& rhs
) {

    // Must be compatible
    checkdim_(*this, rhs);
    
    for (size_t i = 0u; i < nrow(); ++i)
        for (size_t j = 0u; j < ncol(); ++j)
            this->operator()(i, j) += rhs.get_cell(i, j);

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator+=) (
    const Cell_Type& rhs
) {

    for (size_t i = 0u; i < nrow(); ++i) {
        for (size_t j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) += rhs;
        }
    }

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator-=) (
    const BArray<Cell_Type, Data_Type>& rhs
) {

    // Must be compatible
    checkdim_(*this, rhs);
    
    for (size_t i = 0u; i < nrow(); ++i) {
        for (size_t j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) -= rhs.get_cell(i, j);
        }
    }

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator-=) (
    const Cell_Type& rhs
) {

    for (size_t i = 0u; i < nrow(); ++i) {
        for (size_t j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) -= rhs;
        }
    }

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator*=) (
    const Cell_Type& rhs
) {

    for (size_t i = 0u; i < nrow(); ++i) {

        if (ROW(i).size() == 0u)
            continue;

        for (auto col = ROW(i).begin(); col != ROW(i).end(); ++col) {
            this->operator()(i, col->first) *= rhs;
        }
    }

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator/=) (
    const Cell_Type& rhs
) {

    for (size_t i = 0u; i < nrow(); ++i) {

        if (ROW(i).size() == 0u)
            continue;

        for (auto col = ROW(i).begin(); col != ROW(i).end(); ++col) {
            this->operator()(i, col->first) /= rhs;
        }
    }

    return *this;
}

#undef BARRAY_TYPE
#undef BARRAY_TEMPLATE_ARGS
#undef BARRAY_TEMPLATE

#undef ROW
#undef COL

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barray-meat-operators.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barraydense-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_BARRAYDENSE_BONES_HPP 
#define BARRY_BARRAYDENSE_BONES_HPP 1

template<typename Cell_Type, typename Data_Type>
class BArrayDenseRow;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseRow_const;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol_const;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCell;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCell_const;

/**
 * @brief Baseline class for binary arrays.
 * 
 * `BArrayDense` class objects are arbitrary dense-arrays. The data
 * is stored internally in the `el` member, which can be accessed
 * using the member function `get_data()`, by column.
 *
 * @tparam Cell_Type Type of cell (any type).
 * @tparam Data_Type Data type of the array (bool default).
 */
template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDense {
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCol<Cell_Type,Data_Type>;
    friend class BArrayDenseCol_const<Cell_Type,Data_Type>;
    friend class BArrayDenseRow<Cell_Type,Data_Type>;
    friend class BArrayDenseRow_const<Cell_Type,Data_Type>;
    // friend class Support<Cell_Type,Data_Type>;
    // friend class StatsCounter<Cell_Type,Data_Type>;
private:
    size_t N;
    size_t M;
    // size_t NCells = 0u;
    std::vector< Cell_Type > el;
    std::vector< Cell_Type > el_rowsums;
    std::vector< Cell_Type > el_colsums;
    Data_Type * data = nullptr;
    bool delete_data = false;

    static Cell_Type Cell_default;
    static const bool dense = true;

public:
    
    /** 
     * This is as a reference, if we need to iterate through the cells and we need
     * to keep track which were visited, we use this as a reference. So that if
     * cell.visited = true and visited = true, it means that we haven't been here
     * yet. Ideally, any routine using this->visited should switch it at the
     * beginning of the routine.
     */
    bool visited = false;
    

    /**
     * @name Constructors
     * 
     * @param N_ Number of rows
     * @param M_ Number of columns
     * @param source An unsigned vector ranging from 0 to N_
     * @param target An size_t vector ranging from 0 to M_
     * @param target When `true` tries to add repeated observations.
     * @param value Cell_Type defaul fill-in value (zero, by default.)
     */
    ///@{
    
    /** @brief Zero-size array */
    BArrayDense() : N(0u), M(0u), el(0u), el_rowsums(0u), el_colsums(0u) {};
    
    /** @brief Empty array */
    BArrayDense (size_t N_, size_t M_, Cell_Type value = static_cast<Cell_Type>(0)) :
        N(N_), M(M_), el(N_ * M_, value),
        el_rowsums(N_, static_cast<Cell_Type>(value * M_)), el_colsums(M_, static_cast<Cell_Type>(value * N_)) {};
    
    /** @brief Edgelist with data */
    BArrayDense (
        size_t N_,
        size_t M_,
        const std::vector< size_t > & source,
        const std::vector< size_t > & target,
        const std::vector< Cell_Type > & value,
        bool add = true
    );
    
    /** @brief Edgelist with no data (simpler) */
    BArrayDense (
        size_t N_, size_t M_,
        const std::vector< size_t > & source,
        const std::vector< size_t > & target,
        bool add = true
    );
    
    /** @brief Copy constructor */
    BArrayDense(const BArrayDense<Cell_Type,Data_Type> & Array_, bool copy_data = false);
    
    /** @brief Assignment constructor */
    BArrayDense<Cell_Type,Data_Type> & operator=(const BArrayDense<Cell_Type,Data_Type> & Array_);

    /** @brief Move operator */
    BArrayDense(BArrayDense<Cell_Type,Data_Type> && x) noexcept;

    /** @brief Move assignment */
    BArrayDense<Cell_Type,Data_Type> & operator=(BArrayDense<Cell_Type,Data_Type> && x) noexcept;
    ///@}
    
    bool operator==(const BArrayDense<Cell_Type,Data_Type> & Array_);

    ~BArrayDense();
    
    // In principle, copy can be faster by using openmp on the rows
    // since those are independent.
    // BArrayDense(BArrayDense & A);
    
    /**
     * @brief Set the data object
     * 
     * @param data_ 
     * @param delete_data_ 
     */
    ///@{
    void set_data(Data_Type * data_, bool delete_data_ = false);
    Data_Type * D_ptr();
    const Data_Type * D_ptr() const;
    Data_Type & D();
    const Data_Type & D() const;
    ///@}
    
    // Function to access the elements
    // bool check_cell
    void out_of_range(size_t i, size_t j) const;
    Cell_Type get_cell(size_t i, size_t j, bool check_bounds = true) const; 
    std::vector< Cell_Type >      get_col_vec(size_t i, bool check_bounds = true) const;
    std::vector< Cell_Type >      get_row_vec(size_t i, bool check_bounds = true) const;
    void                          get_col_vec(std::vector< Cell_Type > * x, size_t i, bool check_bounds = true) const;
    void                          get_row_vec(std::vector< Cell_Type > * x, size_t i, bool check_bounds = true) const;
    
    BArrayDenseRow<Cell_Type,Data_Type> & row(size_t i, bool check_bounds = true);
    const BArrayDenseRow_const<Cell_Type,Data_Type> row(size_t i, bool check_bounds = true) const;

    BArrayDenseCol<Cell_Type,Data_Type> & col(size_t j, bool check_bounds = true);
    const BArrayDenseCol_const<Cell_Type,Data_Type> col(size_t j, bool check_bounds = true) const;

    /**
     * @brief Get the edgelist
     * 
     * `Entries` is a class with three objects: Two `std::vector` with the row and
     * column coordinates respectively, and one `std::vector` with the corresponding
     * value of the cell.
     * 
     * @return Entries<Cell_Type> 
     */
    Entries<Cell_Type> get_entries() const;

    /**
     * @name Queries
     * @details `is_empty` queries a single cell. `nrow`, `ncol`, and `nnozero`
     * return the number of rows, columns, and non-zero cells respectively.
     * @param i,j Coordinates
     * @param check_bounds If `false` avoids checking bounds.
     */
    ///@{
    bool is_empty(size_t i, size_t j, bool check_bounds = true) const;
    size_t nrow() const noexcept;
    size_t ncol() const noexcept;
    size_t nnozero() const noexcept;
    Cell<Cell_Type> default_val() const;
    ///@}

    /**
     * @name Cell-wise insertion/deletion
     * @param i,j Row,column
     * @param check_bounds When `true` and out of range, the function throws an
     * error.
     * @param check_exists Wither check if the cell exists (before trying to
     * delete/add), or, in the case of `swap_cells`, check if either of both
     * cells exists/don't exist.
     */
    ///@{  
    BArrayDense<Cell_Type,Data_Type> & operator+=(const std::pair<size_t, size_t> & coords);
    BArrayDense<Cell_Type,Data_Type> & operator-=(const std::pair<size_t, size_t> & coords);
    BArrayDenseCell<Cell_Type,Data_Type> operator()(size_t i, size_t j, bool check_bounds = true);
    const Cell_Type operator()(size_t i, size_t j, bool check_bounds = true) const;
    
    void rm_cell(size_t i, size_t j, bool check_bounds = true, bool check_exists = true);
    
    void insert_cell(size_t i, size_t j, const Cell< Cell_Type > & v, bool check_bounds, bool);
    // void insert_cell(size_t i, size_t j, Cell< Cell_Type > && v, bool check_bounds, bool check_exists);
    void insert_cell(size_t i, size_t j, Cell_Type v, bool check_bounds, bool);
    
    void swap_cells(
        size_t i0, size_t j0, size_t i1, size_t j1, bool check_bounds = true,
        int check_exists = CHECK::BOTH,
        int * report     = nullptr
        );
    
    void toggle_cell(size_t i, size_t j, bool check_bounds = true, int check_exists = EXISTS::UKNOWN);
    void toggle_lock(size_t i, size_t j, bool check_bounds = true);
    ///@}
    
    /**@name Column/row wise interchange*/
    ///@{
    void swap_rows(size_t i0, size_t i1, bool check_bounds = true);
    void swap_cols(size_t j0, size_t j1, bool check_bounds = true);
    
    void zero_row(size_t i, bool check_bounds = true);
    void zero_col(size_t j, bool check_bounds = true);
    ///@}
    
    void transpose();
    void clear(bool hard = true);
    void resize(size_t N_, size_t M_);
    void reserve();

    // Advances operators
    // void toggle_iterator
    
    // Misc
    void print(const char * fmt = nullptr, ...) const;

    /**
     * @name Arithmetic operators
     * 
     */
    ///@{
    BArrayDense<Cell_Type,Data_Type>& operator+=(const BArrayDense<Cell_Type,Data_Type>& rhs);
    BArrayDense<Cell_Type,Data_Type>& operator+=(const Cell_Type & rhs);

    BArrayDense<Cell_Type,Data_Type>& operator-=(const BArrayDense<Cell_Type,Data_Type>& rhs);
    BArrayDense<Cell_Type,Data_Type>& operator-=(const Cell_Type & rhs);
    
    BArrayDense<Cell_Type,Data_Type>& operator/=(const Cell_Type & rhs);
    BArrayDense<Cell_Type,Data_Type>& operator*=(const Cell_Type & rhs);
    ///@}
    
    // /**
    //  * @name Casting between types
    //  */
    // ///@{
    // operator BArrayDense<double,bool>() const;
    // operator BArrayDense<int,bool>() const;
    // operator BArrayDense<size_t,bool>() const;
    // operator BArrayDense<bool,bool>() const;
    // ///@}
    
    bool is_dense() const noexcept {return dense;};

    const std::vector< Cell_Type > & get_data() const;
    const Cell_Type rowsum(size_t i) const;
    const Cell_Type colsum(size_t i) const;
};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barraydense-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barraydensecell-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include "typedefs.hpp"

#ifndef BARRY_BARRAYDENSECELL_BONES_HPP
#define BARRY_BARRAYDENSECELL_BONES_HPP 1

#define POS(a, b) (a) + (b) * N

template<typename Cell_Type, typename Data_Type>
class BArrayDense;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol_const;

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCell {
    friend class BArrayDense<Cell_Type,Data_Type>;
    friend class BArrayDenseCol<Cell_Type,Data_Type>;
    friend class BArrayDenseCol_const<Cell_Type,Data_Type>;
private:
  
    BArrayDense<Cell_Type,Data_Type> * dat;
    size_t i;
    size_t j;
  
public:
  
    BArrayDenseCell(
        BArrayDense<Cell_Type,Data_Type> * Array_,
        size_t i_,
        size_t j_,
        bool check_bounds = true
        ) : 
    i(i_), j(j_)
    {

        if (check_bounds)
        {

            if (i >= Array_->nrow())
                throw std::length_error("Row out of range.");
            if (j >= Array_->ncol())
                throw std::length_error("Col out of range.");

        }
        dat = Array_;

    };

    BArrayDenseCell<Cell_Type,Data_Type>& operator=(
        const BArrayDenseCell<Cell_Type,Data_Type> & other
        );

    ~BArrayDenseCell(){};
    void operator=(const Cell_Type & val);
    void operator+=(const Cell_Type & val);
    void operator-=(const Cell_Type & val);
    void operator*=(const Cell_Type & val);
    void operator/=(const Cell_Type & val);

    operator Cell_Type() const;
    bool operator==(const Cell_Type & val) const;
  
};


#undef POS

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barraydensecell-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barraydenserow-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_BARRAYDENSEROW_BONES_HPP
#define BARRY_BARRAYDENSEROW_BONES_HPP

#define POS(a,b) (b) * N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)
#define ZERO_CELL static_cast< Cell_Type >(0.0)

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseRow {
    friend class BArrayDense<Cell_Type,Data_Type>;
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    BArrayDense< Cell_Type,Data_Type > * array;
    Row_type< Cell_Type > row;
    size_t index;
    bool row_filled = false; // after row is filled

    void fill_if_needed()
    {
        if (!row_filled)
        {

            for (size_t j = 0u; j < array->M; ++j)
            {
                
                if (array->el[POS_N(index, j, array->N)] != ZERO_CELL)
                    row[j] = row[POS_N(index, j, array->N)];
                    
            }

            row_filled = true;
            
        }
    }


public:

    BArrayDenseRow(
        BArrayDense< Cell_Type,Data_Type > & array_,
        size_t i
    ) : array(&array_), index(i) {};

    typename Row_type<Cell_Type>::iterator & begin()
    {

        fill_if_needed();
        return row.begin();

    };

    typename Row_type<Cell_Type>::iterator & end()
    {

        fill_if_needed();
        return row.end();

    };

    size_t size() const noexcept
    {

        fill_if_needed();
        return row.size();

    };

    std::pair<size_t,Cell<Cell_Type>> & operator()(size_t i)
    {

        fill_if_needed();
        return row[i];

    }

};

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseRow_const {
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    const BArrayDense< Cell_Type,Data_Type > * array;
    Row_type< Cell_Type > row;
    size_t index;

public:
    BArrayDenseRow_const(
        const BArrayDense< Cell_Type,Data_Type > & array_,
        size_t i
    ) : array(&array_), index(i)
    {

        for (size_t j = 0u; j < array->M; ++j)
        {
            
            if (array->el[POS_N(index, j, array->M)] != ZERO_CELL)
                row[j] = row[POS_N(index, j, array->M)];
                
        }

        return;


    };

    typename Row_type< Cell_Type >::const_iterator begin() const
    {
        return row.begin();
    };

    typename Row_type< Cell_Type >::const_iterator end() const
    {
        return row.end();
    };

    size_t size() const noexcept
    {
        return row.size();
    };

    const std::pair<size_t,Cell<Cell_Type>> operator()(size_t i) const
    {
        return row[i];
    }

};

#undef POS
#undef POS_N
#undef ZERO_CELL

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barraydenserow-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barraydensecol-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_BARRAYDENSECOL_BONES 
#define BARRY_BARRAYDENSECOL_BONES

#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)
#define ZERO_CELL static_cast<Cell_Type>(0.0)

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCol {
    friend class BArrayDense<Cell_Type,Data_Type>;
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    BArrayDense< Cell_Type,Data_Type > * array;
    Col_type<Cell_Type> col;
    size_t index;
    bool col_filled = false;

    void fill_if_needed()
    {
        if (!col_filled)
        {

            for (size_t i = 0u; i < array->N; ++i)
            {
                
                if (array->el[POS_N(i, index, array->N)] != ZERO_CELL)
                    col[i] = col[POS_N(i, index, array->N)];
                    
            }

            col_filled = true;
            
        }
    }

public:
    BArrayDenseCol(
        BArrayDense< Cell_Type,Data_Type > & array_,
        size_t j
    ) : array(&array_), index(j) {};


    typename Col_type<Cell_Type>::iterator & begin()
    {
        fill_if_needed();
        return col.begin();
    };

    typename Col_type<Cell_Type>::iterator & end()
    {
        fill_if_needed();
        return col.end();
    };

    size_t size() const noexcept
    {
        fill_if_needed();
        return col.size();
    };

    std::pair<size_t,Cell_Type*> & operator()(size_t i)
    {
        fill_if_needed();
        return col[i];
    }

};

template <typename Cell_Type = bool, typename Data_Type = bool>
class BArrayDenseCol_const {
    friend class BArrayDenseCell<Cell_Type,Data_Type>;
    friend class BArrayDenseCell_const<Cell_Type,Data_Type>;
private:
    const BArrayDense< Cell_Type,Data_Type > * array;
    size_t index;
    Col_type<Cell_Type> col;

public:
    BArrayDenseCol_const(
        const BArrayDense< Cell_Type,Data_Type > & array_,
        size_t j
    ) : array(&array_), index(j)
    {

        for (size_t i = 0u; i < array->N; ++i)
        {
            
            if (array->el[POS_N(i, index, array->N)] != ZERO_CELL)
                col[i] = col[POS_N(i, index, array->N)];
                
        }

    };

    typename Col_type<Cell_Type>::iterator begin()
    {
        return col.begin();
    };

    typename Col_type<Cell_Type>::iterator end()
    {
        return col.end();
    };


    size_t size() const noexcept
    {
        return col.size();
    };

    const std::pair<size_t,Cell_Type*> operator()(size_t i) const
    {
        return col[i];
    }

};

#undef POS
#undef POS_N
#undef ZERO_CELL

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barraydensecol-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barraydense-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include <stdexcept>
// #include "barraydense-bones.hpp"

#ifndef BARRY_BARRAYDENSE_MEAT_HPP
#define BARRY_BARRAYDENSE_MEAT_HPP 

template<typename Cell_Type, typename Data_Type>
class BArrayDenseRow;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseRow_const;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCol_const;

template<typename Cell_Type, typename Data_Type>
class BArrayDenseCell;


#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]
#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)

template<typename Cell_Type, typename Data_Type>
Cell_Type BArrayDense<Cell_Type,Data_Type>::Cell_default = static_cast< Cell_Type >(1.0); 

#define ZERO_CELL static_cast<Cell_Type>(0.0)

// Edgelist with data
template<typename Cell_Type, typename Data_Type>
inline BArrayDense<Cell_Type, Data_Type>::BArrayDense(
    size_t N_,
    size_t M_,
    const std::vector< size_t > & source,
    const std::vector< size_t > & target,
    const std::vector< Cell_Type > & value,
    bool add
) : N(N_), M(M_), el(N_ * M_, ZERO_CELL), el_rowsums(N_, ZERO_CELL), el_colsums(M_, ZERO_CELL) {
  
    if (source.size() != target.size())
        throw std::length_error("-source- and -target- don't match on length.");
    if (source.size() != value.size())
        throw std::length_error("-sorce- and -value- don't match on length.");
    
    // Writing the data
    for (size_t i = 0u; i < source.size(); ++i)
    {
      
        // Checking range
        bool empty = is_empty(source[i], target[i], true);
        if (add && !empty)
        {

            Cell_Type tmp = el[POS(source[i], target[i])];
            
            el_rowsums[source[i]] += (value[i] - tmp);
            el_colsums[target[i]] += (value[i] - tmp);

            el[POS(source[i], target[i])] += value[i];
            
            continue;

        } 
        
        if (!empty)
            throw std::logic_error("The value already exists. Use 'add = true'.");
          
        el[POS(source[i], target[i])] = value[i];

        el_rowsums[source[i]] += value[i];
        el_colsums[target[i]] += value[i];
        

    }
    
    return;
  
}

// Edgelist without data
template<typename Cell_Type, typename Data_Type>
inline BArrayDense<Cell_Type, Data_Type>:: BArrayDense(
    size_t N_, size_t M_,
    const std::vector< size_t > & source,
    const std::vector< size_t > & target,
    bool add
) : N(N_), M(M_), el(N_ * M_, ZERO_CELL), el_rowsums(N_, ZERO_CELL), el_colsums(M_, ZERO_CELL) {
  
    std::vector< Cell_Type > value(source.size(), static_cast<Cell_Type>(1.0));

    if (source.size() != target.size())
        throw std::length_error("-source- and -target- don't match on length.");
    if (source.size() != value.size())
        throw std::length_error("-sorce- and -value- don't match on length.");

    
    // Writing the data
    for (size_t i = 0u; i < source.size(); ++i)
    {
      
        // Checking range
        bool empty = is_empty(source[i], target[i], true);
        if (add && !empty)
        {

            Cell_Type tmp = el[POS(source[i], target[i])];
            
            el_rowsums[source[i]] += (value[i] - tmp);
            el_colsums[target[i]] += (value[i] - tmp);

            el[POS(source[i], target[i])] += value[i];
            
            continue;

        } 
        
        if (!empty)
            throw std::logic_error("The value already exists. Use 'add = true'.");
          
        el[POS(source[i], target[i])] = value[i];

        el_rowsums[source[i]] += value[i];
        el_colsums[target[i]] += value[i];
        

    }
  
}

template<typename Cell_Type, typename Data_Type>
inline BArrayDense<Cell_Type, Data_Type>:: BArrayDense(
    const BArrayDense<Cell_Type, Data_Type> & Array_,
    bool copy_data
) : N(Array_.N), M(Array_.M){
  
    // Dimensions
    el = Array_.el;
    el_rowsums = Array_.el_rowsums;
    el_colsums = Array_.el_colsums;
    // el.resize(0u);
    // el_rowsums.resize(0u);
    // el_colsums.resize(0u);
    
    // std::copy(Array_.el.begin(), Array_.el.end(), std::back_inserter(el));
    // std::copy(Array_.el_rowsums.begin(), Array_.el_rowsums.end(), std::back_inserter(el_rowsums));
    // std::copy(Array_.el_colsums.begin(), Array_.el_colsums.end(), std::back_inserter(el_colsums));

    // this->NCells  = Array_.NCells;
    this->visited = Array_.visited;
    
    // Data
    if (Array_.data != nullptr)
    {

        if (copy_data)
        {

            data = new Data_Type(*Array_.data);
            delete_data = true;

        } else {

            data = Array_.data;
            delete_data = false;

        }

    }
    
    return;
  
}

template<typename Cell_Type, typename Data_Type>
inline BArrayDense<Cell_Type,Data_Type> & BArrayDense<Cell_Type, Data_Type>::operator=(
    const BArrayDense<Cell_Type, Data_Type> & Array_
) {
  
    // Clearing
    if (this != &Array_)
    {
      
        el = Array_.el;
        el_rowsums = Array_.el_rowsums;
        el_colsums = Array_.el_colsums;
        // el.resize(0u);
        // el_rowsums.resize(0u);
        // el_colsums.resize(0u);
        
        // // Entries
        // std::copy(Array_.el.begin(), Array_.el.end(), std::back_inserter(el));
        // std::copy(Array_.el_rowsums.begin(), Array_.el_rowsums.end(), std::back_inserter(el_rowsums));
        // std::copy(Array_.el_colsums.begin(), Array_.el_colsums.end(), std::back_inserter(el_colsums));


        // this->NCells = Array_.NCells;
        this->N      = Array_.N;
        this->M      = Array_.M;
      
        // Data
        if (data != nullptr)
        {

            if (delete_data)
                delete data;
            data = nullptr;

        }

        if (Array_.data != nullptr)
        {

            data = new Data_Type(*Array_.data);
            delete_data = true;

        }
      
    }
      
    return *this;
  
}

template<typename Cell_Type, typename Data_Type>
inline BArrayDense<Cell_Type, Data_Type>:: BArrayDense(
    BArrayDense<Cell_Type, Data_Type> && x
    ) noexcept :
    N(std::move(x.N)), M(std::move(x.M)),
    // NCells(std::move(x.NCells)),
    el(std::move(x.el)),
    el_rowsums(std::move(x.el_rowsums)),
    el_colsums(std::move(x.el_colsums)),
    data(std::move(x.data)),
    delete_data(std::move(x.delete_data))
{

      x.data        = nullptr;
      x.delete_data = false;

}

template<typename Cell_Type, typename Data_Type>
inline BArrayDense<Cell_Type,Data_Type> & BArrayDense<Cell_Type, Data_Type>::operator=(
    BArrayDense<Cell_Type, Data_Type> && x
) noexcept {
  
    // Clearing
    if (this != &x)
    {
      
        N      = x.N;
        M      = x.M;
        // NCells = x.NCells;
        
        std::swap(el, x.el);
        std::swap(el_rowsums, x.el_rowsums);
        std::swap(el_colsums, x.el_colsums);
              
        // Data
        if (data != nullptr)
        {

            if (delete_data)
                delete data;
            data = nullptr;

        }

        if (x.data != nullptr)
        {

            data        = std::move(x.data);
            delete_data = x.delete_data;

            x.delete_data = false;
            x.data = nullptr;

        }
      
    }
      
    return *this;
  
}

template<typename Cell_Type, typename Data_Type>
inline bool BArrayDense<Cell_Type, Data_Type>::operator== (
    const BArrayDense<Cell_Type, Data_Type> & Array_
) {
    
    // Dimension and number of cells used
    if ( (N != Array_.nrow()) | (M != Array_.ncol()) )
        return false;
    
    // One holds, and the other doesn't.
    if ((!data & Array_.data) | (data & !Array_.data))
        return false;
    
    if (this->el != Array_.el)
        return false;
    
    return true;
}

template<typename Cell_Type, typename Data_Type>
inline BArrayDense<Cell_Type, Data_Type>::~BArrayDense () {
    
    if (delete_data && (data != nullptr))
        delete data;
    
    return;
}

template<typename Cell_Type, typename Data_Type>
inline void BArrayDense<Cell_Type, Data_Type>::set_data (
    Data_Type * data_,
    bool delete_data_
) {  

    if ((data != nullptr) && delete_data)
        delete data;
    
    data        = data_;
    delete_data = delete_data_;
    
    return;
    
}

template<typename Cell_Type, typename Data_Type>
inline Data_Type * BArrayDense<Cell_Type, Data_Type>::D_ptr () {
    return this->data;
}

template<typename Cell_Type, typename Data_Type>
inline const Data_Type * BArrayDense<Cell_Type, Data_Type>::D_ptr () const {
    return this->data;
}

template<typename Cell_Type, typename Data_Type>
 inline Data_Type & BArrayDense<Cell_Type, Data_Type>::D () {
    return *this->data;
}

template<typename Cell_Type, typename Data_Type>
inline const Data_Type & BArrayDense<Cell_Type, Data_Type>::D () const {
    return *this->data;
}

template<typename Cell_Type, typename Data_Type>
inline void BArrayDense<Cell_Type, Data_Type>::out_of_range (
    size_t i,
    size_t j
) const {

    if (i >= N)
    {
        std::string err_msg = "The row is out of range: " + std::to_string(i) + " >= " + std::to_string(N);
        throw std::range_error(err_msg);

    } else if (j >= M)
    {
        std::string err_msg = "The column is out of range: " + std::to_string(j) + " >= " + std::to_string(M);
        throw std::range_error(err_msg);
    }

    return;

}
    
template<typename Cell_Type, typename Data_Type>
inline Cell_Type BArrayDense<Cell_Type, Data_Type>::get_cell (
    size_t i,
    size_t j,
    bool check_bounds
) const {
    
    // Checking boundaries  
    if (check_bounds)
        out_of_range(i,j);
    
    return el[POS(i, j)];
    
}

template<typename Cell_Type, typename Data_Type>
inline std::vector< Cell_Type > BArrayDense<Cell_Type, Data_Type>::get_row_vec (
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    std::vector< Cell_Type > ans;
    ans.reserve(ncol());
    for (size_t j = 0u; j < M; ++j) 
        ans.push_back(el[POS(i, j)]);
    
    return ans;

}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: get_row_vec (
    std::vector<Cell_Type> * x,
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    for (size_t j = 0u; j < M; ++j) 
        x->operator[](j) = el[POS(i, j)];
    
}

template<typename Cell_Type, typename Data_Type> inline std::vector< Cell_Type > BArrayDense<Cell_Type, Data_Type>:: get_col_vec(
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    std::vector< Cell_Type > ans;
    ans.reserve(nrow()); 
    for (size_t j = 0u; j < N; ++j) 
        ans.push_back(el[POS(j, i)]);
    
    return ans;

}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: get_col_vec (
    std::vector<Cell_Type> * x,
    size_t i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    #ifdef __INTEL_LLVM_COMPILER
    #pragma code_align 32
    #endif
    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp simd
    #endif
    for (size_t j = 0u; j < N; ++j) 
        x->operator[](j) = el[POS(j, i)];//this->get_cell(iter->first, i, false);
    
}
template<typename Cell_Type, typename Data_Type>
inline const BArrayDenseRow_const<Cell_Type,Data_Type> BArrayDense<Cell_Type, Data_Type>::row(
    size_t i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, 0u);

    return BArrayDenseRow_const<Cell_Type,Data_Type>(*this, i);

}

template<typename Cell_Type, typename Data_Type>
inline BArrayDenseRow<Cell_Type,Data_Type> & BArrayDense<Cell_Type, Data_Type>::row(
    size_t i,
    bool check_bounds
) {

    if (check_bounds)
        out_of_range(i, 0u);

    return BArrayDenseRow<Cell_Type,Data_Type>(*this, i);

}

template<typename Cell_Type, typename Data_Type>
inline const BArrayDenseCol_const<Cell_Type,Data_Type> 
BArrayDense<Cell_Type,Data_Type>::col(
    size_t j,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(0u, j);

    return BArrayDenseCol_const<Cell_Type,Data_Type>(*this, j);

}

template<typename Cell_Type, typename Data_Type>
inline BArrayDenseCol<Cell_Type,Data_Type> & 
BArrayDense<Cell_Type,Data_Type>::col(
    size_t j,
    bool check_bounds
) {

    if (check_bounds)
        out_of_range(0u, j);

    return BArrayDenseCol<Cell_Type,Data_Type>(*this, j);

}

template<typename Cell_Type, typename Data_Type> inline Entries< Cell_Type > BArrayDense<Cell_Type, Data_Type>:: get_entries() const {
    
    size_t nzero = this->nnozero();

    Entries<Cell_Type> res(nzero);
    
    for (size_t i = 0u; i < N; ++i)
    {
        for (size_t j = 0u; j < M; ++j)
        {

            if (el[POS(i, j)] != BARRY_ZERO_DENSE)
            {

                res.source.push_back(i),
                res.target.push_back(j),
                res.val.push_back(el[POS(i, j)]);

            }
            

        }

    }
    
    return res;

}

template<typename Cell_Type, typename Data_Type> inline bool BArrayDense<Cell_Type, Data_Type>:: is_empty(
    size_t i,
    size_t j,
    bool check_bounds
) const {
    
    if (check_bounds)
        out_of_range(i, j);
    
    return el[POS(i, j)] == ZERO_CELL;
    
}

template<typename Cell_Type, typename Data_Type> inline size_t BArrayDense<Cell_Type, Data_Type>:: nrow() const noexcept {
    return N;
}

template<typename Cell_Type, typename Data_Type> inline size_t BArrayDense<Cell_Type, Data_Type>:: ncol() const noexcept {
    return M;
}

template<typename Cell_Type, typename Data_Type> inline size_t BArrayDense<Cell_Type, Data_Type>:: nnozero() const noexcept {

    size_t nzero = 0u;
    for (auto & v : el)
        if (v != BARRY_ZERO_DENSE)
            nzero++;

    return nzero;
}

template<typename Cell_Type, typename Data_Type>
inline Cell< Cell_Type> BArrayDense<Cell_Type, Data_Type>::default_val() const {
    return this->Cell_default;
}

template<typename Cell_Type, typename Data_Type>
inline BArrayDense<Cell_Type,Data_Type> & BArrayDense<Cell_Type, Data_Type>::operator+=(
    const std::pair<size_t,size_t> & coords
) {
    

    size_t i = coords.first;
    size_t j = coords.second;

    out_of_range(i, j);

    el[POS(i,j)]  += 1;
    el_rowsums[i] += 1;
    el_colsums[j] += 1;
    
    return *this;
    
}

template<typename Cell_Type, typename Data_Type>
inline BArrayDense<Cell_Type,Data_Type> & BArrayDense<Cell_Type, Data_Type>::operator-=(
    const std::pair<size_t,size_t> & coords
) {
    
    size_t i = coords.first;
    size_t j = coords.second;

    out_of_range(i, j);

    Cell_Type old = el[POS(i,j)];

    el[POS(i,j)]   = ZERO_CELL;
    el_rowsums[i] -= old;
    el_colsums[j] -= old;
    
    return *this;
    
}

template<typename Cell_Type, typename Data_Type>
inline BArrayDenseCell<Cell_Type,Data_Type> BArrayDense<Cell_Type, Data_Type>::operator()(  
    size_t i,
    size_t j,
    bool check_bounds
) {
    
    return BArrayDenseCell<Cell_Type,Data_Type>(this, i, j, check_bounds);
    
}

template<typename Cell_Type, typename Data_Type>
inline const Cell_Type BArrayDense<Cell_Type, Data_Type>::operator()(  
    size_t i,
    size_t j,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, j);
    
    return el[POS(i,j)];
    
}

template<typename Cell_Type, typename Data_Type>
inline void BArrayDense<Cell_Type, Data_Type>::rm_cell (
    size_t i,
    size_t j,
    bool check_bounds,
    bool check_exists
) {
    
    // Checking the boundaries
    if (check_bounds)
        out_of_range(i,j);

    // BARRY_UNUSED(check_exists)
        
    // Remove the pointer first (so it wont point to empty)
    el_rowsums[i] -= el[POS(i, j)];
    el_colsums[j] -= el[POS(i, j)];    
    el[POS(i, j)] = BARRY_ZERO_DENSE;
    
    return;

}

template<typename Cell_Type, typename Data_Type>
inline void BArrayDense<Cell_Type, Data_Type>::insert_cell (
    size_t i,
    size_t j,
    const Cell< Cell_Type> & v,
    bool check_bounds,
    bool
) { 
    
    if (check_bounds)
        out_of_range(i,j); 

    if (el[POS(i,j)] == BARRY_ZERO_DENSE)
    {

        el_rowsums[i] += v.value;
        el_colsums[j] += v.value;
        
    } 
    else
    {

        Cell_Type old = el[POS(i,j)];
        el_rowsums[i] += (v.value - old);
        el_colsums[j] += (v.value - old);

    }

    el[POS(i, j)] = v.value;

    return;

    
}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: insert_cell(
    size_t i,
    size_t j,
    Cell_Type v,
    bool check_bounds,
    bool
) {
    
    if (check_bounds)
        out_of_range(i,j);
        
    if (el[POS(i,j)] == BARRY_ZERO_DENSE)
    {

        el_rowsums[i] += v;
        el_colsums[j] += v;
        
    } 
    else
    {

        Cell_Type old = el[POS(i,j)];
        el_rowsums[i] += (v - old);
        el_colsums[j] += (v - old);

    }

    el[POS(i, j)] = v;

}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: swap_cells (
        size_t i0, size_t j0,
        size_t i1, size_t j1,
        bool check_bounds,
        int check_exists,
        int * report
) {
    
    if (check_bounds) {
        out_of_range(i0,j0);
        out_of_range(i1,j1);
    }
    
        
    // Just in case, if this was passed
    if (report != nullptr)
        (*report) = EXISTS::BOTH;
    
    // If source and target coincide, we do nothing
    if ((i0 == i1) && (j0 == j1)) 
        return;

    // Updating rowand col sumns    
    Cell_Type val0 = el[POS(i0,j0)];
    Cell_Type val1 = el[POS(i1,j1)];

    rm_cell(i0, j0, false, false);
    rm_cell(i1, j1, false, false);
    
    // Inserting the cells by reference, these will be deleted afterwards
    insert_cell(i0, j0, val1, false, false);
    insert_cell(i1, j1, val0, false, false);
    
    return;

}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: toggle_cell (
    size_t i,
    size_t j,
    bool check_bounds,
    int check_exists
) {

    if (check_bounds)
        out_of_range(i, j);

    if (el[POS(i,j)] == ZERO_CELL)
        insert_cell(i,j,1,false,false);
    else
        rm_cell(i,j,false,false);
    
    return;
    
}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: swap_rows (
    size_t i0,
    size_t i1,
    bool check_bounds
) {
  
    if (check_bounds)
    {

        out_of_range(i0,0u);
        out_of_range(i1,0u);

    }
     
    // if (NCells == 0u)
    //     return;
    
    // Swapping happens naturally, need to take care of the pointers
    // though
    for (size_t j = 0u; j < M; ++j)
        std::swap(el[POS(i0, j)], el[POS(i1, j)]);

    std::swap(el_rowsums[i0], el_rowsums[i1]);
    
    return;
}

// This swapping is more expensive overall
template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: swap_cols (
    size_t j0,
    size_t j1,
    bool check_bounds
) {

    if (check_bounds)
    {

        out_of_range(0u, j0);
        out_of_range(0u, j1);

    }
    
    if ((el_colsums[j0] == ZERO_CELL) && el_colsums[j1] == ZERO_CELL)
        return;

    // Swapping happens naturally, need to take care of the pointers
    // though
    for (size_t i = 0u; i < N; ++i)
        std::swap(el[POS(i, j0)], el[POS(i, j1)]);

    std::swap(el_colsums[j0], el_colsums[j1]);
    
    return;
}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: zero_row (
    size_t i,
    bool check_bounds
    ) {
  
    if (check_bounds)
        out_of_range(i, 0u);

    if (el_rowsums[i] == ZERO_CELL)
        return;

    // Else, remove all elements
    for (size_t col = 0u; col < M; col++) 
        rm_cell(i, col, false, false);
    
    return;
  
}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: zero_col (
    size_t j,
    bool check_bounds
  ) {
  
    if (check_bounds)
        out_of_range(0u, j);
    
    if (el_colsums[j] == ZERO_CELL)
        return;
    
    // Else, remove all elements
    for (size_t row = 0u; row < N; row++) 
        rm_cell(row, j, false, false);
    
    return;
  
}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: transpose () {
  
    // if (NCells == 0u)
    // {

    //     std::swap(N, M);
    //     return;

    // }

    // Start by flipping the switch 
    visited = !visited;

    // size_t N0 = N, M0 = M;
    std::vector< Cell< Cell_Type > > tmp_el(std::move(el));
    el.resize(N * M, ZERO_CELL);
    for (size_t i = 0u; i < N; ++i) 
        for (size_t j = 0u; j < M; ++j)
            std::swap(tmp_el[POS(i, j)], el[POS_N(j, i, M)]);
    
    // Swapping the values
    std::swap(N, M);
    std::swap(el_rowsums, el_colsums);
    
    return;

}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: clear (
    bool hard
) {
    
    BARRY_UNUSED(hard)
    
    std::fill(el.begin(), el.end(), ZERO_CELL);
    std::fill(el_rowsums.begin(), el_rowsums.end(), ZERO_CELL);
    std::fill(el_colsums.begin(), el_colsums.end(), ZERO_CELL);
    
    return;
    
}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: resize (
    size_t N_,
    size_t M_
) {

    // Moving stuff around
    std::vector< Cell_Type > el_tmp(el);
    el.resize(N_ * M_, ZERO_CELL);
    el_rowsums.resize(N_, ZERO_CELL);
    el_colsums.resize(M_, ZERO_CELL);

    for (size_t i = 0u; i < N; ++i)
    {
        // If reached the end
        if (i >= N_)
            break;

        for (size_t j = 0u; j < M; ++j)
        {

            if (j >= M_)
                break;

            insert_cell(i, j, el_tmp[POS_N(i, j, N_)], false, false);

        }

    }

    N = N_;
    M = M_;
    
    return;

}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: reserve () {

    el.reserve(N * M);
    el_rowsums.reserve(N);
    el_colsums.reserve(M);
    return;
  
}

template<typename Cell_Type, typename Data_Type> inline void BArrayDense<Cell_Type, Data_Type>:: print (
    const char * fmt,
    ...
) const
{
  
    std::va_list args;
    va_start(args, fmt);
    if (fmt != nullptr)
        printf_barry(fmt, args);
    va_end(args);

    for (size_t i = 0u; i < N; ++i)
    {

        printf_barry("[%3i,] ", static_cast<int>(i));

        for (size_t j = 0u; j < M; ++j)
        {

            if (this->is_empty(i, j, false))
            {
                printf_barry("    . ");
            } else {
                printf_barry(" %.2f ", static_cast<double>(this->get_cell(i, j, false)));
            }
            
        }

        printf_barry("\n");

    }
    
    return;
    
}

template<typename Cell_Type, typename Data_Type> inline const std::vector< Cell_Type > & BArrayDense<Cell_Type, Data_Type>:: get_data() const
{
    return el;
}

template<typename Cell_Type, typename Data_Type> inline const Cell_Type BArrayDense<Cell_Type, Data_Type>:: rowsum(size_t i) const
{
    return el_rowsums[i];
}

template<typename Cell_Type, typename Data_Type> inline const Cell_Type BArrayDense<Cell_Type, Data_Type>:: colsum(size_t j) const
{
    return el_colsums[j];
}

#undef ROW
#undef COL
#undef POS
#undef POS_N

#undef ZERO_CELL

#endif

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barraydense-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barraydensecell-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include "barraydensecell-bones.hpp"

#ifndef BARRY_BARRAYDENSECELL_MEAT_HPP
#define BARRY_BARRAYDENSECELL_MEAT_HPP 1

#define POS(a, b) (a) + (b) * dat->N 

template<typename Cell_Type,typename Data_Type>
inline BArrayDenseCell<Cell_Type,Data_Type>& BArrayDenseCell<Cell_Type,Data_Type>::operator=(
    const BArrayDenseCell<Cell_Type,Data_Type> & other
    ) {
    
    Cell_Type val = static_cast<Cell_Type>(other);
    #ifdef BARRY_DEBUG
    Cell_Type old      =  dat->el.at(POS(i,j));
    dat->el.at(POS(i,j))  =  val;
    dat->el_rowsums.at(i) += (val - old);
    dat->el_colsums.at(j) += (val - old);
    #else
    Cell_Type old      =  dat->el[POS(i,j)];
    dat->el[POS(i,j)]  =  val;
    dat->el_rowsums[i] += (val - old);
    dat->el_colsums[j] += (val - old);
    #endif

    return *this;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {

    #ifdef BARRY_DEBUG
    Cell_Type old      =  dat->el.at(POS(i,j));
    dat->el.at(POS(i,j))  =  val;
    dat->el_rowsums.at(i) += (val - old);
    dat->el_colsums.at(j) += (val - old);
    #else
    Cell_Type old      =  dat->el[POS(i,j)];
    dat->el[POS(i,j)]  =  val;
    dat->el_rowsums[i] += (val - old);
    dat->el_colsums[j] += (val - old);
    #endif
    
}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
    
    #ifdef BARRY_DEBUG
    dat->el.at(POS(i,j))  += val;
    dat->el_rowsums.at(i) += val;
    dat->el_colsums.at(j) += val;
    #else
    dat->el[POS(i,j)]  += val;
    dat->el_rowsums[i] += val;
    dat->el_colsums[j] += val;
    #endif

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
    
    #ifdef BARRY_DEBUG
    dat->el.at(POS(i,j))  -= val;
    dat->el_rowsums.at(i) -= val;
    dat->el_colsums.at(j) -= val;
    #else
    dat->el[POS(i,j)]  -= val;
    dat->el_rowsums[i] -= val;
    dat->el_colsums[j] -= val;
    #endif

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
    
    #ifdef BARRY_DEBUG
    Cell_Type old = dat->el.at(POS(i,j));
    dat->el_colsums.at(j) += (old * val - old);
    dat->el_rowsums.at(i) += (old * val - old);
    dat->el.at(POS(i,j)) *= val;
    #else
    Cell_Type old = dat->el[POS(i,j)];
    dat->el_colsums[j] += (old * val - old);
    dat->el_rowsums[i] += (old * val - old);
    dat->el[POS(i,j)] *= val;
    #endif

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
    
    #ifdef BARRY_DEBUG
    Cell_Type old = dat->el.at(POS(i,j));
    dat->el_rowsums.at(i) += (old/val - old);
    dat->el_colsums.at(j) += (old/val - old);
    dat->el.at(POS(i,j))  /= val;
    #else
    Cell_Type old = dat->el[POS(i,j)];
    dat->el_rowsums[i] += (old/val - old);
    dat->el_colsums[j] += (old/val - old);
    dat->el[POS(i,j)]  /= val;
    #endif

}

template<typename Cell_Type,typename Data_Type>
inline BArrayDenseCell<Cell_Type,Data_Type>::operator Cell_Type() const {
        return dat->el[POS(i,j)];
}

template<typename Cell_Type,typename Data_Type>
inline bool BArrayDenseCell<Cell_Type,Data_Type>::operator==(const Cell_Type & val) const {
    return dat->el[POS(i,j)] == val;  
}

#undef POS

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barraydensecell-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/barraydense-meat-operators.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


// #include <stdexcept>
// #include "barraydense-bones.hpp"

#ifndef BARRY_BARRAYDENSE_MEAT_OPERATORS_HPP
#define BARRY_BARRAYDENSE_MEAT_OPERATORS_HPP 1

#define BDENSE_TYPE() BArrayDense<Cell_Type, Data_Type>

#define BDENSE_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BDENSE_TEMPLATE(a,b) \
    template BDENSE_TEMPLATE_ARGS() inline a BDENSE_TYPE()::b

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]
#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)

template BDENSE_TEMPLATE_ARGS()
inline void checkdim_(
    const BDENSE_TYPE()& lhs,
    const BDENSE_TYPE()& rhs
) {

    if (lhs.ncol() != rhs.ncol())
        throw std::length_error("Number of columns do not match.");

    if (lhs.nrow() != rhs.nrow())
        throw std::length_error("Number of rows do not match.");

    return;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator+=) (
    const BDENSE_TYPE()& rhs
) {

    // Must be compatible
    checkdim_(*this, rhs);
    
    for (size_t i = 0u; i < nrow(); ++i)
        for (size_t j = 0u; j < ncol(); ++j)
            this->operator()(i, j) += rhs.get_cell(i, j);

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator+=) (
    const Cell_Type& rhs
) {

    for (size_t i = 0u; i < nrow(); ++i) {
        for (size_t j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) += rhs;
        }
    }

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator-=) (
    const BDENSE_TYPE()& rhs
) {

    // Must be compatible
    checkdim_(*this, rhs);
    
    for (size_t i = 0u; i < nrow(); ++i) {
        for (size_t j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) -= rhs.get_cell(i, j);
        }
    }

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator-=) (
    const Cell_Type& rhs
) {

    for (size_t i = 0u; i < nrow(); ++i) 
        for (size_t j = 0u; j < ncol(); ++j) 
            this->operator()(i, j) -= rhs;
        
    

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator*=) (
    const Cell_Type& rhs
) {

    for (size_t i = 0u; i < nrow(); ++i) 
        for (size_t j = 0u; j < nrow(); ++j)
            el[POS(i, j)] *= rhs;

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator/=) (
    const Cell_Type& rhs
) {

    for (size_t i = 0u; i < nrow(); ++i) 
        for (size_t j = 0u; j < nrow(); ++j)
            el[POS(i, j)] /= rhs;

    return *this;
}

#undef BDENSE_TYPE
#undef BDENSE_TEMPLATE_ARGS
#undef BDENSE_TEMPLATE

#undef ROW
#undef COL
#undef POS
#undef POS_N

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/barraydense-meat-operators.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


    
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/counters-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_COUNTERS_BONES_HPP
#define BARRY_COUNTERS_BONES_HPP 1

/**
 * @defgroup counting
 * @details `barry` includes a flexible way to generate counters based on change
 * statistics. Since most of the time we are counting many motifs in a graph,
 * change statistics make a reasonable (and efficient) way to make such counts.
 * 
 * In particular, let the motif be defined as \f$s(y)\f$, with \f$y\f$ as the
 * binary array. The change statistic when adding cell \f$y_{ij}\f$, i.e. when
 * the cell moves from being emty to have a one, is defined as
 * 
 * \f[
 * \delta(y_{ij}) = s^+_{ij}(y) - s^-_{ij}(y),
 * \f]
 * 
 * where \f$s^+_{ij}(y)\f$ and \f$s^-_{ij}(y)\f$ represent the motif statistic
 * with and without the ij-cell. For example, in the case of networks, the change
 * statistic for the number of edges is always 1. 
 * 
 * To count statistics in an array, the [Counter] class will empty the array, 
 * initialize the counters, and then start counting while adding at each step
 * a single cell, until matching the original array. 
 */

/**
  * @ingroup counting Implementation of motif counting
  * @brief A counter function based on change statistics.
  * 
  * This class is used by `CountStats` and `StatsCounter` as a way to count
  * statistics using change statistics.
  */
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class Counter {
public:
    
    Counter_fun_type<Array_Type,Data_Type> count_fun;
    Counter_fun_type<Array_Type,Data_Type> init_fun;
    Hasher_fun_type<Array_Type,Data_Type> hasher_fun;

    Data_Type data;
    std::string  name = "";
    std::string  desc = "";

    /**
     * @name Creator passing a counter and an initializer
     * 
     * @param count_fun_ The main counter function.
     * @param init_fun_ The initializer function can also be used to check if the
     *  `BArray` as the needed variables (see BArray::data).
     * @param data_ Data to be used with the counter.
     * @param delete_data_ When `true`, the destructor will delete the pointer
     * in the main data.
     */
    ///@{
    Counter() : count_fun(nullptr), init_fun(nullptr), hasher_fun(nullptr) {};
    
    Counter(
        Counter_fun_type<Array_Type,Data_Type> count_fun_,
        Counter_fun_type<Array_Type,Data_Type> init_fun_,
        Hasher_fun_type<Array_Type,Data_Type>  hasher_fun_,
        Data_Type                              data_,
        std::string                            name_        = "",   
        std::string                            desc_        = ""
        ): count_fun(count_fun_), init_fun(init_fun_), hasher_fun(hasher_fun_), data(data_),
            name(name_), desc(desc_) {};
    
    Counter(const Counter<Array_Type,Data_Type> & counter_); ///< Copy constructor
    Counter(Counter<Array_Type,Data_Type> && counter_) noexcept; ///< Move constructor
    Counter<Array_Type,Data_Type> operator=(const Counter<Array_Type,Data_Type> & counter_); ///< Copy assignment
    Counter<Array_Type,Data_Type>& operator=(Counter<Array_Type,Data_Type> && counter_) noexcept; ///< Move assignment
    ///@}

    ~Counter() {};
    
    /***
      * ! Main functions.
      */
    double count(Array_Type & Array, size_t i, size_t j);
    double init(Array_Type & Array, size_t i, size_t j);
    std::string get_name() const;
    std::string get_description() const;
    void set_name(std::string new_name);
    void set_description(std::string new_desc);

    /**
     * @brief Get and set the hasher function
     * 
     * The hasher function is used to characterize the support of the array.
     * This way, if possible, the support enumeration is recycled.
     * 
     * @param fun 
     */
    ///@{
    void set_hasher(Hasher_fun_type<Array_Type,Data_Type> fun);
    Hasher_fun_type<Array_Type,Data_Type> get_hasher();
    ///@}

    /**
     * @brief Print a summary of the counter.
     */
    void print() const;
    
};

/**
  * @brief Vector of counters.
  * 
  * Various functions hold more than one counter, so this class is a helper class
  * that allows managing multiple counters efficiently. The main data is a vector
  * to pointers of counters.
  */
template <typename Array_Type = BArray<>, typename Data_Type = bool>
class Counters {
    
private:
    std::vector< Counter<Array_Type,Data_Type > > data;
    Hasher_fun_type<Array_Type,Data_Type> hasher;
    
public: 
    
    // Constructors
    Counters();
    
    // Destructor needs to deal with the pointers
    ~Counters() {};

    /**
     * @brief Copy constructor
     * @param counter_ 
     */
    Counters(const Counters<Array_Type,Data_Type> & counter_);
    
    /**
     * @brief Move constructor
     * 
     * @param counters_ 
     */
    Counters(Counters<Array_Type,Data_Type> && counters_) noexcept;

    /**
     * @brief Copy assignment constructor
     * 
     * @param counter_ 
     * @return Counters<Array_Type,Data_Type> 
     */
    Counters<Array_Type,Data_Type> operator=(const Counters<Array_Type,Data_Type> & counter_);

    /**
     * @brief Move assignment constructor
     * 
     * @param counter_ 
     * @return Counters<Array_Type,Data_Type>& 
     */
    Counters<Array_Type,Data_Type> & operator=(Counters<Array_Type,Data_Type> && counter_) noexcept;
    
    /**
     * @brief Returns a pointer to a particular counter.
     * 
     * @param idx Id of the counter
     * @return Counter<Array_Type,Data_Type>* 
     */
    Counter<Array_Type,Data_Type> & operator[](size_t idx);

    /**
     * @brief Number of counters in the set.
     * 
     * @return size_t 
     */
    std::size_t size() const noexcept {
        return data.size();
        };
    
    // Functions to add counters
    void add_counter(Counter<Array_Type, Data_Type> counter);
    void add_counter(
        Counter_fun_type<Array_Type,Data_Type> count_fun_,
        Counter_fun_type<Array_Type,Data_Type> init_fun_,
        Hasher_fun_type<Array_Type,Data_Type>  hasher_fun_,
        Data_Type                              data_,
        std::string                            name_        = "",   
        std::string                            desc_        = ""
    );
    
    std::vector< std::string > get_names() const;
    std::vector< std::string > get_descriptions() const;

    /**
     * @brief Generates a hash for the given array according to the counters.
     * 
     * @param array 
     * @param add_dims When `true` (default) the dimmension of the array will
     * be added to the hash.
     * @return std::vector< double > That can be hashed later.
     */
    std::vector< double > gen_hash(
      const Array_Type & array,
      bool add_dims = true
      );

    /**
     * @brief Set the hasher function in addition to
     * the individual hasher functions of each counter.
     * @param fun_ A hasher function that will be appended to the
     * hash generated by the individual counters.
     */
    void add_hash(
      Hasher_fun_type<Array_Type,Data_Type> fun_
    );

    /**
     * @brief Print a summary of the counters in the set.
     * @param max_length_name Maximum length of the name to be printed.
     * @param max_length_desc Maximum length of the description to be printed.
     */
    void print(
      size_t max_length_name = 40,
      size_t max_length_desc = 40
    ) const;
    
};

#endif

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/counters-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/counters-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_COUNTERS_MEAT_HPP
#define BARRY_COUNTERS_MEAT_HPP 1

#define COUNTER_TYPE() Counter<Array_Type,Data_Type>

#define COUNTER_TEMPLATE_ARGS() <typename Array_Type, typename Data_Type>

#define COUNTER_TEMPLATE(a,b) \
    template COUNTER_TEMPLATE_ARGS() inline a COUNTER_TYPE()::b

COUNTER_TEMPLATE(,Counter)(
    const Counter<Array_Type,Data_Type> & counter_
) : count_fun(counter_.count_fun), init_fun(counter_.init_fun), hasher_fun(counter_.hasher_fun) {

    this->data = counter_.data;
    this->name = counter_.name;
    this->desc = counter_.desc;

    return;

}


COUNTER_TEMPLATE(,Counter)(
    Counter<Array_Type,Data_Type> && counter_
    ) noexcept :
    count_fun(std::move(counter_.count_fun)),
    init_fun(std::move(counter_.init_fun)),
    hasher_fun(std::move(counter_.hasher_fun)),
    data(std::move(counter_.data)),
    name(std::move(counter_.name)),
    desc(std::move(counter_.desc))
{

} ///< Move constructor

COUNTER_TEMPLATE(COUNTER_TYPE(),operator=)(
    const Counter<Array_Type,Data_Type> & counter_
)
{

    if (this != &counter_) {

        this->count_fun = counter_.count_fun;
        this->init_fun = counter_.init_fun;
        this->hasher_fun = counter_.hasher_fun;

        
        this->data = counter_.data;
        this->name = counter_.name;
        this->desc = counter_.desc;

    }

    return *this;

}

COUNTER_TEMPLATE(COUNTER_TYPE() &,operator=)(
    Counter<Array_Type,Data_Type> && counter_
) noexcept {

    if (this != &counter_)
    {

        this->data = std::move(counter_.data);

        // Functions
        this->count_fun = std::move(counter_.count_fun);
        this->init_fun = std::move(counter_.init_fun);
        this->hasher_fun = std::move(counter_.hasher_fun);

        // Descriptions
        this->name = std::move(counter_.name);
        this->desc = std::move(counter_.desc);

    }

    return *this;

} ///< Move assignment

COUNTER_TEMPLATE(double, count)(Array_Type & Array, size_t i, size_t j)
{

    if (count_fun == nullptr)
        return 0.0;

    return count_fun(Array, i, j, data);

}

COUNTER_TEMPLATE(double, init)(Array_Type & Array, size_t i, size_t j)
{

    if (init_fun == nullptr)
        return 0.0;

    return init_fun(Array, i, j, data);

}

COUNTER_TEMPLATE(std::string, get_name)() const {
    return this->name;
}

COUNTER_TEMPLATE(std::string, get_description)() const {
    return this->desc;
}

COUNTER_TEMPLATE(void, set_name)(std::string new_name) {
    name = new_name;
}

COUNTER_TEMPLATE(void, set_description)(std::string new_desc) {
    desc = new_desc;
}

COUNTER_TEMPLATE(void, set_hasher)(Hasher_fun_type<Array_Type,Data_Type> fun) {
    hasher_fun = fun;
}

COUNTER_TEMPLATE(void, print)() const {

    printf_barry("Counter:\n");
    printf_barry("  Name       : %s\n", this->get_name().c_str());
    printf_barry("  Description: %s\n", this->get_description().c_str());

    return;

}

#define TMP_HASHER_CALL Hasher_fun_type<Array_Type,Data_Type>
COUNTER_TEMPLATE(TMP_HASHER_CALL, get_hasher)() {
    return hasher_fun;
}
#undef TMP_HASHER_CALL

////////////////////////////////////////////////////////////////////////////////
// Counters
////////////////////////////////////////////////////////////////////////////////

#define COUNTERS_TYPE() Counters<Array_Type,Data_Type>

#define COUNTERS_TEMPLATE_ARGS() <typename Array_Type, typename Data_Type>

#define COUNTERS_TEMPLATE(a,b) \
    template COUNTERS_TEMPLATE_ARGS() inline a COUNTERS_TYPE()::b

COUNTERS_TEMPLATE(, Counters)() : data(0u), hasher(nullptr) {}

COUNTERS_TEMPLATE(COUNTER_TYPE() &, operator[])(size_t idx) {

    return data[idx];

}

COUNTERS_TEMPLATE(, Counters)(const Counters<Array_Type,Data_Type> & counter_) :
    data(counter_.data), hasher(counter_.hasher) {}

COUNTERS_TEMPLATE(, Counters)(Counters<Array_Type,Data_Type> && counters_) noexcept :
    data(std::move(counters_.data)), hasher(std::move(counters_.hasher)) {}

COUNTERS_TEMPLATE(COUNTERS_TYPE(), operator=)(const Counters<Array_Type,Data_Type> & counter_) {

    if (this != &counter_)
    {
        data = counter_.data;
        hasher = counter_.hasher;
    }

    return *this;

}

COUNTERS_TEMPLATE(COUNTERS_TYPE() &, operator=)(Counters<Array_Type,Data_Type> && counters_) noexcept 
{

    if (this != &counters_) {
        data = std::move(counters_.data);
        hasher = std::move(counters_.hasher);
    }

    return *this;

}

COUNTERS_TEMPLATE(void, add_counter)(Counter<Array_Type, Data_Type> counter)
{
    
    data.push_back(counter);
    
    return;
}

COUNTERS_TEMPLATE(void, add_counter)(
    Counter_fun_type<Array_Type,Data_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Type> init_fun_,
    Hasher_fun_type<Array_Type,Data_Type>  hasher_fun_,
    Data_Type                              data_,
    std::string                            name_,
    std::string                            desc_
)
{
  
    data.emplace_back(Counter<Array_Type,Data_Type>(
        count_fun_,
        init_fun_,
        hasher_fun_,
        data_,
        name_,
        desc_
    ));
  
    return;
    
}

COUNTERS_TEMPLATE(std::vector<std::string>, get_names)() const
{

    std::vector< std::string > out;
    out.reserve(this->size());
    for (size_t i = 0u; i < this->size(); ++i)
        out.push_back(this->data.at(i).get_name());

    return out;

}

COUNTERS_TEMPLATE(std::vector<std::string>, get_descriptions)() const
{
    
    std::vector< std::string > out;
    out.reserve(this->size());
    for (size_t i = 0u; i < this->size(); ++i)
        out.push_back(data.at(i).get_description());

    return out;

}

COUNTERS_TEMPLATE(std::vector<double>, gen_hash)(
    const Array_Type & array,
    bool add_dims
)
{
    std::vector<double> res;
    
    // Iterating over the counters
    for (auto & c: data)
    {

        // If there's a hasher function, then use it!
        if (c.get_hasher())
        {

            for (auto v: c.get_hasher()(array, &(c.data)))
                res.push_back(v);

        }

    }

    // Do we need to add the dims?
    if (add_dims)
    {
        res.push_back(array.nrow());
        res.push_back(array.ncol());
    }

    // Ading the global hasher, if one exists
    if (hasher)
    {
        for (auto i: hasher(array, nullptr))
            res.push_back(i);
    }

    // We have to return something...
    if (res.size() == 0u)
        res.push_back(0.0);

    return res;

}

COUNTERS_TEMPLATE(void, add_hash)(
    Hasher_fun_type<Array_Type,Data_Type> fun_
) {

    hasher = fun_;

}

COUNTERS_TEMPLATE(void, print)(
    size_t max_length_name,
    size_t max_length_desc
) const {

    // Iterating through the counters to see the maximum name length
    size_t max_name_length = 0;
    for (const auto & c : data)
    {
        max_name_length = std::max(max_name_length, c.get_name().size());
    }

    max_name_length = std::min(max_name_length, max_length_name);

    // Figuring out the format string so it looks nice
    char fmt[100];
    snprintf(fmt, sizeof(fmt), "  - %%-%zus : %%s\n", max_name_length);

    printf_barry("Counters (%zu):\n", this->size());
    for (size_t i = 0u; i < this->size(); ++i)
    {
        // Figuring out the string to print (if needs truncation)
        auto name_to_print = data.at(i).get_name();
        if (name_to_print.size() > max_name_length)
            name_to_print = name_to_print.substr(0, max_name_length - 3) + "...";
        auto desc_to_print = data.at(i).get_description();
        if (desc_to_print.size() > max_length_desc)
            desc_to_print = desc_to_print.substr(0, max_length_desc - 3) + "...";

        auto c = data.at(i);
        printf_barry(
            fmt,
            // i,
            name_to_print.c_str(),
            desc_to_print.c_str()
        );
    }
}

#undef COUNTER_TYPE
#undef COUNTER_TEMPLATE_ARGS
#undef COUNTER_TEMPLATE
#undef COUNTERS_TYPE
#undef COUNTERS_TEMPLATE_ARGS
#undef COUNTERS_TEMPLATE

#endif 
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/counters-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/statscounter-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_STATSCOUNTER_BONES_HPP 
#define BARRY_STATSCOUNTER_BONES_HPP 1

class NetworkDense;
class NetCounterData;

/**
 * @brief Count stats for a single Array.
 * 
 * Users can a list of functions that can be used with this. The baseline set of
 * arguments is a pointer to a binary array and a dataset to add the counts to.
 */ 
template <typename Array_Type, typename Data_Type>
class StatsCounter {

private:

    // Should receive an array
    const Array_Type *               Array;
    Array_Type                       EmptyArray;
    std::vector< double >            current_stats;
      
    // We will save the data here
    Counters<Array_Type,Data_Type> * counters;
    bool                             counter_deleted  = false;

    std::vector< double > count_all_dense();
    std::vector< double > count_all_sparse();

public:
        
    /**
     * @brief Creator of a `StatsCounter`
     * 
     * @param Array_ A const pointer to a `BArray`.
     */
    StatsCounter(const Array_Type * Array_) :
        Array(Array_), EmptyArray(*Array_),
        counters(new Counters<Array_Type,Data_Type>()) {
        
        // We are removing the entries without freeing the memory. This should
        // make the insertion faster.
        EmptyArray.clear(false);
        
        return;
    }

    /**
     * @brief Copy constructor
     * 
     * @param counter 
     */
    StatsCounter(const StatsCounter<Array_Type,Data_Type> & counter);
    
    /**
     * @brief Can be created without setting the array.
     * 
     */
    StatsCounter() : Array(nullptr), EmptyArray(0u,0u),
        counters(new Counters<Array_Type,Data_Type>()) {};
    ~StatsCounter();
    
    /**
     * @brief Changes the reference array for the counting.
     * 
     * @param Array_ A pointer to an array of class `Array_Type`.
     */
    void reset_array(const Array_Type * Array_);
    
    void add_counter(Counter<Array_Type,Data_Type> f_);
    void set_counters(Counters<Array_Type,Data_Type> * counters_);
    
    /**
     * @brief Counter functions
     * This function recurses through the entries of `Array` and at each step of
     * adding a new cell it uses the functions to list the statistics.
     */
    void count_init(size_t i, size_t j);
    void count_current(size_t i, size_t j);
    std::vector< double > count_all();

    Counters<Array_Type,Data_Type> * get_counters();
    std::vector< std::string > get_names() const;
    std::vector< std::string > get_descriptions() const;

    size_t size() const {return counters->size();};
    
};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/statscounter-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/statscounter-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_STATSCOUNTER_MEAT_HPP
#define BARRY_STATSCOUNTER_MEAT_HPP 1

#define STATSCOUNTER_TYPE() StatsCounter<Array_Type,Data_Type>

#define STATSCOUNTER_TEMPLATE_ARGS() <typename Array_Type, typename Data_Type>

#define STATSCOUNTER_TEMPLATE(a,b) \
    template STATSCOUNTER_TEMPLATE_ARGS() inline a STATSCOUNTER_TYPE()::b

STATSCOUNTER_TEMPLATE(,StatsCounter)(
    const StatsCounter<Array_Type,Data_Type> & counter
)
{

    Array      = counter.Array;
    EmptyArray = *Array;
    EmptyArray.clear();
    current_stats = counter.current_stats;
      
    // We will save the data here
    counters = new Counters<Array_Type,Data_Type>((*counter.counters));
    counter_deleted  = false;

}

STATSCOUNTER_TEMPLATE(,~StatsCounter)()
{
    if (!counter_deleted)
        delete counters;
    return;
}

STATSCOUNTER_TEMPLATE(void, reset_array)(const Array_Type * Array_)
{
    
    Array      = Array_;
    EmptyArray = *Array_;
    EmptyArray.clear();
    
    return;
}

STATSCOUNTER_TEMPLATE(void, add_counter)(Counter<Array_Type,Data_Type> f_)
{
    
    counters->add_counter(f_);
    
    return;
    
}

STATSCOUNTER_TEMPLATE(void, set_counters)(Counters<Array_Type,Data_Type> * counters_)
{
    
    // Cleaning up before replacing the memory
    if (!counter_deleted)
        delete counters;
    counter_deleted = true;
    counters = counters_;
    
    return;
    
}

STATSCOUNTER_TEMPLATE(void, count_init)(size_t i,size_t j)
{
    
    // Do we have any counter?
    if (counters->size() == 0u)
        throw std::logic_error("No counters added: Cannot count without knowning what to count!");
    
    // Iterating through the functions, and updating the set of
    // statistics.
    current_stats.resize(counters->size(), 0.0);
    // change_stats.resize(counters->size(), 0.0);
    for (size_t n = 0u; n < counters->size(); ++n) 
        current_stats[n] = counters->operator[](n).init(EmptyArray, i, j);
    
    return;
}

STATSCOUNTER_TEMPLATE(void, count_current)(size_t i, size_t j)
{
    
    // Iterating through the functions, and updating the set of
    // statistics.
    for (size_t n = 0u; n < counters->size(); ++n) {
        // change_stats[n]   = counters->operator[](n).count(EmptyArray, i, j);
        // current_stats[n] += change_stats[n];
        current_stats[n] += counters->operator[](n).count(EmptyArray, i, j);
    }

    return;
    
}

template<typename Array_Type, typename Data_Type>
inline std::vector< double > StatsCounter<Array_Type,Data_Type>::count_all()
{

    if (Array->is_dense())
    {
        return count_all_dense(); 
    }
    else
    {
        return count_all_sparse();
    }

}

template<typename Array_Type, typename Data_Type>
inline std::vector< double > StatsCounter<Array_Type,Data_Type>::count_all_sparse()
{
    
    // Initializing the counter on the empty array
    count_init(0u, 0u);
    
    // Setting it to zero.
    EmptyArray.clear(false);

    #ifdef BARRY_DEBUG_LEVEL
        #if BARRY_DEBUG_LEVEL > 0
            BARRY_DEBUG_MSG("Initializing -count_all- debug. get_names():")
            BARRY_DEBUG_VEC_PRINT<std::string>(this->get_names());
        #endif
    #endif
    
    // Start iterating through the data
    for (size_t i = 0; i < Array->nrow(); ++i)
    {
        
        const auto & row = Array->row(i, false);

        // Any element?
        if (row.size() == 0u)
            continue;
        
        // If there's one, then update the statistic, by iterating
        for (const auto& col: row)
        {

            // We only insert if it is different from zero
            if (static_cast<int>(col.second.value) == 0)
                continue;
            
            // Adding a cell
            EmptyArray.insert_cell(i, col.first, col.second, false, false);

            #ifdef BARRY_DEBUG_LEVEL
                #if (BARRY_DEBUG_LEVEL >= 1)
                    BARRY_DEBUG_MSG("================================================================================")
                    BARRY_DEBUG_MSG("Debugging Stats counter: current_stats (before)")
                    std::string tmpmgs = "Inserting cell (" +
                        std::to_string(i) + ", " + std::to_string(col.first) + ")";
                    BARRY_DEBUG_MSG(tmpmgs.c_str());
                    BARRY_DEBUG_VEC_PRINT(current_stats);
                    #if (BARRY_DEBUG_LEVEL >= 2)
                        BARRY_DEBUG_MSG("Debugging Stats counter: EmptyArray")
                        EmptyArray.print();
                    #endif
                #endif
            #endif 

            // Computing the change statistics
            count_current(i, col.first);
            #ifdef BARRY_DEBUG_LEVEL
                #if (BARRY_DEBUG_LEVEL >= 1)
                    BARRY_DEBUG_MSG("Debugging Stats counter: current_stats (after)")
                    BARRY_DEBUG_VEC_PRINT(current_stats);
                #endif
            #endif
          
        } 
        
    }
    
    // Adding to the sufficient statistics
    return current_stats;
    
}

template<typename Array_Type, typename Data_Type>
inline std::vector< double > StatsCounter<Array_Type,Data_Type>::count_all_dense()
{
    
    // Initializing the counter on the empty array
    count_init(0u, 0u);
    
    // Setting it to zero.
    EmptyArray.clear(false);

    #ifdef BARRY_DEBUG_LEVEL
        #if BARRY_DEBUG_LEVEL > 0
            BARRY_DEBUG_MSG("Initializing -count_all- debug. get_names():")
            BARRY_DEBUG_VEC_PRINT<std::string>(this->get_names());
        #endif
    #endif
    
    // Start iterating through the data
    for (size_t i = 0u; i < Array->nrow(); ++i)
    {

        for (size_t j = 0u; j < Array->ncol(); ++j)
        {
            // We only insert if it is different from zero
            if (Array->is_empty(i,j))
                continue;
            
            // Adding a cell
            EmptyArray.insert_cell(i, j, 1, false, false);

            #ifdef BARRY_DEBUG_LEVEL
                #if (BARRY_DEBUG_LEVEL >= 1)
                    BARRY_DEBUG_MSG("================================================================================")
                    BARRY_DEBUG_MSG("Debugging Stats counter: current_stats (before)")
                    std::string tmpmgs = "Inserting cell (" +
                        std::to_string(i) + ", " + std::to_string(col.first) + ")";
                    BARRY_DEBUG_MSG(tmpmgs.c_str());
                    BARRY_DEBUG_VEC_PRINT(current_stats);
                    #if (BARRY_DEBUG_LEVEL >= 2)
                        BARRY_DEBUG_MSG("Debugging Stats counter: EmptyArray")
                        EmptyArray.print();
                    #endif
                #endif
            #endif 

            // Computing the change statistics
            count_current(i, j);
            #ifdef BARRY_DEBUG_LEVEL
                #if (BARRY_DEBUG_LEVEL >= 1)
                    BARRY_DEBUG_MSG("Debugging Stats counter: current_stats (after)")
                    BARRY_DEBUG_VEC_PRINT(current_stats);
                #endif
            #endif
        }
        
    }
    
    // Adding to the sufficient statistics
    return current_stats;
    
}

template STATSCOUNTER_TEMPLATE_ARGS()
inline Counters<Array_Type,Data_Type> * STATSCOUNTER_TYPE()::get_counters() {
    return this->counters;
}

STATSCOUNTER_TEMPLATE(std::vector< std::string >, get_names)() const
{
    return this->counters->get_names();
}

STATSCOUNTER_TEMPLATE(std::vector< std::string >, get_descriptions)() const
{
    return this->counters->get_descriptions();
}

#undef STATSCOUNTER_TYPE
#undef STATSCOUNTER_TEMPLATE_ARGS
#undef STATSCOUNTER_TEMPLATE

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/statscounter-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/support-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_SUPPORT_BONES_HPP 
#define BARRY_SUPPORT_BONES_HPP 1

template <typename Cell_Type, typename Data_Type>
class BArray;

template <typename Tdat>
class FreqTable;

template <typename Array_Type, typename Data_Counter_Type>
class Counters;

template <typename Array_Type, typename Data_Rule_Type>
class Rules;

template<typename Array_Type, typename Data_Type>
class Rule;

/**
 * @brief Compute the support of sufficient statistics
 * 
 * Given an array and a set of counters, this object iterates throughout the
 * support set of the Array while at the same time computing the support of
 * the sufficient statitics.
 * 
 * The members `rule` and `rule_dyn` allow constraining the support. The first
 * will establish which cells of the array will be used to iterate, for example,
 * in the case of social networks, self-loops are not allowed, so the entire
 * diagonal would be fixed to zero, reducing the size of the support.
 * 
 * In the case of `rule_dyn`, the function will stablish dynamically whether
 * the current state will be included in the counts or not. For example, this
 * set of rules can be used to constrain the support to networks that have a
 * prescribed degree sequence. 
 */ 
template <
    typename Array_Type         = BArray<bool, bool>,
    typename Data_Counter_Type  = bool,
    typename Data_Rule_Type     = bool,
    typename Data_Rule_Dyn_Type = bool 
    >
class Support {
    
private:
    void calc_backend_sparse(
        size_t pos = 0u,
        std::vector< Array_Type > * array_bank = nullptr,
        std::vector< double > * stats_bank = nullptr
    );

    void calc_backend_dense(
        size_t pos = 0u,
        std::vector< Array_Type > * array_bank = nullptr,
        std::vector< double > * stats_bank = nullptr
    );

    /**
     * @brief Reference array to generate the support.
     */
    Array_Type                               EmptyArray; ///< Temp array used to iterate through the support.
    FreqTable<>                              data;       ///< Table with the support.
    Counters<Array_Type,Data_Counter_Type> * counters;   ///< Vector of couter functions.
    Rules<Array_Type,Data_Rule_Type> *       rules;      ///< Vector of static rules (cells to iterate).
    Rules<Array_Type,Data_Rule_Dyn_Type> *   rules_dyn;  ///< Vector of dynamic rules (to include/exclude a realizaton).
    
    size_t iter_counter = 0u; ///< Number of iterations (internal use only).

public:
    
    size_t N, M;
    bool delete_counters  = true;
    bool delete_rules     = true;
    bool delete_rules_dyn = true;
    size_t max_num_elements = BARRY_MAX_NUM_ELEMENTS;
    
    // Temp variables to reduce memory allocation
    std::vector< double > current_stats;
    std::vector< size_t > coordinates_free;
    std::vector< size_t > coordinates_locked;
    size_t coordiantes_n_free;
    size_t coordiantes_n_locked;
    std::vector< double > change_stats;
    std::vector< size_t > hashes;
    std::vector< bool   > hashes_initialized;
    size_t n_counters;
    
    /**@brief Constructor passing a reference Array.
      */
    Support(const Array_Type & Array_) :
        EmptyArray(Array_),
        counters(new Counters<Array_Type,Data_Counter_Type>()),
        rules(new Rules<Array_Type,Data_Rule_Type>()),
        rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
        N(Array_.nrow()), M(Array_.ncol()), current_stats() {};
    
    /**@brief Constructor specifying the dimensions of the array (empty).
      */
    Support(size_t N_, size_t M_) :
        EmptyArray(N_, M_),
        counters(new Counters<Array_Type,Data_Counter_Type>()),
        rules(new Rules<Array_Type,Data_Rule_Type>()),
        rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
        N(N_), M(M_), current_stats() {};
    
    Support() :
        EmptyArray(0u, 0u),
        counters(new Counters<Array_Type,Data_Counter_Type>()),
        rules(new Rules<Array_Type,Data_Rule_Type>()),
        rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
        N(0u), M(0u), current_stats() {};
    
    ~Support() {
        
        if (delete_counters)
            delete counters;
        if (delete_rules)
            delete rules;
        if (delete_rules_dyn)
            delete rules_dyn;

    };
    
    void init_support(
        std::vector< Array_Type > * array_bank = nullptr,
        std::vector< double > * stats_bank = nullptr
    );
    
    /**
     * @name Resets the support calculator
     * 
     * If needed, the counters of a support object can be reused.
     * 
     * @param Array_ New array over which the support will be computed.
     */
    ///@{
    void reset_array();
    void reset_array(const Array_Type & Array_);
    ///@}
    
    /**
     * @name Manage counters 
     * 
     * @param f_ A counter to be added.
     * @param counters_ A vector of counters to be added.
     */
    ///@{
    void add_counter(Counter<Array_Type,Data_Counter_Type> f_);
    void set_counters(Counters<Array_Type,Data_Counter_Type> * counters_);
    ///@}
    
    /**
     * @name Manage rules 
     * 
     * @param f_ A rule to be added.
     * @param counters_ A vector of rules to be added.
     */
    void add_rule(Rule<Array_Type, Data_Rule_Type> * f_);
    void add_rule(Rule<Array_Type,Data_Rule_Type> f_);
    void set_rules(Rules<Array_Type,Data_Rule_Type> * rules_);
    void add_rule_dyn(Rule<Array_Type, Data_Rule_Dyn_Type> * f_);
    void add_rule_dyn(Rule<Array_Type,Data_Rule_Dyn_Type> f_);
    void set_rules_dyn(Rules<Array_Type,Data_Rule_Dyn_Type> * rules_);
    bool eval_rules_dyn(const std::vector<double> & counts, const size_t & i, const size_t & j);
    // bool eval_rules_dyn(const double * counts, const size_t & i, const size_t & j);
    ///@}

    /**
     * @brief Computes the entire support
     * 
     * Not to be used by the user. Sets the starting point in the array
     * (column-major).
     *  
     * @param array_bank If specified, the counter will add to the vector each 
     * possible state of the array, as it counts.
     * 
     * @param stats_bank If specified, the counter will add to the vector each
     * possible set of statistics, as it counts.
     * 
     */
    void calc(
        std::vector< Array_Type > * array_bank = nullptr,
        std::vector< double > * stats_bank = nullptr,
        size_t max_num_elements_ = 0u
    );
    
    const std::vector< double > & get_counts() const;
    std::vector< double > * get_current_stats(); ///< List current statistics.
    void print() const;
    
    const FreqTable< double > &              get_data() const;
    Counters<Array_Type,Data_Counter_Type> * get_counters();   ///< Vector of couter functions.
    Rules<Array_Type,Data_Rule_Type> *       get_rules();      ///< Vector of static rules (cells to iterate).
    Rules<Array_Type,Data_Rule_Dyn_Type> *   get_rules_dyn();  ///< Vector of dynamic rules (to include/exclude a realizaton).
    
};


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/support-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/support-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_SUPPORT_MEAT
#define BARRY_SUPPORT_MEAT_HPP 1

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::init_support(
    std::vector< Array_Type > * array_bank,
    std::vector< double > * stats_bank
) {

    // Resetting the counter
    this->iter_counter = 0u;
    
    // Computing the locations
    coordinates_free.clear();
    coordinates_locked.clear();
    rules->get_seq(EmptyArray, &coordinates_free, &coordinates_locked);

    coordiantes_n_free   = coordinates_free.size() / 2u;
    coordiantes_n_locked = coordinates_locked.size() / 2u;
    n_counters           = counters->size();

    hashes.resize(coordiantes_n_free, 0u);
    hashes_initialized.resize(coordiantes_n_free, false);
    
    // Computing initial statistics
    if (EmptyArray.nnozero() > 0u)
    {

        for (size_t i = 0u; i < coordiantes_n_free; ++i)
            EmptyArray.rm_cell(
                coordinates_free[i * 2u],
                coordinates_free[i * 2u + 1u],
                false, true
                );
                
    }

    // Looked coordinates should still be removed if these are
    // equivalent to zero
    for (size_t i = 0u; i < coordiantes_n_locked; ++i)
    {

        if (static_cast<int>(EmptyArray(
            coordinates_locked[i * 2u], coordinates_locked[i * 2u + 1u]
            )) == 0)

            EmptyArray.rm_cell(
                coordinates_locked[i * 2u],
                coordinates_locked[i * 2u + 1u],
                false, true
                );

    }

    // Do we have any counter?
    if (n_counters == 0u)
        throw std::logic_error("No counters added: Cannot compute the support without knowning what to count!");

    // Initial count (including constrains)
    if (coordiantes_n_locked)
    {

        StatsCounter<Array_Type,Data_Counter_Type> tmpcount(&EmptyArray);
        tmpcount.set_counters(counters);
        current_stats = tmpcount.count_all();

    }
    else
    {

        current_stats.resize(n_counters, 0.0);

        // Initialize counters
        for (size_t n = 0u; n < n_counters; ++n)
        {

            current_stats[n] = counters->operator[](n).init(
                EmptyArray,
                coordinates_free[0u],
                coordinates_free[1u]
                );

        }

    }

    // Resizing support
    data.reserve(
        pow(2.0, static_cast<double>(coordiantes_n_free)),
        counters->size()
        ); 

    // Adding to the overall count
    bool include_it = rules_dyn->operator()(EmptyArray, 0u, 0u);
    if (include_it)
        data.add(current_stats, nullptr);

    change_stats.resize(coordiantes_n_free * n_counters, 0.0);
        
    if (include_it && (array_bank != nullptr)) 
        array_bank->push_back(EmptyArray);
    
    if (include_it && (stats_bank != nullptr))
        std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

    return;

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::reset_array() {
    
    data.clear();
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::reset_array(const Array_Type & Array_) {
    
    data.clear();
    EmptyArray = Array_;
    N = Array_.nrow();
    M = Array_.ncol();
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::calc_backend_sparse(
        size_t pos,
        std::vector< Array_Type > * array_bank,
        std::vector< double > * stats_bank
    ) {
    
    #ifdef BARRY_USER_INTERRUPT
    if (++iter_counter % 1000u == 0u)
    {
        BARRY_USER_INTERRUPT
    }
    #endif

    // Did we reached the end??
    if (pos >= coordiantes_n_free)
        return;
            
    // We will pass it to the next step, if the iteration makes sense.
    calc_backend_sparse(pos + 1u, array_bank, stats_bank);
    
    // Once we have returned, everything will be back as it used to be, so we
    // treat the data as if nothing has changed.
    const size_t & coord_i = coordinates_free[pos * 2u];
    const size_t & coord_j = coordinates_free[pos * 2u + 1u];

    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray.insert_cell(
        coord_i,
        coord_j,
        EmptyArray.default_val().value,
        false, false
        );

    // Counting
    // std::vector< double > change_stats(counters.size());
    double tmp_chng;
    size_t change_stats_different = hashes_initialized[pos] ? 0u : 1u;
    for (size_t n = 0u; n < n_counters; ++n)
    {

        tmp_chng = counters->operator[](n).count(
            EmptyArray,
            coord_i,
            coord_j
            );
        
        if ((tmp_chng < DBL_MIN) & (tmp_chng > -DBL_MIN))
        {

            change_stats[pos * n_counters + n] = 0.0;

        }
        else
        {

            change_stats_different++;
            current_stats[n] += tmp_chng;
            change_stats[pos * n_counters + n] = tmp_chng;

        }

    }
    
    // Adding to the overall count
    BARRY_CHECK_SUPPORT(data, max_num_elements)
    if (rules_dyn->size() > 0u)
    {
        
        if (rules_dyn->operator()(
            EmptyArray,
            coord_i,
            coord_j
            ))
        {

            if (change_stats_different > 0u)
                hashes[pos] = data.add(current_stats, nullptr);
            else
                (void) data.add(current_stats, &hashes[pos]);

            // Need to save?
            if (array_bank != nullptr)
                array_bank->push_back(EmptyArray);
            
            if (stats_bank != nullptr)
                std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

        }
            

    } else {

        if (change_stats_different > 0u)
            hashes[pos] = data.add(current_stats, nullptr);
        else
            (void) data.add(current_stats, &hashes[pos]);

        // Need to save?
        if (array_bank != nullptr)
            array_bank->push_back(EmptyArray);
        
        if (stats_bank != nullptr)
            std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

    }
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc_backend_sparse(pos + 1u, array_bank, stats_bank);
    
    // We need to restore the state of the cell
    EmptyArray.rm_cell(
        coord_i,
        coord_j,
        false, false
        );
    
    if (change_stats_different > 0u)
    {
        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp simd
        #endif
        for (size_t n = 0u; n < n_counters; ++n) 
            current_stats[n] -= change_stats[pos * n_counters + n];
    }
        
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::calc_backend_dense(
        size_t pos,
        std::vector< Array_Type > * array_bank,
        std::vector< double > * stats_bank
    ) {

    #ifdef BARRY_USER_INTERRUPT
    if (++iter_counter % 1000u == 0u)
    {
        BARRY_USER_INTERRUPT
    }
    #endif
    
    // Did we reached the end??
    if (pos >= coordiantes_n_free)
        return;
            
    // We will pass it to the next step, if the iteration makes sense.
    calc_backend_dense(pos + 1u, array_bank, stats_bank);
    
    // Once we have returned, everything will be back as it used to be, so we
    // treat the data as if nothing has changed.
    const size_t & coord_i = coordinates_free[pos * 2u];
    const size_t & coord_j = coordinates_free[pos * 2u + 1u];

    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray.insert_cell(coord_i, coord_j, 1, false, false);

    // Counting
    // std::vector< double > change_stats(counters.size());
    double tmp_chng;
    size_t change_stats_different = hashes_initialized[pos] ? 0u : 1u;
    for (size_t n = 0u; n < n_counters; ++n)
    {

        tmp_chng = counters->operator[](n).count(
            EmptyArray,
            coord_i,
            coord_j
            );

        if ((tmp_chng < DBL_MIN) & (tmp_chng > -DBL_MIN))
        {

            change_stats[pos * n_counters + n] = 0.0;

        }
        else
        {
            if (std::isnan(tmp_chng))
                throw std::domain_error("Undefined number.");

            change_stats_different++;
            current_stats[n] += tmp_chng;
            change_stats[pos * n_counters + n] = tmp_chng;

        }

    }
    
    // Adding to the overall count
    BARRY_CHECK_SUPPORT(data, max_num_elements)
    if (rules_dyn->size() > 0u)
    {
        
        if (rules_dyn->operator()(EmptyArray, coord_i, coord_j))
        {

            if (change_stats_different > 0u)
                hashes[pos] = data.add(current_stats, nullptr);
            else
                (void) data.add(current_stats, &hashes[pos]);

            // Need to save?
            if (array_bank != nullptr)
                array_bank->push_back(EmptyArray);
            
            if (stats_bank != nullptr)
                std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

        }
            

    }
    else
    {

        if (change_stats_different > 0u)
            hashes[pos] = data.add(current_stats, nullptr);
        else
            (void) data.add(current_stats, &hashes[pos]);

        // Need to save?
        if (array_bank != nullptr)
            array_bank->push_back(EmptyArray);
        
        if (stats_bank != nullptr)
            std::copy(current_stats.begin(), current_stats.end(), std::back_inserter(*stats_bank));

    }
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc_backend_dense(pos + 1u, array_bank, stats_bank);
    
    // We need to restore the state of the cell
    EmptyArray.rm_cell(coord_i, coord_j, false, false);
    
    if (change_stats_different > 0u)
    {
        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp simd
        #endif
        for (size_t n = 0u; n < n_counters; ++n) 
            current_stats[n] -= change_stats[pos * n_counters + n];
    }
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void
Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::calc(
    std::vector< Array_Type > * array_bank,
    std::vector< double > * stats_bank,
    size_t max_num_elements_
) {

    if (max_num_elements_ != 0u)
        this->max_num_elements = max_num_elements_;

    // Generating sequence
    this->init_support(array_bank, stats_bank);

    // Recursive function to count
    if (EmptyArray.is_dense())
        calc_backend_dense(0u, array_bank, stats_bank);
    else
        calc_backend_sparse(0u, array_bank, stats_bank);

    change_stats.clear();

    if (max_num_elements_ != 0u)
        this->max_num_elements = BARRY_MAX_NUM_ELEMENTS;

    if (this->data.size() == 0u)
    {
        throw std::logic_error("The array has support of size 0 (i.e., empty support). This could be a problem in the rules (constraints).\n");
    }


    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::add_counter(
        Counter<Array_Type,Data_Counter_Type> f_
) {
    
    counters->add_counter(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::set_counters(
        Counters<Array_Type,Data_Counter_Type> * counters_
) {
    
    // Cleaning up before replacing the memory
    if (delete_counters)
        delete counters;
    delete_counters = false;
    counters = counters_;
    
    return;
    
}

/////////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::add_rule(
        Rule<Array_Type, Data_Rule_Type> * f_
) {
    
    rules->add_rule(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::add_rule(
        Rule<Array_Type,Data_Rule_Type> f_
) {
    
    rules->add_rule(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::set_rules(
        Rules<Array_Type,Data_Rule_Type> * rules_
) {
    
    // Cleaning up before replacing the memory
    if (delete_rules)
        delete rules;
    delete_rules = false;
    rules = rules_;
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::add_rule_dyn(
        Rule<Array_Type, Data_Rule_Dyn_Type> * f_
) {
    
    rules_dyn->add_rule(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::add_rule_dyn(
        Rule<Array_Type,Data_Rule_Dyn_Type> f_
) {
    
    rules_dyn->add_rule(f_);
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::set_rules_dyn(
        Rules<Array_Type,Data_Rule_Dyn_Type> * rules_
) {
    
    // Cleaning up before replacing the memory
    if (delete_rules_dyn)
        delete rules_dyn;
    delete_rules_dyn = false;
    rules_dyn = rules_;
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline bool Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::eval_rules_dyn(
    const std::vector< double > & counts,
    const size_t & i,
    const size_t & j
) {

    if (rules_dyn->size() == 0u)
        return true;

    // Swapping pointers for a while
    std::vector< double > tmpstats = current_stats;
    current_stats = counts;

    bool rule_res = rules_dyn->operator()(EmptyArray, i, j);
    current_stats = tmpstats;

    return rule_res;

}

// template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
//inline bool Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::eval_rules_dyn(
//     const double * counts,
//     const size_t & i,
//     const size_t & j
// ) {

//     if (rules_dyn->size() == 0u)
//         return true;

//     // Swapping pointers for a while
//     std::vector< double > tmpstats = current_stats;
//     current_stats = counts;

//     bool rule_res = rules_dyn->operator()(EmptyArray, i, j);
//     current_stats = tmpstats;

//     return rule_res;

// }

//////////////////////////
template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const std::vector< double > & Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::get_counts() const {
    
    return data.get_data(); 
    
}

// template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
// inline const MapVec_type<> * Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::get_counts_ptr() const {
    
//     return data.get_data_ptr();
      
// }

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< double > * Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::get_current_stats() {
    return &this->current_stats;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::print() const {

    // Starting from the name of the stats
    printf_barry("Position of variables:\n");
    for (size_t i = 0u; i < n_counters; ++i)
    {
        printf_barry(
            "[% 2i] %s\n",
            static_cast<int>(i),
            counters->operator[](i).name.c_str()
        );
    }

    data.print();
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const FreqTable<double> & Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::get_data() const {
    return this->data;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Counters<Array_Type,Data_Counter_Type> * Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::get_counters() {
    return this->counters;
}   
    
template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Rules<Array_Type,Data_Rule_Type> * Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::get_rules() {
    return this->rules;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Rules<Array_Type,Data_Rule_Dyn_Type> * Support<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::get_rules_dyn() {
    return this->rules_dyn;
}


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/support-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/powerset-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_POWERSET_BONES_HPP 
#define BARRY_POWERSET_BONES_HPP 1

/**
 * @brief Powerset of a binary array
 * 
 * @tparam Array_Type 
 * @tparam Data_Rule_Type 
 */
template <typename Array_Type = BArray<>, typename Data_Rule_Type = bool> 
class PowerSet {
    
private:
    void calc_backend_sparse(size_t pos = 0u);  
    void calc_backend_dense(size_t pos = 0u);  

public:
    Array_Type                         EmptyArray;
    std::vector< Array_Type >          data;
    Rules<Array_Type,Data_Rule_Type> * rules;

    size_t N, M;
    bool rules_deleted   = false;

    // Tempvars
    std::vector< size_t >  coordinates_free;
    std::vector< size_t >  coordinates_locked;
    size_t n_free;
    size_t n_locked;
    
    /**
     * @name Construct and destroy a PowerSet object
     * 
     */
    ///@{
    PowerSet() : 
    EmptyArray(), data(0u), rules(new Rules<Array_Type,Data_Rule_Type>()), N(0u), M(0u) {};
    PowerSet(size_t N_, size_t M_) :
        EmptyArray(N_, M_), data(0u),
        rules(new Rules<Array_Type,Data_Rule_Type>()), N(N_), M(M_) {};
    PowerSet(const Array_Type & array);

    ~PowerSet();
    ///@}
    
    void init_support();
    void calc();
    void reset(size_t N_, size_t M_);
    
    /**
     * @name Wrappers for the `Rules` member. 
     * @details These will add rules to the model, which are shared by the
     * support and the actual counter function.
     */
    ///@{
    void add_rule(Rule<Array_Type, Data_Rule_Type> rule);
    void add_rule(
        Rule_fun_type<Array_Type,Data_Rule_Type> count_fun_,
        Data_Rule_Type data_
    );
    ///@}
    

    /** @name Getter functions */
    ///@{
    const std::vector< Array_Type > * get_data_ptr() const {return &data;};
    std::vector< Array_Type > get_data() const {return data;};
    typename std::vector< Array_Type >::iterator begin() {return data.begin();};
    typename std::vector< Array_Type >::iterator end() {return data.end();};
    std::size_t size() const noexcept {return data.size();};
    const Array_Type& operator[](const size_t & i) const {return data.at(i);};
    ///@}
    
};

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/powerset-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/powerset-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_POWERSET_MEAT_HPP
#define BARRY_POWERSET_MEAT_HPP 1

template <typename Array_Type, typename Data_Rule_Type>
inline PowerSet<Array_Type,Data_Rule_Type>::PowerSet(
    const Array_Type & array
) : EmptyArray(array), data(0u),
        rules(new Rules<Array_Type,Data_Rule_Type>()), N(array.nrow()), M(array.ncol()) {

}

template <typename Array_Type, typename Data_Rule_Type>
inline PowerSet<Array_Type,Data_Rule_Type>::~PowerSet() {
    if (!this->rules_deleted)
        delete rules;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::init_support()
{
    
    // Computing the locations
    coordinates_free.clear();
    coordinates_locked.clear();
    rules->get_seq(EmptyArray, &coordinates_free, &coordinates_locked);

    n_free   = coordinates_free.size() / 2u;
    n_locked = coordinates_locked.size() / 2u;
    
    // Computing initial statistics
    if (EmptyArray.nnozero() > 0u)
    {

        if (EmptyArray.is_dense())
        {

            for (size_t i = 0u; i < n_free; ++i) 
                EmptyArray(
                    coordinates_free[i * 2u],
                    coordinates_free[i * 2u + 1u]
                    ) = 0;

        }
        else
        {

            for (size_t i = 0u; i < n_free; ++i) 
                EmptyArray.rm_cell(
                    coordinates_free[i * 2u],
                    coordinates_free[i * 2u + 1u],
                    false,
                    true
                );


        }
            
    }

    // EmptyArray.clear(true);
    // EmptyArray.reserve();
    
    // Resizing support
    data.reserve(pow(2.0, n_free)); 

    // Adding the empty array to the set
    data.push_back(EmptyArray);
    
    return;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type, Data_Rule_Type>::calc_backend_sparse(
    size_t pos
)
{
    
    // Did we reached the end??
    if (pos >= n_free)
        return;
            
    // We will pass it to the next step, if the iteration makes sense.
    calc_backend_sparse(pos + 1u);
        
    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray.insert_cell(
        coordinates_free[pos * 2u],
        coordinates_free[pos * 2u + 1u],
        EmptyArray.default_val().value,
        false, false
        );

    data.push_back(EmptyArray);

    #ifdef BARRY_USER_INTERRUPT
    if (data.size() % 1000u == 0u)
    {
        BARRY_USER_INTERRUPT
    }
    #endif
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc_backend_sparse(pos + 1u);
    
    // We need to restore the state of the cell
    EmptyArray.rm_cell(
        coordinates_free[pos * 2u],
        coordinates_free[pos * 2u + 1u],
        false, false
        );  
    
    return;
    
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type, Data_Rule_Type>::calc_backend_dense(
    size_t pos
)
{
    
    // Did we reached the end??
    if (pos >= n_free)
        return;
            
    // We will pass it to the next step, if the iteration makes sense.
    calc_backend_dense(pos + 1u);
        
    // Toggle the cell (we will toggle it back after calling the counter)
    EmptyArray(coordinates_free[pos * 2u], coordinates_free[pos * 2u + 1u]) = 1;

    data.push_back(EmptyArray);
    
    // Again, we only pass it to the next level iff the next level is not
    // passed the last step.
    calc_backend_dense(pos + 1u);
    
    // We need to restore the state of the cell
    EmptyArray(coordinates_free[pos * 2u], coordinates_free[pos * 2u + 1u]) = 0;
    
    return;
    
}


/***
  * Function to generate the powerset of the 
  */
template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type, Data_Rule_Type>::calc() {

    // Generating sequence
    this->init_support();

    // Recursive function to count
    if (EmptyArray.is_dense())
        calc_backend_dense(0u);
    else
        calc_backend_sparse(0u);

    return;
    
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::reset(
        size_t N_,
        size_t M_
) {
    
    data.empty();
    N = N_, M = M_;
    
    return;

}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
        Rule<Array_Type, Data_Rule_Type> rule
) {
    
    rules->add_rule(rule);
    return;
}

template <typename Array_Type, typename Data_Rule_Type>
inline void PowerSet<Array_Type,Data_Rule_Type>::add_rule(
        Rule_fun_type<Array_Type,Data_Rule_Type> rule_fun_,
        Data_Rule_Type data_
) {
    
    rules->add_rule(
        rule_fun_,
        data_
    );
    
    return;
    
}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/powerset-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/model-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_MODEL_BONES_HPP 
#define BARRY_MODEL_BONES_HPP 1

/**
 * @ingroup stat-models
 * @brief General framework for discrete exponential models.
 * This class allows generating discrete exponential models in the form of a linear
 * exponential model:
 * \f[
 * \frac{
 *    \exp{\left(\theta^{\mbox{t}}c(A)\right)}
 *  }{
 *    \sum_{A'\in \mathcal{A}}\exp{\left(\theta^{\mbox{t}}c(A')\right)}
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
     * P{Y(i,j) = | Y^C, theta}, i.e., the probability of observing the entry Y(i,j) equal
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
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/model-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/model-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_MODEL_MEAT_HPP 
#define BARRY_MODEL_MEAT_HPP 1

/**
 * @defgroup stat-models Statistical Models
 * @brief Statistical models available in `barry`.
 */

inline double update_normalizing_constant(
    const std::vector<double> & params,
    const double * support,
    size_t k,
    size_t n
)
{
    double res = 0.0;

    std::vector< double > resv(n, 0.0);

    for (size_t j = 0u; j < (k - 1u); ++j)
    {

        const double p = params[j];
        
        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp simd 
        #elif defined(__GNUC__) && !defined(__clang__)
            #pragma GCC ivdep
        #endif
        for (size_t i = 0u; i < n; ++i)
            resv[i] += (*(support + i * k + 1u + j)) * p;

    }

    // Accumulate resv to a double res        
    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp simd reduction(+:res)
    #elif defined(__GNUC__) && !defined(__clang__)
        #pragma GCC ivdep
    #endif
    for (size_t i = 0u; i < n; ++i)
    {
        res += std::exp(resv[i] BARRY_SAFE_EXP) * (*(support + i * k));
    }




    #ifdef BARRY_DEBUG
    if (std::isnan(res))
        throw std::overflow_error(
            std::string("NaN in update_normalizing_constant. ") +
            std::string("res = ") + std::to_string(res) +
            std::string(", k = ") + std::to_string(k) +
            std::string(", n = ") + std::to_string(n)
            );
    if (std::isinf(res))
        throw std::overflow_error(
            std::string("Inf in update_normalizing_constant. ") +
            std::string("res = ") + std::to_string(res) +
            std::string(", k = ") + std::to_string(k) +
            std::string(", n = ") + std::to_string(n)
            );

    #endif

    return res;
    
}

inline double likelihood_(
        const double * stats_target,
        const std::vector< double > & params,
        const double normalizing_constant,
        size_t n_params,
        bool log_ = false
) {
    
    if (n_params != params.size())
        throw std::length_error("-stats_target- and -params- should have the same length.");
        
    double numerator = 0.0;
    
    // Computing the numerator
    #ifdef __INTEL_LLVM_COMPILER
    #pragma code_align 32
    #endif
    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp simd reduction(+:numerator)
    #endif
    for (size_t j = 0u; j < n_params; ++j)
        numerator += *(stats_target + j) * params[j];

    if (!log_)
        numerator = std::exp(numerator BARRY_SAFE_EXP);
    else
        return numerator BARRY_SAFE_EXP - std::log(normalizing_constant);

    double ans = numerator/normalizing_constant;

    #ifdef BARRY_DEBUG
    if (std::isnan(ans))
        throw std::overflow_error(
            std::string("NaN in likelihood_. ") +
            std::string("numerator = ") + std::to_string(numerator) +
            std::string(", normalizing_constant = ") +
            std::to_string(normalizing_constant)
            );
    if (std::isinf(ans))
        throw std::overflow_error(
            std::string("Inf in likelihood_. ") +
            std::string("numerator = ") + std::to_string(numerator) +
            std::string(", normalizing_constant = ") +
            std::to_string(normalizing_constant)
            );

    if (ans > 1.0)
        throw std::overflow_error(
            std::string("Likelihood > 1.0") +
            std::string("numerator = ") + std::to_string(numerator) +
            std::string(", normalizing_constant = ") +
            std::to_string(normalizing_constant)
            );
    #endif

    return ans;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline void Model<Array_Type, Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::update_normalizing_constants(
    const std::vector< double > & params,
    size_t ncores,
    int i
) {

    const size_t n = stats_support_sizes.size();

    // Barrier to make sure paralelization makes sense
    if ((ncores > 1u) && (n < 128u))
        ncores = 1u;

    
    if (i >= 0)
        ncores = 1u;
    
    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp parallel for firstprivate(params) num_threads(ncores) \
        shared(n, normalizing_constants, first_calc_done, \
            stats_support, stats_support_sizes, stats_support_sizes_acc, i) \
        default(none)
    #endif
    for (size_t s = 0u; s < n; ++s)
    {

        if ((i > -1) && (i != static_cast<int>(s)))
            continue;

        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[s];

        first_calc_done[s] = true;
        normalizing_constants[s] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[s] * k
                ], k, n
        );

    }

    return;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::update_likelihoods(
    const std::vector< double > & params,
    size_t ncores
) {
    
    update_normalizing_constants(params, ncores);

    size_t n_params = params.size();

    if (stats_likelihood.size() != stats_target.size())
        stats_likelihood.resize(stats_target.size());

    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp parallel for simd num_threads(ncores) \
        shared(n_params, stats_target, normalizing_constants, arrays2support, \
            params) \
        default(none)
    #endif
    for (size_t s = 0u; s < stats_target.size(); ++s)
    {
        stats_likelihood[s] = 0.0;
        for (size_t j = 0u; j < n_params; ++j)
            stats_likelihood[s] += stats_target[s][j] * params[j];

        stats_likelihood[s] =
            std::exp(stats_likelihood[s] BARRY_SAFE_EXP)/
            normalizing_constants[arrays2support[s]];
    }
    
    return;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::update_pset_probs(
    const std::vector< double > & params,
    size_t ncores,
    int i
) {

    update_normalizing_constants(params, ncores, i);
    
    if (i > -1)
        params_last[i] = params;

    size_t n_params = params.size();
    pset_probs.resize(
        pset_locations.back() + 
        pset_sizes.back()
        );

    // No need to paralelize if there is only one core
    if (i >= 0)
       ncores = 1u; 

    #if defined(__OPENMP) || defined(_OPENMP)
    #pragma omp parallel for num_threads(ncores) collapse(1) \
        shared(n_params, pset_stats, pset_probs, normalizing_constants, pset_sizes, \
            params, i) \
        default(none)
    #endif
    for (size_t s = 0u; s < pset_sizes.size(); ++s)
    {

        if ((i >= 0) && (i != static_cast<int>(s)))
            continue;

        // When does the pset starts
        size_t pset_start = pset_locations[s];

        // Looping over observations of the pset
        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp simd 
        #endif
        for (size_t a = 0u; a < pset_sizes[s]; ++a)
        {

            // Start location in the array
            size_t start_loc = pset_start * n_params + a * n_params;
            
            pset_probs[pset_start + a] = 0.0;

            // Looping over the parameters
            for (size_t j = 0u; j < n_params; ++j)
                pset_probs[pset_start + a] +=
                    pset_stats[start_loc + j] * params[j];

            // Now turning into a probability
            pset_probs[pset_start + a] =
                std::exp(pset_probs[pset_start + a] BARRY_SAFE_EXP)/
                normalizing_constants[s];
        }

        #ifdef BARRY_DEBUG
        // Making sure the probabilities add to one
        double totprob = 0.0;
        for (size_t i_ = 0u; i_ < pset_sizes[s]; ++i)
            totprob =+ pset_probs[pset_start + i_];

        if (std::abs(totprob - 1) > 1e-6)
            throw std::runtime_error(
                std::string("Probabilities do not add to one! ") +
                std::string("totprob = ") + std::to_string(totprob)
            );

        #endif
    }
    
    return;

}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::Model() :
    stats_support(0u),
    stats_support_sizes(0u),
    stats_support_sizes_acc(0u),
    stats_support_n_arrays(0u),
    stats_target(0u),
    stats_likelihood(0u),
    arrays2support(0u),
    keys2support(0u),
    pset_arrays(0u), pset_stats(0u),
    counters(new Counters<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
    support_fun(), counter_fun(), delete_counters(true),
    delete_rules(true),
    delete_rules_dyn(true),
    transform_model_fun(nullptr),
    transform_model_term_names(0u)
{  

    // Counters are shared
    support_fun.set_counters(counters);
    counter_fun.set_counters(counters);
    
    // Rules are shared
    support_fun.set_rules(rules);
    support_fun.set_rules_dyn(rules_dyn);
    
    return;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::Model(
    size_t size_
    ) :
    stats_support(0u),
    stats_support_sizes(0u),
    stats_support_sizes_acc(0u),
    stats_support_n_arrays(0u),
    stats_target(0u),
    stats_likelihood(0u),
    arrays2support(0u), keys2support(0u), 
    pset_arrays(0u), pset_stats(0u),
    counters(new Counters<Array_Type,Data_Counter_Type>()),
    rules(new Rules<Array_Type,Data_Rule_Type>()),
    rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>()),
    support_fun(), counter_fun(), delete_counters(true),
    delete_rules(true),
    delete_rules_dyn(true),
    transform_model_fun(nullptr),
    transform_model_term_names(0u)
{
    
    stats_target.reserve(size_);
    arrays2support.reserve(size_);

    // Counters are shared
    support_fun.set_counters(counters);
    counter_fun.set_counters(counters);
    
    // Rules are shared
    support_fun.set_rules(rules);
    support_fun.set_rules_dyn(rules_dyn);
        
    return;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::Model(
    const Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type> & Model_
    ) : 
    stats_support(Model_.stats_support),
    stats_support_sizes(Model_.stats_support_sizes),
    stats_support_sizes_acc(Model_.stats_support_sizes_acc),
    stats_support_n_arrays(Model_.stats_support_n_arrays),
    stats_target(Model_.stats_target),
    stats_likelihood(Model_.stats_likelihood),
    arrays2support(Model_.arrays2support),
    keys2support(Model_.keys2support),
    pset_arrays(Model_.pset_arrays),
    pset_stats(Model_.pset_stats),
    pset_probs(Model_.pset_probs),
    pset_sizes(Model_.pset_sizes),
    pset_locations(Model_.pset_locations),
    counters(new Counters<Array_Type,Data_Counter_Type>(*(Model_.counters))),
    rules(new Rules<Array_Type,Data_Rule_Type>(*(Model_.rules))),
    rules_dyn(new Rules<Array_Type,Data_Rule_Dyn_Type>(*(Model_.rules_dyn))),
    support_fun(),
    counter_fun(),
    params_last(Model_.params_last),
    normalizing_constants(Model_.normalizing_constants),
    first_calc_done(Model_.first_calc_done),
    delete_counters(true),
    delete_rules(true),
    delete_rules_dyn(true),
    transform_model_fun(Model_.transform_model_fun),
    transform_model_term_names(Model_.transform_model_term_names)
    {
    
    // Counters are shared
    support_fun.set_counters(counters);
    counter_fun.set_counters(counters);
    
    // Rules are shared
    support_fun.set_rules(rules);
    support_fun.set_rules_dyn(rules_dyn);

    return;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type> & 
    Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::operator=(
    const Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type> & Model_
) {
    
    // Clearing
    if (this != &Model_) {

        if (delete_counters)
            delete counters;

        if (delete_rules)
            delete rules;
        
        if (delete_rules_dyn)
            delete rules_dyn;
        
        stats_support              = Model_.stats_support;
        stats_support_sizes        = Model_.stats_support_sizes;
        stats_support_sizes_acc    = Model_.stats_support_sizes_acc;
        stats_support_n_arrays     = Model_.stats_support_n_arrays;
        stats_target               = Model_.stats_target;
        stats_likelihood           = Model_.stats_likelihood;
        arrays2support             = Model_.arrays2support;
        keys2support               = Model_.keys2support;
        pset_arrays                = Model_.pset_arrays;
        pset_stats                 = Model_.pset_stats;
        pset_probs                 = Model_.pset_probs;
        pset_sizes                 = Model_.pset_sizes;
        pset_locations             = Model_.pset_locations;
        counters                   = new Counters<Array_Type,Data_Counter_Type>(*(Model_.counters));
        rules                      = new Rules<Array_Type,Data_Rule_Type>(*(Model_.rules));
        rules_dyn                  = new Rules<Array_Type,Data_Rule_Dyn_Type>(*(Model_.rules_dyn));
        delete_counters            = true;
        delete_rules               = true;
        delete_rules_dyn           = true;
        params_last                = Model_.params_last;
        normalizing_constants      = Model_.normalizing_constants;
        first_calc_done            = Model_.first_calc_done;
        transform_model_fun        = Model_.transform_model_fun;
        transform_model_term_names = Model_.transform_model_term_names;

        // Counters are shared
        support_fun.set_counters(counters);
        counter_fun.set_counters(counters);
        
        // Rules are shared
        support_fun.set_rules(rules);
        support_fun.set_rules_dyn(rules_dyn);
        
    }
        
    return *this;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: store_psets() noexcept {
    with_pset = true;
    return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< double > Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: gen_key(
    const Array_Type & Array_
) {
    return this->counters->gen_hash(Array_);   
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_counter(
        Counter<Array_Type, Data_Counter_Type> & counter
) {
    
    counters->add_counter(counter, Data_Counter_Type());
    return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_counter(
    Counter_fun_type<Array_Type,Data_Counter_Type> count_fun_,
    Counter_fun_type<Array_Type,Data_Counter_Type> init_fun_,
    Data_Counter_Type                              data_
) {
    
    counters->add_counter(
        count_fun_,
        init_fun_,
        data_
    );
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: set_counters(
    Counters<Array_Type,Data_Counter_Type> * counters_
) {

    if (delete_counters) {
        delete counters;
        delete_counters = false;
    }
    
    this->counters = counters_;
    support_fun.set_counters(counters);
    counter_fun.set_counters(counters);
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_hasher(
    Hasher_fun_type<Array_Type,Data_Counter_Type> fun_
) {

    counters->add_hash(fun_);

}

////////////////////////////////////////////////////////////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_rule(
    Rule<Array_Type, Data_Rule_Type> & rules
) {
    
    rules->add_rule(rules, Data_Rule_Type());
    return;
}


template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: set_rules(
    Rules<Array_Type,Data_Rule_Type> * rules_
) {

    if (delete_rules)
        delete rules;

    this->rules = rules_;
    this->delete_rules = false;

    support_fun.set_rules(rules);
    return;

}

////////////////////////////////////////////////////////////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_rule_dyn(
    Rule<Array_Type, Data_Rule_Dyn_Type> & rules_
) {
    
    rules_dyn->add_rule(rules_, Data_Rule_Dyn_Type());
    return;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: add_rule_dyn(
    Rule_fun_type<Array_Type,Data_Rule_Dyn_Type> rule_fun_,
    Data_Rule_Dyn_Type                           data_
) {
    
    rules_dyn->add_rule(
        rule_fun_,
        data_
    );
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: set_rules_dyn(
    Rules<Array_Type,Data_Rule_Dyn_Type> * rules_
) {

    if (delete_rules_dyn)
        delete rules_dyn;

    this->rules_dyn = rules_;
    this->delete_rules_dyn = false;
    support_fun.set_rules_dyn(rules_dyn);
    return;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::add_array(
    const Array_Type & Array_,
    bool force_new
) {
    
    // Array counts (target statistics)
    counter_fun.reset_array(&Array_);
    
    if (transform_model_fun)
    {
        
        auto tmpcounts = counter_fun.count_all();
        stats_target.emplace_back(
            transform_model_fun(&tmpcounts[0u], tmpcounts.size())
            );

    } else
        stats_target.push_back(counter_fun.count_all());
    
    // If the data hasn't been analyzed earlier, then we need to compute
    // the support
    std::vector< double > key = counters->gen_hash(Array_);
    MapVec_type< double, size_t >::const_iterator locator = keys2support.find(key);
    if (force_new | (locator == keys2support.end()))
    {

        // Current size of the support stats
        size_t stats_support_size = stats_support.size();
        
        // Adding to the map
        keys2support[key] = stats_support_sizes.size();
        stats_support_n_arrays.push_back(1u);       // How many elements now
        arrays2support.push_back(stats_support_sizes.size()); // Map of the array id to the support
        
        // Computing support using the counters included in the model
        support_fun.reset_array(Array_);
        
        /** When computing with the powerset, we need to grow the corresponding
            * vectors on the fly */
        if (with_pset)
        {
            
            // Making space for storing the support
            pset_arrays.resize(pset_arrays.size() + 1u);

            // Current size of the powerset
            size_t pset_stats_size = pset_stats.size();
            
            try
            {
                
                support_fun.calc(
                    &(pset_arrays[pset_arrays.size() - 1u]),
                    &pset_stats
                );
                
            }
            catch (const std::exception& e)
            {
                
                printf_barry(
                    "A problem ocurred while trying to add the array (and recording the powerset). "
                );
                printf_barry("with error %s\n", e.what());
                printf_barry("Here is the array that generated the error.\n");
                Array_.print();
                throw std::logic_error("");
                
            }

            // Recording the number of elements
            pset_locations.push_back(
                pset_locations.size() == 0u ?
                    0u :
                    pset_locations.back() + pset_sizes.back()
                );

            pset_sizes.push_back(
                (pset_stats.size() - pset_stats_size) / (counter_fun.size())
                );
                
                

            
        }
        else
        {
            try
            {

                support_fun.calc();
                
            }
            catch (const std::exception& e)
            {

                printf_barry(
                    "A problem ocurred while trying to add the array (and recording the powerset). "
                );
                printf_barry("with error %s\n", e.what());
                printf_barry("Here is the array that generated the error.\n");
                Array_.print();
                throw std::logic_error("");

            }
        }
        
        if (transform_model_fun)
        {
            auto tmpsupport = support_fun.get_counts();
            size_t k = counter_fun.size();
            size_t n = tmpsupport.size() / (k + 1);

            std::vector< double > s_new(0u);            
            s_new.reserve(tmpsupport.size());

            for (size_t i = 0u; i < n; ++i)
            {

                // Appending size
                s_new.push_back(tmpsupport[i * (k + 1u)]);

                // Applying transformation and adding to the new set
                auto res = transform_model_fun(&tmpsupport[i * (k + 1u) + 1u], k);
                std::copy(res.begin(), res.end(), std::back_inserter(s_new));

            }

            for (auto & s : s_new)
                stats_support.push_back(s);
            

        } else {
            for (const auto & s: support_fun.get_counts())
                stats_support.push_back(s);
        }
        
        // Making room for the previous parameters. This will be used to check if
        // the normalizing constant has been updated or not.
        params_last.push_back(stats_target[0u]);
        normalizing_constants.push_back(0.0);
        first_calc_done.push_back(false);

        // Incrementing the size of the support set
        if (stats_support_sizes.size() == 0u)
        {
            stats_support_sizes_acc.push_back(0u);    
        } else {
            stats_support_sizes_acc.push_back(
                stats_support_sizes.back() + 
                stats_support_sizes_acc.back()
            );
        }


        stats_support_sizes.push_back(
            
            (stats_support.size() - stats_support_size)/
                (counter_fun.size() + 1u)

            );
        
        return arrays2support.size() - 1u;
        
    }
    
    // Increasing the number of arrays in that stat
    ++stats_support_n_arrays[locator->second];
    
    // Adding the corresponding map
    arrays2support.push_back(locator->second);
    
    return arrays2support.size() - 1u;

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood(
    const std::vector<double> & params,
    const size_t & i,
    bool as_log,
    bool no_update_normconst
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    size_t idx = arrays2support[i];

    // Checking if this actually has a change of happening
    if (this->stats_support_sizes[idx] == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Checking if we have updated the normalizing constant or not
    if (!no_update_normconst && (!first_calc_done[idx] || !vec_equal_approx(params, params_last[idx])))
    {
        
        first_calc_done[idx] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[idx];

        normalizing_constants[idx] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[idx] * k
                ], k, n
        );
        
        params_last[idx] = params;
        
    }
    
    return likelihood_(
        &stats_target[i],
        params,
        normalizing_constants[idx],
        nterms(),
        as_log
    );
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood(
    const std::vector<double> & params,
    const Array_Type & Array_,
    int i,
    bool as_log,
    bool no_update_normconst
) {
    
    // Key of the support set to use
    int loc;

    if (i < 0)
    {

        std::vector< double > key = counters->gen_hash(Array_);
        MapVec_type< double, size_t >::const_iterator locator = keys2support.find(key);
        if (locator == keys2support.end()) 
            throw std::range_error(
                "This type of array has not been included in the model."
                );

        loc = locator->second;

    }
    else
    {

        if (static_cast<size_t>(i) >= arrays2support.size())
            throw std::range_error(
                "This type of array has not been included in the model."
                );

        loc = arrays2support[i];

    }

    // Checking if this actually has a change of happening
    if (this->stats_support_sizes[loc] == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Counting stats_target
    StatsCounter< Array_Type, Data_Counter_Type> tmpstats(&Array_);

    tmpstats.set_counters(this->counters);
    
    std::vector< double > target_ = tmpstats.count_all();

    if (transform_model_fun)
        target_ = transform_model_fun(&target_[0u], target_.size());

    // Checking if we have updated the normalizing constant or not
    if (!no_update_normconst && (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc])) )
    {
        
        first_calc_done[loc] = true;

        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[loc];
        
        normalizing_constants[loc] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[loc] * k
                ], k, n
        );
        
        params_last[loc] = params;
        
    }

    // Checking if passes the rules
    if (!support_fun.eval_rules_dyn(target_, 0u, 0u))
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    return likelihood_(
        &target_[0u],
        params,
        normalizing_constants[loc],
        nterms(),
        as_log
    );
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood(
    const std::vector<double> & params,
    const std::vector<double> & target_,
    const size_t & i,
    bool as_log,
    bool no_update_normconst
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    size_t loc = arrays2support[i];

    // Checking if passes the rules
    if (!support_fun.eval_rules_dyn(target_, 0u, 0u))
    {

        // Concatenating the elements of target_ into aa single string
        std::string target_str = "";
        for (size_t i = 0u; i < target_.size(); ++i)
            target_str += std::to_string(target_[i]) + " ";

        throw std::range_error(
            "The array is not in the support set. The array's statistics are: " +
            target_str +
            std::string(".")
            );
    }
        

    // Checking if this actually has a change of happening
    if (this->stats_support_sizes[loc] == 0u)
    {
        throw std::logic_error("The support set for this array is empty.");
    }
    
    // Checking if we have updated the normalizing constant or not
    if (!no_update_normconst && (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc])) ) {
        
        first_calc_done[loc] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[loc];

        normalizing_constants[loc] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[loc] * k
                ], k, n
        );
        
        params_last[loc] = params;
        
    }
    
    return likelihood_(
        &target_[0u],
        params,
        normalizing_constants[loc],
        nterms(),
        as_log
    );
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood(
    const std::vector<double> & params,
    const double * target_,
    const size_t & i,
    bool as_log,
    bool no_update_normconst
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    size_t loc = arrays2support[i];

    // Checking if passes the rules
    if (support_fun.get_rules_dyn()->size() > 0u)
    {

        std::vector< double > tmp_target;
        tmp_target.reserve(nterms());
        for (size_t t = 0u; t < nterms(); ++t)
            tmp_target.push_back(*(target_ + t));

        if (!support_fun.eval_rules_dyn(tmp_target, 0u, 0u))
        {
            // Concatenating the elements of target_ into aa single string
            std::string target_str = "";
            for (size_t i = 0u; i < nterms(); ++i)
                target_str += std::to_string((*target_ + i)) + " ";

            throw std::range_error(
                "The array is not in the support set. The array's statistics are: " + target_str + std::string(".")
                );
        }

    }

    // Checking if this actually has a change of happening
    if (this->stats_support_sizes[loc] == 0u)
    {
        throw std::logic_error("The support set for this array is empty.");
    }
    
    // Checking if we have updated the normalizing constant or not
    if (!no_update_normconst && (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) )) {
        
        first_calc_done[loc] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support_sizes[loc];

        normalizing_constants[loc] = update_normalizing_constant(
            params, &stats_support[
                stats_support_sizes_acc[loc] * k
            ], k, n
        );
        
        params_last[loc] = params;
        
    }
    
    return likelihood_(
        target_,
        params,
        normalizing_constants[loc],
        nterms(),
        as_log
    );
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::likelihood_total(
    const std::vector<double> & params,
    bool as_log,
    BARRY_NCORES_ARG(),
    bool no_update_normconst
) {
    
    size_t params_last_size = params_last.size();

    if (!no_update_normconst)
    {
        #if defined(__OPENMP) || defined(_OPENMP)
        #pragma omp parallel for num_threads(ncores) \
            shared(normalizing_constants, params_last, first_calc_done, \
                stats_support, stats_support_sizes, stats_support_sizes_acc) \
            firstprivate(params)
        #endif
        for (size_t i = 0u; i < params_last_size; ++i)
        {

            if (!first_calc_done[i] || !vec_equal_approx(params, params_last[i]) )
            {

                size_t k = params.size() + 1u;
                size_t n = stats_support_sizes[i];
                
                first_calc_done[i] = true;
                normalizing_constants[i] = update_normalizing_constant(
                    params, &stats_support[
                        stats_support_sizes_acc[i] * k
                    ], k, n
                );
                
                params_last[i] = params;
                
            }

        }
    }
    
    double res = 0.0;
    if (as_log)
    {

        for (size_t i = 0; i < stats_target.size(); ++i) 
            res += vec_inner_prod(
                &stats_target[i][0u],
                &params[0u],
                params.size()
                ) BARRY_SAFE_EXP;
        
        #if defined(__OPENMP) || defined(_OPENMP) 
        #pragma omp simd reduction(-:res)
        #endif
        for (size_t i = 0u; i < params_last_size; ++i)
            res -= (std::log(normalizing_constants[i]) * this->stats_support_n_arrays[i]);

    } else {
        
        res = 1.0;
        size_t stats_target_size = stats_target.size();
        #if defined(__OPENMP) || defined(_OPENMP) 
        #pragma omp simd reduction(*:res)
        #endif
        for (size_t i = 0; i < stats_target_size; ++i)
            res *= std::exp(
                vec_inner_prod(
                    &stats_target[i][0u],
                    &params[0u],
                    params.size()
                ) BARRY_SAFE_EXP) / 
                normalizing_constants[arrays2support[i]];
        
    }
    
    return res;
    
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline const std::vector< double > &
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_normalizing_constants() const {
    
    return normalizing_constants;
    
}

template<
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline const std::vector< double > &
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_likelihoods() const {
    
    return stats_likelihood;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const std::vector< Array_Type > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset(
    const size_t & i
) {

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");


    return &pset_arrays[arrays2support[i]];

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const double *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_stats(
    const size_t & i
) {

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    return &pset_stats[pset_locations[arrays2support[i]] * counter_fun.size()];

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: print_stats(size_t i) const
{
    
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    size_t k       = nterms();
    size_t nunique = stats_support_sizes.size();

    for (size_t l = 0u; l < nunique; ++l)
    {

        printf_barry("% 5i ", static_cast<int>(l));

        printf_barry("counts: %.0f motif: ", stats_support[
            stats_support_sizes_acc[l] * (k + 1u) 
            // l * (k + 1u)
            ]);
        
        for (size_t j = 0u; j < k; ++j)
        {
            printf_barry(
                "%.2f, ",
                stats_support[
                    stats_support_sizes_acc[l] * (k + 1u) + j + 1u
                    ]);
        }

        printf_barry("\n");

    }
    
    return;
    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline void Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::print() const
{

    // Relevant information:
    // - Number of arrays involved
    // - Size of the support
    // - Terms involved

    int min_v = std::numeric_limits<int>::max();
    int max_v = 0;

    for (const auto & stat : this->stats_support_sizes)
    {

        if (static_cast<int>(stat) > max_v)
            max_v = static_cast<int>(stat);
        
        if (static_cast<int>(stat) < min_v)
            min_v = static_cast<int>(stat);

    }  

    // The vectors in the support reflec the size of nterms x entries
    // max_v /= static_cast<int>(nterms() + 1);
    // min_v /= static_cast<int>(nterms() + 1);

    if (this->size() > 0u)
    {
        printf_barry(
            "Num. of Arrays       : %i\n",
            static_cast<int>(this->size())
        );
        printf_barry(
            "Support size         : %i\n",
            static_cast<int>(this->size_unique())
        );
        printf_barry("Support size range   : [%i, %i]\n", min_v, max_v);
    }
    else 
    {
        printf_barry("Num. of Arrays       : 0\n");
        printf_barry("Support size         : -\n");
        printf_barry("Support size range   : -\n");
    }
    

    if (with_pset)
    {
        printf_barry("Arrays in powerset   : %i\n",
            static_cast<int>(std::accumulate(pset_sizes.begin(), pset_sizes.end(), 0u))
        );
    }


    printf_barry("Transform. Fun.      : %s\n", transform_model_fun ? "yes": "no");
    printf_barry("Model terms (%i)    :\n", static_cast<int>(this->nterms()));
    for (auto & cn : this->colnames())
    {
        printf_barry(" - %s\n", cn.c_str());
    }

    if (this->nrules() > 0u)
    {
        printf_barry(
            "Model rules (%i)     :\n",
            static_cast<int>(this->nrules())
        );
    
        for (auto & rn : rules->get_names())
        {
            printf_barry(" - %s\n", rn.c_str());
        }
    }

    if (this->nrules_dyn() > 0u)
    {
        printf_barry(
            "Model rules dyn (%i):\n",
            static_cast<int>(this->nrules_dyn())
        );
    
        for (auto & rn : rules_dyn->get_names())
        {
            printf_barry(" - %s\n", rn.c_str());
        }
    }

    return;

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: size() const noexcept
{
    // INITIALIZED()
    return this->stats_target.size();

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: size_unique() const noexcept
{

    // INITIALIZED()
    return this->stats_support_sizes.size();

} 

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: nterms() const noexcept
{
 
    if (transform_model_fun)
        return transform_model_term_names.size();
    else
        return this->counters->size();

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: nrules() const noexcept
{
 
    return this->rules->size();

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: nrules_dyn() const noexcept
{
 
    return this->rules_dyn->size();

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline size_t Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: support_size() const noexcept
{

    // INITIALIZED()
    return stats_support_sizes_acc.back();
    // size_t tot = 0u;
    // for (auto& a : stats_support)
    //     tot += a.size();

    // return tot;

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< std::string > Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: colnames() const
{
    
    if (transform_model_fun)
        return transform_model_term_names;
    else
        return counters->get_names();

}
    
template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline Array_Type
Model<Array_Type,Data_Counter_Type,Data_Rule_Type, Data_Rule_Dyn_Type>::sample(
    const size_t & i,
    const std::vector<double> & params
) {

    // Are we recording this?
    if (!this->with_pset)
        throw std::logic_error("Sampling is only available when store_pset() is active.");

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    // Getting the index
    size_t a = arrays2support[i];
    
    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    // Updating the current pset
    if (pset_probs.size() == 0u)
        update_pset_probs(params, 1u, static_cast<int>(a));

    // Sampling an array
    size_t j = 0u;
    if (vec_equal_approx(params, params_last[a]))
    // If precomputed, then no need to recalc support
    {

        const double * probs = &pset_probs[pset_locations[a]];
        while (cumprob < r)
            cumprob += *(probs + j++);

        if (j > 0u)
            j--;

    } else { 
       
        update_pset_probs(params, 1u, static_cast<int>(a));

        const double * probs = &pset_probs[pset_locations[a]];
        while (cumprob < r)
            cumprob += *(probs + j++);

        if (j > 0u)
            j--;

        #ifdef BARRY_DEBUG
        if (j > pset_arrays.at(a).size())
            throw std::logic_error(
                std::string(
                    "Something went wrong when sampling from a different set of.") +
                std::string("parameters. Please report this bug: ") +
                std::string(" cumprob: ") + std::to_string(cumprob) +
                std::string(" r: ") + std::to_string(r)
                );
        #endif
        
    }
    
    #ifdef BARRY_DEBUG
    return this->pset_arrays.at(a).at(j);   
    #else
    return this->pset_arrays[a][j];   
    #endif

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Array_Type Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: sample(
    const Array_Type & Array_,
    const std::vector<double> & params
) {

    // Are we recording this?
    if (!this->with_pset)
        throw std::logic_error("Sampling is only available when store_pset() is active.");

    size_t i;

    // If the data hasn't been analyzed earlier, then we need to compute
    // the support
    std::vector< double > key = counters->gen_hash(Array_);
    MapVec_type< double, size_t >::const_iterator locator = keys2support.find(key);
    if (locator == keys2support.end())
    {
        size_t stats_support_size = stats_support.size();

        // Adding to the map
        keys2support[key] = stats_support_sizes.size();
        stats_support_n_arrays.push_back(1u);       // How many elements now
        arrays2support.push_back(stats_support_sizes.size()); // Map of the array id to the support
        
        // Computing support using the counters included in the model
        support_fun.reset_array(Array_);
        
        /** When computing with the powerset, we need to grow the corresponding
            * vectors on the fly */
        if (with_pset)
        {
            
            // Current size of the powerset
            size_t pset_stats_size = pset_stats.size();

            // Making space for storing the support
            pset_arrays.resize(pset_arrays.size() + 1u);
            // pset_stats.resize(pset_stats.size() + 1u);
            // pset_probs.resize(pset_probs.size() + 1u);
            
            try
            {
                
                support_fun.calc(
                    &(pset_arrays[pset_arrays.size() - 1u]),
                    &pset_stats
                );
                
            }
            catch (const std::exception& e)
            {
                
                printf_barry(
                    "A problem ocurred while trying to add the array (and recording the powerset). "
                );
                printf_barry("with error %s\n", e.what());
                throw std::logic_error("");
                
            }

            // Recording the number of elements
            pset_locations.push_back(
                pset_locations.size() == 0u ?
                    0u :
                    pset_locations.back() + pset_sizes.back()
                );

            pset_sizes.push_back(
                (pset_stats.size() - pset_stats_size) / (counter_fun.size())
                );

            // Increasing the space to store probabilities
            pset_probs.resize(pset_probs.size() + pset_sizes.back());

            
        }
        else
        {
            support_fun.calc();
        }
        
        if (transform_model_fun)
        {
            auto tmpsupport = support_fun.get_counts();
            size_t k = counter_fun.size();
            size_t n = tmpsupport.size() / (k + 1);

            std::vector< double > s_new(0u);            
            s_new.reserve(tmpsupport.size());

            for (size_t i = 0u; i < n; ++i)
            {

                // Appending size
                s_new.push_back(tmpsupport[i * (k + 1u)]);

                // Applying transformation and adding to the new set
                auto res = transform_model_fun(&tmpsupport[i * (k + 1u) + 1u], k);
                std::copy(res.begin(), res.end(), std::back_inserter(s_new));

            }

            for (auto & s : s_new)
                stats_support.push_back(s);
            // stats_support.push_back(s_new);

        } else {
            for (auto & s : support_fun.get_counts())
                stats_support.push_back(s);

            // stats_support.push_back(support_fun.get_counts());
        }
        
        // Making room for the previous parameters. This will be used to check if
        // the normalizing constant has been updated or not.
        params_last.push_back(stats_target[0u]);
        normalizing_constants.push_back(0.0);
        first_calc_done.push_back(false);

        // Incrementing the size of the support set
        if (stats_support_sizes.size() == 0u)
        {
            stats_support_sizes_acc.push_back(0u);    
        } else {
            stats_support_sizes_acc.push_back(
                stats_support_sizes.back() + 
                stats_support_sizes_acc.back()
            );
        }


        stats_support_sizes.push_back(
            
            (stats_support.size() - stats_support_size)/
                (counter_fun.size() + 1u)

            );

        
        i = arrays2support.size() - 1u;
    } else
        // Retrieving the corresponding position in the support
        i = locator->second;

    // Getting the index
    size_t a = arrays2support[i];
    
    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    size_t k = params.size();

    // Sampling an array
    size_t j = 0u;
    double * probs = &pset_probs[ pset_locations[a] ];
    if (first_calc_done[a] && (vec_equal_approx(params, params_last[a])))
    // If precomputed, then no need to recalc support
    {

        while (cumprob < r)
            cumprob += *(probs + j++);

        if (j > 0u)
            j--;

    } else { 
       
        // probs.resize(pset_arrays[a].size());
        std::vector< double > temp_stats(params.size());
        const double * stats = &pset_stats[pset_locations[a] * k];

        int i_matches = -1;
        for (size_t array = 0u; array < pset_sizes[a]; ++array)
        {

            // Filling out the parameters
            for (auto p = 0u; p < params.size(); ++p)
                temp_stats[p] = stats[array * k + p];

            *(probs + array) = this->likelihood(params, temp_stats, i, false);
            cumprob += *(probs + array);

            if (i_matches == -1 && cumprob >= r)
                i_matches = array;
        }

        #ifdef BARRY_DEBUG
        if (i_matches < 0)
            throw std::logic_error(
                std::string(
                    "Something went wrong when sampling from a different set of.") +
                std::string("parameters. Please report this bug: ") +
                std::string(" cumprob: ") + std::to_string(cumprob) +
                std::string(" r: ") + std::to_string(r)
                );
        #endif

        j = i_matches;
        first_calc_done[a] = true;
    }
    

    #ifdef BARRY_DEBUG
    return this->pset_arrays.at(a).at(j);   
    #else
    return this->pset_arrays[a][j];   
    #endif

}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline double Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: conditional_prob(
    const Array_Type & Array_,
    const std::vector< double > & params,
    size_t i,
    size_t j
) {

    // Generating a copy of the array so we can update
    Array_Type A(Array_, true);

    // Making sure we add it first
    A.insert_cell(i, j, A.default_val(), true, false);

    // Computing the change stats_target
    std::vector< double > tmp_counts;
    tmp_counts.reserve(counters->size());
    for (size_t ii = 0u; ii < counters->size(); ++ii)
        tmp_counts.push_back(counters->operator[](ii).count(A, i, j));

    // If there is a transformation function, it needs to be
    // applied before dealing with the likelihood.
    if (transform_model_fun)
        tmp_counts = transform_model_fun(&tmp_counts[0u], tmp_counts.size());

    return 1.0/
        (1.0 + std::exp(-vec_inner_prod<double>(
            &params[0u], &tmp_counts[0u], params.size()
            )));

    
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline const std::mt19937 * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_rengine() const {
    return this->rengine;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Counters<Array_Type,Data_Counter_Type> * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_counters() {
    return this->counters;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Rules<Array_Type,Data_Rule_Type> * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_rules() {
    return this->rules;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Rules<Array_Type,Data_Rule_Dyn_Type> * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_rules_dyn() {
    return this->rules_dyn;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_support_fun() {
    return &this->support_fun;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< std::vector< double > > * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_stats_target()
{
    return &stats_target;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< double > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_stats_support()
{
    return &stats_support;
}

// Implementation of get_stats_support_sizes()
template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< size_t > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_stats_support_sizes()
{
    return &stats_support_sizes;
}

// Implementation of get_stats_support_sizes_acc()
template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< size_t > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_stats_support_sizes_acc()
{
    return &stats_support_sizes_acc;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< size_t > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_arrays2support()
{
    return &arrays2support;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector< std::vector< Array_Type > > *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_arrays() {
    return &pset_arrays;
}

template <typename Array_Type, typename Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>
inline std::vector<double> *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_stats() {
    return &pset_stats;
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline std::vector<double> *
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::get_pset_probs() {
    return &pset_probs;
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline std::vector< size_t > * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_pset_sizes()
{
    return &pset_sizes;
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline std::vector< size_t > * Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>:: get_pset_locations()
{
    return &pset_locations;
}

template <
    typename Array_Type,
    typename Data_Counter_Type,
    typename Data_Rule_Type,
    typename Data_Rule_Dyn_Type
    >
inline void
Model<Array_Type,Data_Counter_Type, Data_Rule_Type, Data_Rule_Dyn_Type>::set_transform_model(
    std::function<std::vector<double>(double *,size_t)> fun,
    std::vector< std::string > names
    )
{

    if (transform_model_fun)
        throw std::logic_error("A transformation function for the model has already been established.");
    
    transform_model_fun = fun;
    transform_model_term_names = names;

    size_t k = counters->size(); 

    auto stats_support_old = stats_support;

    // Applying over the support
    for (size_t nsupport = 0u; nsupport < stats_support_sizes.size(); ++nsupport)
    {

        // How many observations in the support
        size_t n = stats_support_sizes[nsupport];

        // Iterating through each observation in the nsupport'th 
        for (size_t i = 0; i < n; ++i)
        {

            // Applying transformation and adding to the new set
            auto res = transform_model_fun(
                &stats_support_old[
                    stats_support_sizes_acc[nsupport] * (k + 1u) +
                    i * (k + 1u) + 1u
                    ],
                k
                );

            if (res.size() != transform_model_term_names.size())
                throw std::length_error(
                    std::string("The transform vector from -transform_model_fun- ") +
                    std::string(" does not match the size of ") + 
                    std::string("-transform_model_term_names-.")
                    );

            // Resizing stats_support if the transform stats do not match the
            // previous size
            if ((nsupport == 0u) && (i == 0u) && (res.size() != k))
                stats_support.resize(
                    (res.size() + 1) * (
                        stats_support_sizes_acc.back() +
                        stats_support_sizes.back()
                        )
                );

            // Weigth
            stats_support[
                stats_support_sizes_acc[nsupport] * (res.size() + 1u) +
                (res.size() + 1u) * i
                ] = stats_support_old[
                    stats_support_sizes_acc[nsupport] * (k + 1u) +
                    i * (k + 1u)
                ];

            // Copying the rest of the elements
            for (size_t j = 0u; j < res.size(); ++j)
                stats_support[
                    stats_support_sizes_acc[nsupport] * (res.size() + 1u) +
                    (res.size() + 1u) * i + j + 1u
                    ] = res[j];

        }

    }

    // Applying over the target statistics
    for (auto & s : stats_target)
        s = transform_model_fun(&s[0u], k);

    // Checking if there is support included
    if (with_pset)
    {

        // Applying it to the support
        for (auto s = 0u; s < pset_arrays.size(); ++s)
        {
            std::vector< double > new_stats;
            size_t pset_stats_loc = pset_locations[s] * k;

            for (auto a = 0u; a < pset_arrays[s].size(); ++a)
            {
                // Computing the transformed version of the data
                auto tmpstats = transform_model_fun(
                    &pset_stats[pset_stats_loc + a * k], k
                    );

                // Storing the new values
                for (auto p = 0u; p < k; ++p)
                    new_stats.push_back(tmpstats[p]);
            }

            // Updating the dataset
            for (size_t stat = 0u; stat < new_stats.size(); ++stat)
                pset_stats[pset_stats_loc + stat] = new_stats[stat];

        }

    }

    // And, resizing the last set of parameters
    for (auto & p : params_last)
        p.resize(transform_model_term_names.size());

    return;

}

#undef MODEL_TEMPLATE
#undef MODEL_TEMPLATE_ARGS
#undef MODEL_TYPE

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/model-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


    
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/rules-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_RULES_BONES_HPP
#define BARRY_RULES_BONES_HPP 1

template <typename Array_Type, typename Data_Type>
bool rule_fun_default(const Array_Type * array, size_t i, size_t j, Data_Type * dat) {
    return false;
}

/**
  * @brief
  * Rule for determining if a cell should be included in a sequence
  * @details
  * Rules can be used together with `Support` and `PowerSet` to determine
  * which cells should be included when enumerating all possible realizations of
  * a binary array.
  * @tparam Array_Type An object of class `BArray`.
  * @tparam Data_Type Any type.
  */
template<typename Array_Type = BArray<>, typename Data_Type = bool>
class Rule {
    
private:
    Rule_fun_type<Array_Type,Data_Type> fun;
    Data_Type dat;
    
    std::string  name = "";
    std::string  desc = "";
    
public:

    /**
     * @name Construct a new Rule object
     * @brief Construct a new Rule object
     * 
     * @param fun_ A function of type `Rule_fun_type`.
     * @param dat_ Data pointer to be passed to `fun_`
     * @param delete_dat_ When `true`, the `Rule` destructor will delete the
     * pointer, if defined.
     */
    ///@{
    Rule() : fun(rule_fun_default<Array_Type,Data_Type>) {};
    Rule(
        Rule_fun_type<Array_Type,Data_Type> fun_,
        Data_Type dat_,
        std::string name_        = "",   
        std::string desc_        = ""
        ) : fun(fun_), dat(dat_), name(name_), desc(desc_) {};
    ///@}

    ~Rule() {};

    Data_Type & D(); ///< Read/Write access to the data.
    
    bool operator()(const Array_Type & a, size_t i, size_t j);

    std::string & get_name();
    std::string & get_description();

    std::string get_name() const;
    std::string get_description() const;
    
};

/**
  * @brief Vector of objects of class Rule
  * 
  * @tparam Array_Type An object of class `BArray`
  * @tparam Data_Type Any type.
  */
template<typename Array_Type, typename Data_Type>
class Rules {

private:
    std::vector< Rule<Array_Type,Data_Type> > data;
    
public:
    Rules() {};

    Rules(const Rules<Array_Type,Data_Type> & rules_);
    Rules<Array_Type,Data_Type> operator=(const Rules<Array_Type,Data_Type> & rules_);

    ~Rules() {};

    size_t size() const noexcept {
        return data.size();
    };
    
    /**
      * @name Rule adding
      * 
      * @param rule 
      */
    ///@{
    void add_rule(Rule<Array_Type, Data_Type> rule);
    void add_rule(
        Rule_fun_type<Array_Type,Data_Type> rule_,
        Data_Type data_,
        std::string name_ = "",
        std::string description_ = ""
    );
    ///@}

    /**
     * @brief Check whether a given cell is free or locked
     * 
     * @param a A `BArray` object
     * @param i row position
     * @param j col position
     * @return true If the cell is locked
     * @return false If the cell is free
     */

    bool operator()(const Array_Type & a, size_t i, size_t j);
    
    /**
     * @brief Computes the sequence of free and locked cells in an BArray
     * 
     * @param a An object of class `BArray`.
     * @param free Pointer to a vector of pairs (i, j) listing the free cells.
     * @param locked (optional) Pointer to a vector of pairs (i, j) listing the
     * locked cells.
     * @return Nothing.
     */
    void get_seq(
        const Array_Type & a,
        std::vector< size_t > * free,
        std::vector< size_t > * locked = nullptr
    );

    std::vector< std::string > get_names() const;
    std::vector< std::string > get_descriptions() const;

    // Iterator
    typename std::vector< Rule<Array_Type,Data_Type> >::iterator begin() {
        return data.begin();
    };
    typename std::vector< Rule<Array_Type,Data_Type> >::iterator end() {
        return data.end();
    };
    
};


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/rules-bones.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/rules-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_RULES_MEAT_HPP
#define BARRY_RULES_MEAT_HPP 1

template <typename Array_Type, typename Data_Type>
inline Rules<Array_Type,Data_Type>::Rules(
    const Rules<Array_Type,Data_Type> & rules_
) {

    // Copy all rules, if a rule is tagged as 
    // to be deleted, then copy the value
    for (auto i = 0u; i != rules_.size(); ++i)
        this->add_rule(rules_.data[i]);

    return;

}

template <typename Array_Type, typename Data_Type>
Rules<Array_Type,Data_Type> Rules<Array_Type,Data_Type>::operator=(
    const Rules<Array_Type,Data_Type> & rules_
) {

    if (this != &rules_) {

        // Copy all rules, if a rule is tagged as 
        // to be deleted, then copy the value
        for (auto i = 0u; i != rules_.size(); ++i)
            this->add_rule(rules_.data[i]);

    }

    return *this;

}

template<typename Array_Type, typename Data_Type>
inline Data_Type & Rule<Array_Type,Data_Type>::D()
{
    return dat;
}

template<typename Array_Type, typename Data_Type>
inline bool Rule<Array_Type,Data_Type>::operator()(const Array_Type & a, size_t i, size_t j) {
    return fun(a, i, j, dat);
}

template<typename Array_Type, typename Data_Type>
inline std::string & Rule<Array_Type,Data_Type>::get_name()
{
    return name;
}

template<typename Array_Type, typename Data_Type>
inline std::string & Rule<Array_Type,Data_Type>::get_description()
{
    return desc;
}

template<typename Array_Type, typename Data_Type>
inline std::string Rule<Array_Type,Data_Type>::get_name() const
{
    return name;
}

template<typename Array_Type, typename Data_Type>
inline std::string Rule<Array_Type,Data_Type>::get_description() const
{
    return desc;
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::add_rule(
        Rule<Array_Type, Data_Type> rule
) {
    
    data.push_back(rule);
    
    return;
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::add_rule(
        Rule_fun_type<Array_Type,Data_Type> rule_,
        Data_Type                           data_,
        std::string name_,
        std::string description_
) {
       
    data.push_back(Rule<Array_Type,Data_Type>(
        rule_,
        data_,
        name_,
        description_
    ));
    
    return;
    
}

template <typename Array_Type, typename Data_Type>
inline bool Rules<Array_Type,Data_Type>::operator()(
    const Array_Type & a, size_t i, size_t j
) {
    
    if (data.size()==0u)
        return true;
    
    for (auto & f: data)
        if (!f.operator()(a, i, j))
            return false;
    
    return true;
    
}

template <typename Array_Type, typename Data_Type>
inline void Rules<Array_Type,Data_Type>::get_seq(
    const Array_Type & a,
    std::vector< size_t > * free,
    std::vector< size_t > * locked
) {

    
    size_t N = a.nrow();
    size_t K = a.ncol();
    
    // Reserving some space
    (void) free->empty();
    (void) free->reserve(2u * N * K);
    
    for (size_t i = 0u; i < N; ++i)
    {

        for (size_t j = 0u; j < K; ++j)
        {

            // Locked cells are skipped
            if (!this->operator()(a, i, j))
            {

                if (locked != nullptr)
                {

                    locked->push_back(i);
                    locked->push_back(j);

                }

                continue;

            }

            free->push_back(i);
            free->push_back(j);
                
        }

    }
    
    free->shrink_to_fit();

    return;

}

template<typename Array_Type, typename Data_Type>
inline std::vector<std::string> Rules<Array_Type, Data_Type>::get_names() const
{

    std::vector< std::string > out;
    out.reserve(this->size());
    for (size_t i = 0u; i < this->size(); ++i)
        out.push_back(this->data.at(i).get_name());

    return out;

}

template<typename Array_Type, typename Data_Type>
inline std::vector<std::string> Rules<Array_Type, Data_Type>::get_descriptions() const
{
    
    std::vector< std::string > out;
    out.reserve(this->size());
    for (size_t i = 0u; i < this->size(); ++i)
        out.push_back(this->data.at(i).get_description());

    return out;

}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/rules-meat.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


    
    namespace counters {
        namespace network {
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/counters/network.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRAY_NETWORK_H
#define BARRAY_NETWORK_H 1

/**
 * @ingroup counting 
 * @details Details on the available counters for `NetworkData` can be found in
 * the \ref counters-network section.
 * 
 */
///@{

/**
 * @brief Data class for Networks.
 * 
 * This holds information about whether the graph is directed or not, and,
 * if defined, vectors of node (vertex) attributes (`vertex_attr`).
 * 
 */
class NetworkData {
public:
    
    bool directed = true;
    std::vector< std::vector< double > > vertex_attr;
    
    NetworkData() : vertex_attr(0u) {};
    
    /**
     * @brief Constructor using a single attribute
     * @param vertex_attr_ Double vector of length equal to the number of vertices
     * in the data.
     * @param directed_ When `true` the graph as treated as directed.
     */
    NetworkData(
        std::vector< double >  vertex_attr_,
        bool directed_ = true
    ) : directed(directed_), vertex_attr(1u, vertex_attr_) {};
    
    /**
     * @brief Constructor using multiple attributes
     * @param vertex_attr_ Vector of double vectors. The size equals to the number
     * of attributes to be created. Each individual vector should be of length
     * equal to the number of vertices.
     * @param directed_ When `true` the graph as treated as directed.
     */
    NetworkData(
        std::vector< std::vector< double > > vertex_attr_,
        bool directed_ = true
    ) : directed(directed_), vertex_attr(vertex_attr_) {};
    
    
    ~NetworkData() {};
};

/**
  * @brief Data class used to store arbitrary size_t or double vectors */
class NetCounterData {
public:
    
    std::vector< size_t > indices;
    std::vector< double > numbers;
    
    NetCounterData() : indices(0u), numbers(0u) {};
    NetCounterData(
        const std::vector< size_t > & indices_,
        const std::vector< double > & numbers_
    ): indices(indices_), numbers(numbers_) {};
    
    ~NetCounterData() {};
    
    // const size_t get_size_t
    
};

#define NET_C_DATA_IDX(i) (data.indices[i])
#define NET_C_DATA_NUM(i) (data.numbers[i])


/**
 * @name Convenient typedefs for network objects.
 */
///@{
typedef BArray<double, NetworkData> Network;
typedef BArrayDense<int, NetworkData> NetworkDense;

#define BARRY_ZERO_NETWORK 0.0
#define BARRY_ZERO_NETWORK_DENSE 0

template <typename Tnet = Network>
using NetCounter =  Counter<Tnet, NetCounterData >;

template <typename Tnet = Network>
using NetCounters =  Counters<Tnet, NetCounterData>;

template <typename Tnet = Network>
using NetSupport =  Support<Tnet, NetCounterData >;

template <typename Tnet = Network>
using NetStatsCounter =  StatsCounter<Tnet, NetCounterData>;

template <typename Tnet>
using NetModel =  Model<Tnet, NetCounterData>;

template <typename Tnet = Network>
using NetRule =  Rule<Tnet, bool>;

template <typename Tnet = Network>
using NetRules =  Rules<Tnet, bool>;
///@}

/**@name Macros for defining counters
  */
///@{
/**Function for definition of a network counter function*/
#define NETWORK_COUNTER(a) \
template<typename Tnet = Network>\
inline double (a) (const Tnet & Array, size_t i, size_t j, NetCounterData & data)

/**Lambda function for definition of a network counter function*/
#define NETWORK_COUNTER_LAMBDA(a) \
Counter_fun_type<Tnet, NetCounterData> a = \
    [](const Tnet & Array, size_t i, size_t j, NetCounterData & data)

#define NETWORKDENSE_COUNTER_LAMBDA(a) \
Counter_fun_type<NetworkDense, NetCounterData> a = \
    [](const NetworkDense & Array, size_t i, size_t j, NetCounterData & data)
///@}


/**@name Macros for defining rules
  */
///@{
/**Function for definition of a network counter function*/
#define NETWORK_RULE(a) \
template<typename Tnet = Network>\
inline bool (a) (const Tnet & Array, size_t i, size_t j, bool & data)

/**Lambda function for definition of a network counter function*/
#define NETWORK_RULE_LAMBDA(a) \
Rule_fun_type<Tnet, bool> a = \
[](const Tnet & Array, size_t i, size_t j, bool & data)
///@}

/**
  * @weakgroup  counters-network Network counters
  * @brief Counters for network models
  * @param counters A pointer to a `NetCounters` object (`Counters`<`Network`, `NetCounterData`>).
  */
///@{
// -----------------------------------------------------------------------------
/**@brief Number of edges */
template<typename Tnet = Network>
inline void counter_edges(NetCounters<Tnet> * counters)
{
    
    NETWORK_COUNTER_LAMBDA(count_edges)
    {
        return 1.0;
    };
    
    counters->add_counter(
        count_edges, nullptr, nullptr,
        NetCounterData(), 
        "Edge counts", 
        "Number of edges"
        );
    
    return;

}


// -----------------------------------------------------------------------------
/**@brief Number of isolated vertices */
template<typename Tnet = Network>
inline void counter_isolates(NetCounters<Tnet> * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        if (i == j)
            return 0.0;
        
        double res = 0.0;
        
        // i is sending its first tie
        if (Array.row(i).size() == 1u && Array.col(i).size() == 0u)
            res -= 1.0;
        
        // j is receiving its first tie, meaning that he
        // has no other tie but i's?
        if (Array.row(j).size() == 0u && Array.col(j).size() == 1u)
            res -= 1.0;
        
        return res;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        return static_cast<double>(Array.nrow());
    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        NetCounterData(),
        "Isolates",
        "Number of isolate vertices"
        );

    return;
}

template<>
inline void counter_isolates(NetCounters<NetworkDense> * counters)
{
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {

        if (i == j)
            return 0.0;
        
        double res = 0.0;
        
        // Checking the in and out degree
        if (Array.rowsum(i) == 1u && Array.colsum(i) == 0u)
            res -= 1.0;

        // Now looking at j
        if (Array.rowsum(j) == 0u && Array.colsum(j) == 1u)
            res -= 1.0;
        
        return res;

    };
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_init)
    {
        return static_cast<double>(Array.nrow());
    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        NetCounterData(),
        "Isolates", "Number of isolate vertices"
        );

    return;

}

// -----------------------------------------------------------------------------
/**@brief Number of mutual ties */
template<typename Tnet = Network>
inline void counter_mutual(NetCounters<Tnet> * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        // Is there any tie at ji? If not, then we have a new mutual!
        // but this only makes sence if the jth row and ith column exists
        // if ((Array.nrow() > j) && (Array.ncol() > i)) 
        if (i == j)
            return 0.0;
        
        // printf_barry("Checking if it is empty or not at (%i, %i)... ", i, j);
        if (!Array.is_empty(j, i, false))
        {
            // printf_barry("Yes, mutual.\n");
            return 1.0;
        }
        // printf_barry("No, no mutual.\n");
        
        return 0.0;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {

        if (Array.nrow() != Array.ncol())
            throw std::logic_error("The -mutual- counter only works on square arrays.");
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!Array.D_ptr()->directed)
            throw std::logic_error(
                "The -mutual- counter only works on directed (non-symmetric) arrays."
                );
        
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        NetCounterData(),
        "Reciprocity",
        "Number of mutual ties"
        );

    return;

}


// 2-istars --------------------------------------------------------------------
template<typename Tnet = Network>
inline void counter_istar2(NetCounters<Tnet> * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        // Need to check the receiving, if he/she is getting a new set of stars
        // when looking at triads
        
        if (Array.col(j).size() == 1u)
            return 0.0;
        
        return static_cast<double>(Array.col(j).size() - 1.0);

    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr,
        NetCounterData(),
        "Istar 2",
        "Indegree 2-star"
        );
    
    return ;
}

template<>
inline void counter_istar2(NetCounters<NetworkDense> * counters)
{
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {
        // Need to check the receiving, if he/she is getting a new set of stars
        // when looking at triads
        // int indeg = 1;
        // for (size_t k = 0u; k < Array.nrow(); ++k)
        // {
        //     if (i == k)
        //         continue;

        //     if (Array(k,j) != BARRY_ZERO_NETWORK_DENSE)
        //         indeg++;
        // }

        // if (indeg == 1)
        //     return 0.0;
        
        // return static_cast<double>(indeg - 1);
        return static_cast<double>(Array.colsum(j) - 1);

    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr,
        NetCounterData(),
        "Istar 2",
        "Indegree 2-star"
        );
    
    return ;
}


// 2-ostars --------------------------------------------------------------------
template<typename Tnet = Network>
inline void counter_ostar2(NetCounters<Tnet> * counters)
{
   
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        // Need to check the receiving, if he/she is getting a new set of stars
        // when looking at triads
        
        if (Array.row(i).size() == 1u)
            return 0.0;
        
        return static_cast<double>( Array.row(i).size() - 1.0);

    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr,
        NetCounterData(),
        "Ostar 2",
        "Outdegree 2-star"
        );

    return ;
    
}

template<>
inline void counter_ostar2(NetCounters<NetworkDense> * counters)
{
   
    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {

        // Need to check the receiving, if he/she is getting a new set of stars
        // when looking at triads
        // int nties = 0;
        // for (size_t k = 0u; k < Array.ncol(); ++k)
        // {
        //     if (Array(i, k) != BARRY_ZERO_NETWORK_DENSE)
        //         ++nties;
        // }

        // if (nties == 1u)
        //     return 0.0;
        
        // return static_cast<double>(nties - 1.0);
        return static_cast<double>(Array.rowsum(i) - 1);

    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr,
        NetCounterData(),
        "Ostar 2",
        "Outdegree 2-star"
        );

    return ;
    
}


// ttriads ---------------------------------------------------------------------
template<typename Tnet = Network>
inline void counter_ttriads(NetCounters<Tnet> * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        // Self ties do not count
        if (i == j)
            return 0.0;
        
        double ans = 0.0;
        
        // Case 1: i-j, i-k, j-k
        if (Array.row(j).size() < Array.row(i).size())
        {
            
            for (auto j_row = Array.row(j).begin(); j_row != Array.row(j).end(); ++j_row) 
                if ((j != j_row->first) && (i != j_row->first) && !Array.is_empty(i, j_row->first, false))
                    ans += 1.0;
                
        } else {
            
            for (auto i_row = Array.row(i).begin(); i_row != Array.row(i).end(); ++i_row) 
                if ((i != i_row->first) && (i_row->first != j) && !Array.is_empty(j, i_row->first, false))
                    ans += 1.0;
                
        }
        
        // Case 2: i-j, i-k, k-j  
        if (Array.row(i).size() > Array.col(j).size())
        {
            
            for (auto j_col = Array.col(j).begin(); j_col != Array.col(j).end(); ++j_col)
                if ((j != j_col->first) && (i != j_col->first) && !Array.is_empty(i, j_col->first, false))
                    ans += 1.0;
                
        } else {
            
            for (auto i_row = Array.row(i).begin(); i_row != Array.row(i).end(); ++i_row) 
                if ((i != i_row->first) && (j != i_row->first) && !Array.is_empty(i_row->first, j, false))
                    ans += 1.0;
                
        }
        
        // Case 3: i->j, k->j, k->i
        if (Array.col(i).size() > Array.col(j).size())
        {
            
            for (auto j_col = Array.col(j).begin(); j_col != Array.col(j).end(); ++j_col)
                if ((j != j_col->first) && (i != j_col->first) && !Array.is_empty(j_col->first, i, false))
                    ans += 1.0;
                
        } else {
            
            for (auto i_col = Array.col(i).begin(); i_col != Array.col(i).end(); ++i_col) 
                if ((i != i_col->first) && (j != i_col->first) && !Array.is_empty(i_col->first, j, false))
                    ans += 1.0;
                
        }
        
        // The regular counter double counts
        return ans;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!(Array.D_ptr()->directed))
            throw std::invalid_argument("The ttriads counter is only valid for directed networks. This is undirected.");

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        NetCounterData(),
        "Balance",
        "Number of directed triangles"
    );
    
    return;

}

template<>
inline void counter_ttriads(NetCounters<NetworkDense> * counters)
{
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {

        const auto & dat = Array.get_data();
        size_t N = Array.nrow();

        // Self ties do not count
        if (i == j)
            return 0.0;

        // This is the first i sends, so nothing will change
        if (Array.rowsum(i) == BARRY_ZERO_NETWORK_DENSE)
            return 0.0;

        
        double ans = 0.0;
        for (size_t k = 0u; k < N; ++k)
        {

            // In all cases k receives, so if not, then continue
            if ((Array.colsum(k) == BARRY_ZERO_NETWORK_DENSE) && (Array.rowsum(k) == BARRY_ZERO_NETWORK_DENSE))
                continue;

            if ((j != k) & (i != k))
            {

                if (dat[k * N + i] != BARRY_ZERO_NETWORK_DENSE)
                {
                    // Case 1: i-j, i-k, j-k
                    if (dat[k * N + j])
                        ans += 1.0;

                    // Case 2: i-j, i-k, k-j 
                    if (dat[j * N + k] != BARRY_ZERO_NETWORK_DENSE)
                        ans += 1.0;
                }
                
                // Case 3: i-j, k-i, k-j
                if ((dat[i * N + k] != BARRY_ZERO_NETWORK_DENSE) && (dat[j * N + k] != BARRY_ZERO_NETWORK_DENSE))
                    ans += 1.0;

            }
        }
        
        // The regular counter double counts
        return ans;

    };
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!(Array.D_ptr()->directed))
            throw std::invalid_argument("The ttriads counter is only valid for directed networks. This is undirected.");

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        NetCounterData(),
        "Balance",
        "Number of directed triangles"
    );
    
    return;

}


// Cycle triads --------------------------------------------------------------
template<typename Tnet = Network>
inline void counter_ctriads(NetCounters<Tnet> * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        if (i == j)
            return 0.0;
        
        double ans = 0.0;
        if (Array.col(i).size() < Array.row(j).size())
        {
            
            for (auto i_col = Array.col(i).begin(); i_col != Array.col(i).end(); ++i_col) 
                if ((i != i_col->first) && (j != i_col->first) && !Array.is_empty(j, i_col->first, false))
                    ans += 1.0;
                
        } else {
            
            for (auto j_row = Array.row(j).begin(); j_row != Array.row(j).end(); ++j_row) 
                if ((j != j_row->first) && (i != j_row->first) && !Array.is_empty(j_row->first, i, false))
                    ans += 1.0;
                
        }
        
        return ans;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {

        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!(Array.D_ptr()->directed))
            throw std::invalid_argument(
                "The ctriads counter is only valid for directed networks. This is undirected."
                );

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        NetCounterData(),
        "Cyclical triads"
    );

    return;
    
}

template<>
inline void counter_ctriads(NetCounters<NetworkDense> * counters)
{
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {

        if (i == j)
            return 0.0;
        
        // i->j->k->i
        double ans = 0.0;
        #if defined(__OPENMP) || defined(_OPENMP) 
        #pragma omp simd reduction(+:ans)
        #endif
        for (size_t k = 0u; k < Array.nrow(); ++k)
        {

            // If isolated, then next
            if (Array.colsum(k) == BARRY_ZERO_NETWORK_DENSE)
                continue;

            if (Array.rowsum(k) == BARRY_ZERO_NETWORK_DENSE)
                continue;

            if (i != k && j != k)
            {

                if ((Array(j, k) != BARRY_ZERO_NETWORK_DENSE) && (Array(k, i) != BARRY_ZERO_NETWORK_DENSE))
                    ans += 1.0;

            }
    }
        
        return ans;

    };
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_init)
    {

        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!(Array.D_ptr()->directed))
            throw std::invalid_argument(
                "The ctriads counter is only valid for directed networks. This is undirected."
                );

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        NetCounterData(),
        "Cyclical triads"
    );

    return;
    
}
    
// Density --------------------------------------------------------------
template<typename Tnet = Network>
inline void counter_density(NetCounters<Tnet> * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return
            1.0/(Array.nrow() * (Array.ncol() - 1.0)) / (
                (Array.D_ptr()->directed)? 1.0 : 2.0
            );
        
    };
    
    // Preparing the counter data and returning. We make sure that the memory is 
    // released so we set delete_data = true.
    counters->add_counter(
        tmp_count, nullptr, nullptr,
        NetCounterData(),
        "Density",
        "Proportion of present ties"
    );

    return ;
    
}

// idegree1.5  -------------------------------------------------------------
template<typename Tnet = Network>
inline void counter_idegree15(NetCounters<Tnet> * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        // In case of the first, we need to add
        if (Array.col(j).size() == 1u)
            return 1.0;
        
        return 
            pow(static_cast<double> (Array.col(j).size()), 1.5) -
            pow(static_cast<double> (Array.col(j).size() - 1), 1.5)
            ;
        
    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr,
        NetCounterData(),
        "Indegree^(1.5)"
    );

    return;
    
}

template<>
inline void counter_idegree15(NetCounters<NetworkDense> * counters)
{
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {
        
        // In case of the first, we need to add
        int ideg = 0;
        for (size_t k = 0u; k < Array.nrow(); ++k)
        {
            if (k == j)
                continue;

            if (Array(k, j) != BARRY_ZERO_NETWORK_DENSE)
                ideg++;

        }

        if (ideg == 0)
            return 0.0;
        
        if (ideg == 1)
            return 1.0;

        double res = std::pow(static_cast<double> (ideg), 1.5) -
            std::pow(static_cast<double> (ideg - 1.0), 1.5);

        if (std::isnan(res))
            throw std::domain_error("Resulting indeg is undefined.");
        
        return 
            std::pow(static_cast<double> (ideg), 1.5) -
            std::pow(static_cast<double> (ideg - 1.0), 1.5)
            ;
        
    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr,
        NetCounterData(),
        "Indegree^(1.5)"
    );

    return;
    
}

// odegree1.5  -------------------------------------------------------------
template<typename Tnet = Network>
inline void counter_odegree15(NetCounters<Tnet> * counters)
{
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        // In case of the first, we need to add
        if (Array.row(i).size() == 1u)
            return 1.0;
        
        return 
            pow(static_cast<double>(Array.row(i).size()), 1.5) -
            pow(static_cast<double>(Array.row(i).size() - 1), 1.5)
            ;
        
    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr,
        NetCounterData(),
        "Outdegree^(1.5)"
    );

    return;
    
}

template<>
inline void counter_odegree15(NetCounters<NetworkDense> * counters)
{
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {
        
        // In case of the first, we need to add
        int odeg = 0;
        for (size_t k = 0u; k < Array.ncol(); ++k)
        {

            if (k == i)
                continue;

            if (Array(i, k) != BARRY_ZERO_NETWORK_DENSE)
                odeg++;

        }

        if (odeg == 0)
            return 0.0;

        if (odeg == 1)
            return 1.0;
        
        return 
            pow(static_cast<double>(odeg), 1.5) -
            pow(static_cast<double>(odeg - 1), 1.5)
            ;
        
    };
    
    counters->add_counter(
        tmp_count, nullptr, nullptr,
        NetCounterData(),
        "Outdegree^(1.5)"
    );

    return;
    
}


// -----------------------------------------------------------------------------
/**@brief Sum of absolute attribute difference between ego and alter */
template<typename Tnet = Network>
inline void counter_absdiff(
    NetCounters<Tnet> * counters,
    size_t attr_id,
    double alpha = 1.0
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return std::pow(std::fabs(
                Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][i] - 
                    Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][j]
        ), NET_C_DATA_NUM(0u));
        
    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (Array.D_ptr()->vertex_attr.size() == 0u)
            throw std::range_error("No attributes in the Array.");
        
        if ((NET_C_DATA_IDX(0u) != 0u) && (Array.D_ptr()->vertex_attr.size() <= (NET_C_DATA_IDX(0u) - 1u)))
            throw std::range_error("Attribute index out of range.");
        
        return 0.0;
        
    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        NetCounterData({attr_id}, {alpha}),
        "Absdiff"
    );
    
    return;
    
}
    
// -----------------------------------------------------------------------------
/**@brief Sum of attribute difference between ego and alter to pow(alpha)*/
template<typename Tnet = Network>
inline void counter_diff(
    NetCounters<Tnet> * counters,
    size_t attr_id,
    double alpha     = 1.0,
    double tail_head = true
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return std::pow(NET_C_DATA_NUM(1u) * (
                Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][i] - 
                    Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][j]
        ), NET_C_DATA_NUM(0u));
        
    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (Array.D_ptr()->vertex_attr.size() == 0u)
            throw std::range_error("No attributes in the Array.");
        
        if ((NET_C_DATA_IDX(0u) != 0u) && (Array.D_ptr()->vertex_attr.size() <= (NET_C_DATA_IDX(0u) - 1u)))
            throw std::range_error("Attribute index out of range.");
        
        return 0.0;
        
    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        NetCounterData({attr_id}, {alpha, tail_head ? 1.0: -1.0}),
        "Absdiff^(" + std::to_string(alpha) + ")"
    );
    
    return;
    
}

// Nodeicov, nodeocov, and Nodematch -------------------------------------------
NETWORK_COUNTER(init_single_attr)
{
    
    if (Array.D_ptr() == nullptr)
        throw std::logic_error("The array data has not been initialized");
    
    if (Array.D_ptr()->vertex_attr.size() == 0u)
        throw std::range_error("No attributes in the Array.");
    
    if ((NET_C_DATA_IDX(0u) != 0u) && (Array.D_ptr()->vertex_attr.size() <= (NET_C_DATA_IDX(0u) - 1u)))
        throw std::range_error("Attribute index out of range.");
    
    return 0.0;
    
}

// -----------------------------------------------------------------------------
//*@brief Attribute sum over receiver nodes */
template<typename Tnet = Network>
inline void counter_nodeicov(
    NetCounters<Tnet> * counters,
    size_t attr_id
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][j];
        
    };
    
    counters->add_counter(
        tmp_count, init_single_attr<Tnet>, nullptr,
        NetCounterData({attr_id}, {}),
        "nodeicov", "Sum of ego attribute"
    );
      
    return;

}

// -----------------------------------------------------------------------------
//*@brief Attribute sum over sender nodes */
template<typename Tnet = Network>
inline void counter_nodeocov(
    NetCounters<Tnet> * counters,
    size_t attr_id
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][i];
        
    };
    
    counters->add_counter(
        tmp_count, init_single_attr<Tnet>, nullptr,
        NetCounterData({attr_id}, {}),
        "nodeocov", "Sum of alter attribute"
    );
    
    return;

}

// -----------------------------------------------------------------------------
//*@brief Attribute sum over receiver and sender nodes */
template<typename Tnet = Network>
inline void counter_nodecov(
    NetCounters<Tnet> * counters,
    size_t attr_id
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][i] +
            Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][j];
        
    };
    
    counters->add_counter(
        tmp_count, init_single_attr<Tnet>, nullptr,
        NetCounterData({attr_id}, {}),
        "nodecov", "Sum of nodes covariates"
    );
    
    return;
}

// -----------------------------------------------------------------------------
//* @brief Number of homophililic ties */
template<typename Tnet = Network>
inline void counter_nodematch(
    NetCounters<Tnet> * counters,
    size_t attr_id
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        return 
        (
                Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][i] == 
                    Array.D_ptr()->vertex_attr[NET_C_DATA_IDX(0u)][j]
        )? 1.0 : 0.0;
        
    };
    
    // Preparing the counter data and returning. We make sure that the memory is 
    // released so we set delete_data = true.
    counters->add_counter(
        tmp_count, init_single_attr<Tnet>, nullptr,
        NetCounterData({attr_id}, {}),
        "Homophily",
        "Number of homophilic ties"
    );
    
    return ;
    
}

// -----------------------------------------------------------------------------
/** @brief Counts number of vertices with a given in-degree */
template<typename Tnet = Network>
inline void counter_idegree(
    NetCounters<Tnet> * counters,
    std::vector< size_t > d
) {

    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        size_t d = Array.col(j).size();
        if (d == NET_C_DATA_IDX(0u))
            return 1.0;
        else if (d == (NET_C_DATA_IDX(0u) + 1))
            return -1.0;
        
        return 0.0;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!Array.D_ptr()->directed)
            throw std::logic_error("-odegree- counter is only valid for directed graphs");
        
        if (NET_C_DATA_IDX(0u) == 0u)
            return static_cast<double>(Array.nrow());
        
        return 0.0;

    };
    
    for (auto iter = d.begin(); iter != d.end(); ++iter)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            NetCounterData({*iter}, {}),
            "Nodes indeg " + std::to_string(*iter),
            "Number of nodes with indigree " + std::to_string(*iter)
        );
    
    return;  

}

template<>
inline void counter_idegree(
    NetCounters<NetworkDense> * counters,
    std::vector< size_t > d
) {

    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {
        
        size_t indeg = 0u;
        for (size_t k = 0u; k < Array.nrow(); ++k)
            if (Array(k, j) != BARRY_ZERO_NETWORK_DENSE)
                indeg++;

        if (indeg == NET_C_DATA_IDX(0u))
            return 1.0;
        else if (indeg == (NET_C_DATA_IDX(0u) + 1))
            return -1.0;
        
        return 0.0;

    };
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!Array.D_ptr()->directed)
            throw std::logic_error("-odegree- counter is only valid for directed graphs");
        
        if (NET_C_DATA_IDX(0u) == 0u)
            return static_cast<double>(Array.nrow());
        
        return 0.0;

    };
    
    for (auto iter = d.begin(); iter != d.end(); ++iter)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            NetCounterData({*iter}, {}),
            "Nodes indeg " + std::to_string(*iter),
            "Number of nodes with indigree " + std::to_string(*iter)
        );
    
    return;  

}

// -----------------------------------------------------------------------------
/**@brief Counts number of vertices with a given out-degree */
template<typename Tnet = Network>
inline void counter_odegree(
    NetCounters<Tnet> * counters,
    std::vector<size_t> d
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        size_t d = Array.row(i).size();
        if (d == NET_C_DATA_IDX(0u))
            return 1.0;
        else if (d == (NET_C_DATA_IDX(0u) + 1))
            return -1.0;
        
        return 0.0;

    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!Array.D_ptr()->directed)
            throw std::logic_error("-odegree- counter is only valid for directed graphs");
        
        if (NET_C_DATA_IDX(0u) == 0u)
            return static_cast<double>(Array.nrow());
        
        return 0.0;

    };
        
        
    for (auto iter = d.begin(); iter != d.end(); ++iter) 
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            NetCounterData({*iter}, {}),
            "Nodes w/ outdeg " + std::to_string(*iter),
            "Number of nodes with outdegree " + std::to_string(*iter)
        );
    
    return;  
    
}
    
template<>
inline void counter_odegree(
    NetCounters<NetworkDense> * counters,
    std::vector<size_t> d
) {
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {
        
        size_t d = 0;
        for (size_t k = 0u; k < Array.ncol(); ++k)
            if (Array(i, k) != BARRY_ZERO_NETWORK_DENSE)
                d++;
        
        if (d == NET_C_DATA_IDX(0u))
            return 1.0;
        else if (d == (NET_C_DATA_IDX(0u) + 1))
            return -1.0;
        
        return 0.0;

    };
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_init)
    {
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (!Array.D_ptr()->directed)
            throw std::logic_error("-odegree- counter is only valid for directed graphs");
        
        if (NET_C_DATA_IDX(0u) == 0u)
            return static_cast<double>(Array.nrow());
        
        return 0.0;

    };
        
        
    for (auto iter = d.begin(); iter != d.end(); ++iter) 
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            NetCounterData({*iter}, {}),
            "Nodes w/ outdeg " + std::to_string(*iter),
            "Number of nodes with outdegree " + std::to_string(*iter)
        );
    
    return;  
    
}


// -----------------------------------------------------------------------------
/** @brief Counts number of vertices with a given out-degree */
template<typename Tnet = Network>
inline void counter_degree(
    NetCounters<Tnet> * counters,
    std::vector<size_t> d
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {
        
        size_t d = Array.row(i).size();
        if (d == NET_C_DATA_IDX(0u))
            return 1.0;
        else if (d == (NET_C_DATA_IDX(0u) + 1))
            return -1.0;
        
        return 0.0;
    };
    
    NETWORK_COUNTER_LAMBDA(tmp_init) {
        
        if (Array.D_ptr() == nullptr)
            throw std::logic_error("The array data has not been initialized");
        
        if (Array.D_ptr()->directed)
            throw std::logic_error("-degree- counter is only valid for undirected graphs");
        
        if (NET_C_DATA_IDX(0u) == 0u)
            return static_cast<double>(Array.nrow());
        
        return 0.0;
    };
    
    
    for (auto iter = d.begin(); iter != d.end(); ++iter)
    {
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            NetCounterData({*iter}, {})
        );
    }
    
    return;  
}

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry//counters/network-css.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_CSS_COUNTERS
#define BARRY_CSS_COUNTERS

// n: Net size, 
// s: Start of the i-th network
// e: end of the i-th network
// ego_id: Ego of the cell (i, j)
#define CSS_SIZE() \
    size_t n      = data.indices[0u]; \
    size_t s      = data.indices[1u]; \
    size_t e      = data.indices[2u]; \
    size_t ctype  = data.indices[3u]; \
    size_t ego_id = data.indices[4u]; \
    if (ctype > 2) \
        throw std::range_error("Counter type should be 0, 1, or 2.");

// Check whether ego_id is involved in the current cell
// ctype: Type of counter
// 0: All cells
// 1: Only if perceiver
// 2: Only if not perceiver
#define CSS_MATCH_TYPE() \
    if (ctype != 0u) { /* all counts */ \
        if (ctype == 1u) { /* Only if perceiver */ \
            if ((i_ != ego_id) && (j_ != ego_id)) return 0.0; \
        } else if (ctype == 2u) { /* Only if not perceiver */ \
            if ((i_ == ego_id) || (j_ == ego_id)) return 0.0; \
        } \
    };

// Variables in case that the current cell corresponds to the True
#define CSS_CASE_TRUTH() if ((i < n) && (j < n)) 

// i_: i-th index of the cell
// j_: j-th index of the cell
// tji: True value of the cell (i, j)
// pij: Perceived value of the cell (i, j)
// pji: Perceived value of the cell (j, i)
#define CSS_TRUE_CELLS() \
    size_t i_ = i; \
    size_t j_ = j; \
    CSS_MATCH_TYPE() \
    double tji = static_cast<double>(Array(j, i, false)); \
    double pij = static_cast<double>(Array(i + s, j + s, false)); \
    double pji = static_cast<double>(Array(j + s, i + s, false));

// Variables in case that the current cell corresponds to the Perceived
#define CSS_CASE_PERCEIVED() else if (((i >= s) && (i < e)) & ((j >= s) && (j < e)))

// i_: i-th index of the cell
// j_: j-th index of the cell
// tji: True value of the cell (i, j)
// pji: Perceived value of the cell (i, j)
// tij: True value of the cell (j, i)
#define CSS_PERCEIVED_CELLS() \
    size_t i_ = i - s; \
    size_t j_ = j - s; \
    CSS_MATCH_TYPE() \
    double tji = static_cast<double>(Array(j - s, i - s, false)); \
    double pji = static_cast<double>(Array(j, i, false)); \
    double tij = static_cast<double>(Array(i - s, j - s, false));



// Nothing for else (for now)
#define CSS_CASE_ELSE()

// Checks whether the start and end of the node (perceived network) falls within
// the boundaries of the graph.
#define CSS_CHECK_SIZE_INIT() \
    /* The indices fall within the network */ \
    if ((data.indices.at(0) > Array.ncol()) \
    | (data.indices.at(2) > Array.ncol())) \
        throw std::range_error("The network does not match the prescribed size."); 

#define CSS_CHECK_SIZE() for (size_t i = 0u; i < end_.size(); ++i) {\
    if (i == 0u) continue; \
    else if (end_[i] < end_[i-1u]) \
        throw std::logic_error("Endpoints should be specified in order.");}

#define CSS_APPEND(name) std::string name_ = (name);\
    for (size_t i = 0u; i < end_.size(); ++i) { \
    std::string tmpname = name_ + " (" + std::to_string(i) + ")" + \
    ((counter_type == 1u) ? " (only perceiver)" : ((counter_type == 2u)? " (only alters)": ""));\
    counters->add_counter(tmp_count, tmp_init, nullptr, \
            NetCounterData({netsize, i == 0u ? netsize : end_[i-1], end_[i], counter_type, i}, {}),\
            tmpname);}

#define CSS_NET_COUNTER_LAMBDA_INIT() NETWORK_COUNTER_LAMBDA(tmp_init) {\
        CSS_CHECK_SIZE_INIT() \
        return 0.0; \
    };


/**
 * @brief Counts errors of commission 
 * @param netsize Size of the reference (true) network 
 * @param end_ Vector indicating one past the ending index of each network. (see details)
 * @param counter_type Size_t indicating the type of counter to use. Possible
 *  values are: 0: Count all, 1: Only count if perceiver is involved, and 
 *  2: Only count if perceiver is not involved.
 * @details 
 * The `end_` parameter should be of length `N of networks` - 1. It is
 * assumed that the first network ends at `netsize`.
 */
template<typename Tnet = Network>
inline void counter_css_partially_false_recip_commi(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {
        

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            // Checking change stat of the true net
            CSS_TRUE_CELLS()
            return pij * pji * (1.0 - 2.0 * tji) - (1.0 - tji)*(
                pij * (1.0 - pji) + (1.0 - pij) * pji
            );

        } CSS_CASE_PERCEIVED() {

            // Checking change stat of the percieved net
            CSS_PERCEIVED_CELLS()
            return pji * (tij * (1.0 - tji) + (1.0 - tij) * tji) +
                (1.0 - tij) * (1.0 - tji) * (1 - 2.0 * pji)
            ;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };

    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Partially false recip (comission)")

    return;
    
}

/** @brief Counts errors of omission */
template<typename Tnet = Network>
inline void counter_css_partially_false_recip_omiss(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {
        
        // Getting the network size
        CSS_SIZE()
        
        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * ((1.0 - pij) * pji + pij * (1.0 - pji)) +
                (1.0 - 2.0 * tji) * (1.0 - pij) * (1.0 - pji)
            ;

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return tji * tij * (1.0 - 2.0 * pji) - 
                (1.0 - pji) * ((1.0 - tij) * tji + tij * (1.0 - tji))
            ;

        } CSS_CASE_ELSE()
            return 0.0;

    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Partially false recip (omission)")
        
    return;
    
}

/** @brief Counts completely false reciprocity (comission) */
template<typename Tnet = Network>
inline void counter_css_completely_false_recip_comiss(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return -(1.0 - tji) * pij * pji;

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return (1.0 - tij) * (1.0 - tji) * pji;

        } CSS_CASE_ELSE()
            return 0.0;

    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Completely false recip (comission)")
        
    return;
    
}

/** @brief Counts completely false reciprocity (omission) */
template<typename Tnet = Network>
inline void counter_css_completely_false_recip_omiss(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * (1.0 - pij) * (1.0 - pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return - tij * tji * (1.0 - pji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Completely false recip (omission)")
        
    return;
    
}

/** @brief Counts mixed reciprocity errors */
template<typename Tnet = Network>
inline void counter_css_mixed_recip(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return (1.0 - tji) * (1.0 - pij) * pji - tji * pij * (1.0 - pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return (1.0 - tij) * tji * (1.0 - pji) - tij * (1.0 - tij) * pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("Mixed reciprocity errors")
        
    return;
    
}

/////////////////////////// CENSUS

template<typename Tnet = Network>
inline void counter_css_census01(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count)
    {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return -(1.0 - pij) * (1.0 - pji) * (1.0 - tji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return -(1.0 - tij) * (1.0 - tji) * (1.0 - pji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    // CSS_NET_COUNTER_LAMBDA_INIT()
    NETWORK_COUNTER_LAMBDA(tmp_init)
    {

        CSS_CHECK_SIZE_INIT()
        double n_dbl = static_cast<double>(data.indices[0u]);

        // Discount in case of the type of counter
        size_t ctype = data.indices[3u];

        if (ctype == 1u) /* Only perceiver */
        {

            return (n_dbl - 1.0); // * (Array.D().directed ? 2.0 : 1.0);

        } else if (ctype == 2u) /* All but the perceiver */
        {
            // We remove the perceiver from the eq.
            n_dbl -= 1.0;
        }

        // At the beginning is all zero
        return n_dbl * (n_dbl - 1.0) / 2.0; // / (Array.D().directed ? 1.0 : 2.0);

    };
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(01) Accurate null")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census02(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return -(1.0 - tji) * ((1.0 - pij) * pji + pij * (1.0 - pji));

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return (1.0 - tij) * (1.0 - tji) * (1 - 2.0 * pji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(02) Partial false positive (null)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census03(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return -(1.0 - tji) * pij * pji;

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return (1.0 - tij) * (1.0 - tji) *pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(03) Complete false positive (null)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census04(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return (1.0 - pij) * (1.0 - pji) * (1.0 - 2.0 * tji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return -(1.0 - pji) * ((1.0 - tij) * tji + tij * (1.0 - tji));

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(04) Partial false negative (assym)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census05(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return pij * (1.0 - tji) * (1.0 - pji) - (1.0 - pij) * tji * pji;

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return tij * (1.0 - tji) * (1.0 - pji) - (1.0 - tij) * tji * pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(05) Accurate assym")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census06(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return (1.0 - pij) * (1.0 - tji) * pji - pij * tji * (1.0 - pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return (1.0 - tij) * tji * (1.0 - pji) - tij * (1.0 - tji) * pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(06) Mixed assym")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census07(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return pij * pji * (1.0 - 2.0 * tji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return pji * (tij * (1.0 - tji) + (1.0 - tij) * tji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(07) Partial false positive (assym)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census08(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * (1.0 - pij) * (1.0 - pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return - tij * tji * (1.0 - pji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(08) Complete false negative (full)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census09(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * (pij * (1.0 - pji) + (1.0 - pij) * pji);

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return tij * tji * (1.0 - 2.0 * pji);

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(09) Partial false negative (full)")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census10(
    NetCounters<Tnet> * counters,
    size_t netsize,
    const std::vector< size_t > & end_,
    size_t counter_type = 0u
) {

    NETWORK_COUNTER_LAMBDA(tmp_count) {

        // Getting the network size
        CSS_SIZE()

        // True network
        CSS_CASE_TRUTH()
        {

            CSS_TRUE_CELLS()
            return tji * pij * pji;

        } CSS_CASE_PERCEIVED() {

            CSS_PERCEIVED_CELLS()
            return tij * tji * pji;

        } CSS_CASE_ELSE()
            return 0.0;
        
    };
    
    CSS_NET_COUNTER_LAMBDA_INIT()
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(10) Accurate full")
        
    return;

}

#undef CSS_APPEND
#undef CSS_CASE_TRUTH
#undef CSS_TRUE_CELLS
#undef CSS_CASE_PERCEIVED
#undef CSS_PERCEIVED_CELLS
#undef CSS_CASE_ELSE
#undef CSS_CHECK_SIZE_INIT
#undef CSS_CHECK_SIZE
#undef CSS_NET_COUNTER_LAMBDA_INIT
#undef CSS_MATCH_TYPE
#undef CSS_SIZE
#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry//counters/network-css.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



///@}


/**
 * @name Rules for network models
 * @param rules A pointer to a `NetRules` object (`Rules`<`Network`, `bool`>).
 */
///@{
// -----------------------------------------------------------------------------
/**@brief Number of edges */
template<typename Tnet = Network>
inline void rules_zerodiag(NetRules<Tnet> * rules) {
    
    NETWORK_RULE_LAMBDA(no_self_tie) {
        return i != j;
    };
    
    rules->add_rule(
        no_self_tie, false,
        "No self-ties",
        "No self-ties"
        );
    
    return;
}

///@}

///@}

#undef NET_C_DATA_IDX
#undef NET_C_DATA_NUM

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/counters/network.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


        }
    }
    
}

namespace netcounters = barry::counters::network;

#define COUNTER_FUNCTION(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
    inline double (a) (const Array_Type & Array, size_t i, size_t j, Data_Type & data)\

#define COUNTER_LAMBDA(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
    Counter_fun_type<Array_Type, Data_Type> a = \
    [](const Array_Type & Array, size_t i, size_t j, Data_Type & data)

#define RULE_FUNCTION(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
    inline bool (a) (const Array_Type & Array, size_t i, size_t j, Data_Type & data)\

#define RULE_LAMBDA(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
    Rule_fun_type<Array_Type, Data_Type> a = \
    [](const Array_Type & Array, size_t i, size_t j, Data_Type & data)

#endif
