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

#ifdef BARRY_USE_OMP
#include <omp.h>
#endif

#ifndef BARRY_HPP
#define BARRY_HPP 

#define BARRY_VERSION_MAYOR 0
#define BARRY_VERSION_MINOR 1
#define BARRY_VERSION BARRY_VERSION_MAYOR ## . ## BARRY_VERSION_MINOR

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
    #define BARRY_MAX_NUM_ELEMENTS static_cast< size_t >(UINT_MAX/2u)
#endif

#ifdef BARRY_USE_OMP
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
        for(const auto & iter : (a)) 
            printf_barry("%.4f ", static_cast<double>(iter));
        printf_barry("]\n");
        return;
    }

    // Specialization for the printer
    template<>
    inline void BARRY_DEBUG_VEC_PRINT(const std::vector< int > & a) {
        printf_barry("%s  [", BARRY_DEBUG_HEADER);
        for(const auto & iter : (a)) 
            printf_barry("%i ", iter);
        printf_barry("]\n");
        return;
    }

    template<>
    inline void BARRY_DEBUG_VEC_PRINT(const std::vector< std::string > & a) {
        printf_barry("%s \n", BARRY_DEBUG_HEADER);
        for(const auto & iter : (a)) 
            printf_barry("%s %s\n", BARRY_DEBUG_HEADER, iter.c_str());
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
        printf_barry("|");

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
 * @param unit, uint Focal cell
 * @param Data_Type Data associated with the function, for example, id of the attribute
 *  in the Array.
 * @return `Counter_fun_type` a double (the change statistic)
 * @return `Rule_fun_type` a bool. True if the cell is blocked.
 */
///@{
template <typename Array_Type, typename Data_Type>
using Counter_fun_type = std::function<double(const Array_Type &, uint, uint, Data_Type &)>;

template <typename Array_Type, typename Data_Type>
using Rule_fun_type = std::function<bool(const Array_Type &, uint, uint, Data_Type &)>;
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
    #pragma GCC ivdep
    #endif
    for (unsigned int i = 0u; i < n; ++i)
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
    #pragma GCC ivdep
    #endif
    for (unsigned int i = 0u; i < n; ++i)
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
    // uint ncols;
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

    for (unsigned int i = 0u; i < n; ++i)
    {
        
        std::vector< double > tmp(k, 0.0);

        for (unsigned j = 1u; j < (k + 1u); ++j)
            tmp[j - 1u] = data[i * (k + 1) + j];
        
        ans.push_back(
            std::make_pair<std::vector<double>,unsigned int>(
                std::move(tmp),
                static_cast<unsigned int>(data[i * (k + 1u)])
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

    unsigned int grand_total = 0u;

    printf_barry("%7s | %s\n", "Counts", "Stats");

    for (unsigned int i = 0u; i < n; ++i)
    {

        printf_barry("%7i | ", static_cast<int>(data[i * (k + 1u)]));

        for (unsigned int j = 1u; j < (k + 1u); ++j)
            printf_barry(" %.2f", data[i * (k + 1) + j]);
        printf_barry("\n");

        grand_total += static_cast<unsigned int>(data[i * (k + 1u)]);

    }

    printf_barry("Grand total: %i\n", grand_total);

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
        for (unsigned int i = 1u; i < x.size(); ++i)
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

template <> inline void Cell<unsigned int>::add(unsigned int x) {
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
template<> inline Cell< uint >::Cell() : value(1u), visited(false), active(true) {}
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
 * `std::vector< std::unordered_map<unsigned int,Cell_Type> >`.
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
    uint N;
    uint M;
    uint NCells = 0u;
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
     * @param target An unsigned int vector ranging from 0 to M_
     * @param target When `true` tries to add repeated observations.
     */
    ///@{
    
    /** @brief Zero-size array */
    BArray() : N(0u), M(0u), NCells(0u), el_ij(0u), el_ji(0u) {};
    
    /** @brief Empty array */
    BArray (uint N_, uint M_) : N(N_), M(M_), NCells(0u), el_ij(N_), el_ji(M_) {};
    
    /** @brief Edgelist with data */
    BArray (
        uint N_, uint M_,
        const std::vector< uint > & source,
        const std::vector< uint > & target,
        const std::vector< Cell_Type > & value,
        bool add = true
    );
    
    /** @brief Edgelist with no data (simpler) */
    BArray (
        uint N_, uint M_,
        const std::vector< uint > & source,
        const std::vector< uint > & target,
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
    void out_of_range(uint i, uint j) const;
    Cell_Type get_cell(uint i, uint j, bool check_bounds = true) const; 
    std::vector< Cell_Type >      get_col_vec(uint i, bool check_bounds = true) const;
    std::vector< Cell_Type >      get_row_vec(uint i, bool check_bounds = true) const;
    void                          get_col_vec(std::vector< Cell_Type > * x, uint i, bool check_bounds = true) const;
    void                          get_row_vec(std::vector< Cell_Type > * x, uint i, bool check_bounds = true) const;
    const Row_type< Cell_Type > & row(uint i, bool check_bounds = true) const;
    const Col_type< Cell_Type > & col(uint i, bool check_bounds = true) const;

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
    bool is_empty(uint i, uint j, bool check_bounds = true) const;
    uint nrow() const noexcept;
    uint ncol() const noexcept;
    uint nnozero() const noexcept;
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
    BArray<Cell_Type,Data_Type> & operator+=(const std::pair<uint, uint> & coords);
    BArray<Cell_Type,Data_Type> & operator-=(const std::pair<uint, uint> & coords);
    BArrayCell<Cell_Type,Data_Type> operator()(uint i, uint j, bool check_bounds = true);
    const Cell_Type operator()(uint i, uint j, bool check_bounds = true) const;
    
    void rm_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
    
    void insert_cell(uint i, uint j, const Cell< Cell_Type > & v, bool check_bounds, bool check_exists);
    void insert_cell(uint i, uint j, Cell< Cell_Type > && v, bool check_bounds, bool check_exists);
    void insert_cell(uint i, uint j, Cell_Type v, bool check_bounds, bool check_exists);
    
    void swap_cells(
        uint i0, uint j0, uint i1, uint j1, bool check_bounds = true,
        int check_exists = CHECK::BOTH,
        int * report     = nullptr
        );
    
    void toggle_cell(uint i, uint j, bool check_bounds = true, int check_exists = EXISTS::UKNOWN);
    void toggle_lock(uint i, uint j, bool check_bounds = true);
    ///@}
    
    /**@name Column/row wise interchange*/
    ///@{
    void swap_rows(uint i0, uint i1, bool check_bounds = true);
    void swap_cols(uint j0, uint j1, bool check_bounds = true);
    
    void zero_row(uint i, bool check_bounds = true);
    void zero_col(uint j, bool check_bounds = true);
    ///@}
    
    void transpose();
    void clear(bool hard = true);
    void resize(uint N_, uint M_);
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
    // operator BArray<uint,bool>() const;
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
    uint i;
    uint j;
  
public:
  
    BArrayCell(BArray<Cell_Type,Data_Type> * Array_, uint i_, uint j_, bool check_bounds = true) : 
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
    uint i;
    uint j;
    
public:
  
    BArrayCell_const(const BArray<Cell_Type,Data_Type> * Array_, uint i_, uint j_, bool check_bounds = true) : 
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

#define BARRAY_TYPE() BArray<Cell_Type, Data_Type>

#define BARRAY_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BARRAY_TEMPLATE(a,b) \
    template BARRAY_TEMPLATE_ARGS() inline a BARRAY_TYPE()::b

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]


template<typename Cell_Type, typename Data_Type>
Cell<Cell_Type> BArray<Cell_Type,Data_Type>::Cell_default = Cell<Cell_Type>(static_cast<Cell_Type>(1.0)); 


// Edgelist with data
BARRAY_TEMPLATE(,BArray) (
    uint N_, uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
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
    for (uint i = 0u; i < source.size(); ++i) {
      
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
BARRAY_TEMPLATE(,BArray) (
    uint N_, uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
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
    for (uint i = 0u; i < source.size(); ++i) {
      
        // Checking range
        if ((source[i] >= N_) | (target[i] >= M_))
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
            std::pair<uint, Cell< Cell_Type> >(
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

BARRAY_TEMPLATE(,BArray) (
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
    for (uint i = 0u; i < N; ++i)
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

BARRAY_TEMPLATE(BARRAY_TYPE() &, operator=) (
    const BArray<Cell_Type,Data_Type> & Array_
) {
  
    // Clearing
    if (this != &Array_)
    {
      
        this->clear(true);
        this->resize(Array_.N, Array_.M);
        
        // Entries
        for (uint i = 0u; i < N; ++i)
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

BARRAY_TEMPLATE(,BArray) (
    BARRAY_TYPE() && x
  ) noexcept :
  N(0u), M(0u), NCells(0u),
  data(nullptr),
  delete_data(x.delete_data)
  {

    this->clear(true);
    this->resize(x.N, x.M);
    
    // Entries
    for (uint i = 0u; i < N; ++i) {
      
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

BARRAY_TEMPLATE(BARRAY_TYPE() &, operator=) (
    BARRAY_TYPE() && x
) noexcept {
  
    // Clearing
    if (this != &x) {
      
        this->clear(true);
        this->resize(x.N, x.M);
        
        // Entries
        for (uint i = 0u; i < N; ++i) {
          
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

BARRAY_TEMPLATE(bool, operator==) (
    const BARRAY_TYPE() & Array_
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

BARRAY_TEMPLATE(,~BArray) () {
    
    if (delete_data && (data != nullptr))
        delete data;
    
    return;
}

BARRAY_TEMPLATE(void, set_data) (
    Data_Type * data_, bool delete_data_
) {  

    if ((data != nullptr) && delete_data)
        delete data;
    
    data        = data_;
    delete_data = delete_data_;
    
    return;
    
}

BARRAY_TEMPLATE(Data_Type *, D_ptr) ()
{
    return this->data;
}

template<typename Cell_Type, typename Data_Type>
inline const Data_Type * BArray<Cell_Type,Data_Type>::D_ptr() const
{
    return this->data;
}

BARRAY_TEMPLATE(Data_Type &, D) ()
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

BARRAY_TEMPLATE(void, out_of_range) (
    uint i,
    uint j
) const {

    if (i >= N)
        throw std::range_error("The row is out of range.");
    else if (j >= M)
        throw std::range_error("The column is out of range.");
    return;

}
    
BARRAY_TEMPLATE(Cell_Type, get_cell) (
    uint i,
    uint j,
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

BARRAY_TEMPLATE(std::vector< Cell_Type >, get_row_vec) (
    uint i,
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

BARRAY_TEMPLATE(void, get_row_vec) (
    std::vector< Cell_Type > * x,
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    for (const auto & iter : row(i, false)) 
        x->at(iter.first) = iter.second.value; // this->get_cell(i, iter->first, false);
    
}

BARRAY_TEMPLATE(std::vector< Cell_Type >, get_col_vec) (
    uint i,
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

BARRAY_TEMPLATE(void, get_col_vec) (
    std::vector<Cell_Type> * x,
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    for (const auto & iter : col(i, false)) 
        x->at(iter.first) = iter.second->value;//this->get_cell(iter->first, i, false);
    
}

BARRAY_TEMPLATE(const Row_type< Cell_Type > &, row) (
    uint i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, 0u);

    return this->el_ij[i];

}

BARRAY_TEMPLATE(const Col_type< Cell_Type > &, col) (
    uint i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(0u, i);

    return this->el_ji[i];
    
}

BARRAY_TEMPLATE(Entries< Cell_Type >, get_entries) () const {
    
    Entries<Cell_Type> res(NCells);
    
    for (uint i = 0u; i < N; ++i) {
        
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

BARRAY_TEMPLATE(bool, is_empty) (
    uint i,
    uint j,
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


BARRAY_TEMPLATE(unsigned int, nrow) () const noexcept {
    return N;
}


BARRAY_TEMPLATE(unsigned int, ncol) () const noexcept {
    return M;
}


BARRAY_TEMPLATE(unsigned int, nnozero) () const noexcept {
    return NCells;
}

BARRAY_TEMPLATE(Cell< Cell_Type >, default_val) () const {
    return this->Cell_default;
}

BARRAY_TEMPLATE(BARRAY_TYPE() &, operator+=) (
    const std::pair<uint,uint> & coords
) {
    
    this->insert_cell(
            coords.first,
            coords.second,
            this->Cell_default,
            true, true
    );
    
    return *this;
    
}

BARRAY_TEMPLATE(BARRAY_TYPE() &, operator-=) (
    const std::pair<uint,uint> & coords
) {
    
    this->rm_cell(
            coords.first,
            coords.second,
            true, true
    );
    
    return *this;
    
}

template BARRAY_TEMPLATE_ARGS()
inline BArrayCell<Cell_Type,Data_Type> BARRAY_TYPE()::operator()(  
    uint i,
    uint j,
    bool check_bounds
) {
    
    return BArrayCell<Cell_Type,Data_Type>(this, i, j, check_bounds);
    
}

template BARRAY_TEMPLATE_ARGS()
inline const Cell_Type BARRAY_TYPE()::operator() (  
    uint i,
    uint j,
    bool check_bounds
) const {
    
    return get_cell(i, j, check_bounds);
    
}

BARRAY_TEMPLATE(void, rm_cell) (
    uint i,
    uint j,
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

BARRAY_TEMPLATE(void, insert_cell) (
        uint i,
        uint j,
        const Cell< Cell_Type> & v,
        bool check_bounds,
        bool check_exists
    ) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    if (check_exists) {
        
        // Checking if nothing here, then we move along
        if (ROW(i).size() == 0u) {
            
            ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v));
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            return;
            
        }
        
        // In this case, the row exists, but we are checking that the value is empty  
        if (ROW(i).find(j) == ROW(i).end()) {
            
            ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v)); 
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            
        } else {
            throw std::logic_error("The cell already exists.");
        }
        
        
    } else {
        
        ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v));
        COL(j).emplace(i, &ROW(i)[j]);
        NCells++;
        
    }
    
    return;
    
}

BARRAY_TEMPLATE(void, insert_cell) (
        uint i,
        uint j,
        Cell< Cell_Type> && v,
        bool check_bounds,
        bool check_exists
    ) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    if (check_exists) {
        
        // Checking if nothing here, then we move along
        if (ROW(i).size() == 0u) {
            
            ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v));
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            return;
            
        }
        
        // In this case, the row exists, but we are checking that the value is empty  
        if (ROW(i).find(j) == ROW(i).end()) {
            
            ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v)); 
            COL(j).emplace(i, &ROW(i)[j]);
            NCells++;
            
        } else {
            throw std::logic_error("The cell already exists.");
        }
        
        
    } else {
        
        ROW(i).insert(std::pair< uint, Cell<Cell_Type>>(j, v));
        COL(j).emplace(i, &ROW(i)[j]);
        NCells++;
        
    }
    
    return;
    
}

BARRAY_TEMPLATE(void, insert_cell) (
    uint i,
    uint j,
    Cell_Type v,
    bool check_bounds,
    bool check_exists
) {
        
    return insert_cell(i, j, Cell<Cell_Type>(v, visited), check_bounds, check_exists);

}

BARRAY_TEMPLATE(void, swap_cells) (
    uint i0, uint j0,
    uint i1, uint j1,
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

BARRAY_TEMPLATE(void, toggle_cell) (
    uint i,
    uint j,
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

BARRAY_TEMPLATE(void, swap_rows) (
    uint i0,
    uint i1,
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
BARRAY_TEMPLATE(void, swap_cols) (
    uint j0,
    uint j1,
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

BARRAY_TEMPLATE(void, zero_row) (
    uint i,
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

BARRAY_TEMPLATE(void, zero_col) (
    uint j,
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

BARRAY_TEMPLATE(void, transpose) () {
  
    // Start by flipping the switch 
    visited = !visited;
    
    // Do we need to resize (increase) either?
    if      (N > M) el_ji.resize(N);
    else if (N < M) el_ij.resize(M);
    
    // uint N0 = N, M0 = M;
    int status;
    for (uint i = 0u; i < N; ++i)
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

BARRAY_TEMPLATE(void, clear) (
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
        
        for (unsigned int i = 0u; i < N; ++i)
            zero_row(i, false);
        
    }
    
    return;
    
}

BARRAY_TEMPLATE(void, resize) (
    uint N_,
    uint M_
) {
  
    // Removing rows
    if (N_ < N)
        for (uint i = N_; i < N; ++i)
            zero_row(i, false);
    
    // Removing cols
    if (M_ < M)
        for (uint j = M_; j < M; ++j)
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

BARRAY_TEMPLATE(void, reserve) () {
#ifdef BARRAY_USE_UNORDERED_MAP
    for (uint i = 0u; i < N; i++)
        ROW(i).reserve(M);
    
    for (uint i = 0u; i < M; i++)
        COL(i).reserve(N);
#endif
    return;
  
}

BARRAY_TEMPLATE(void, print) (
    const char * fmt,
    ...
) const {
  
    std::va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);

    for (uint i = 0u; i < N; ++i)
    {

        #ifdef BARRY_DEBUG_LEVEL
            #if BARRY_DEBUG_LEVEL > 1
                printf_barry("%s [%3i,]", BARRY_DEBUG_HEADER, i);
            #endif
        #else
        printf_barry("[%3i,] ", i);
        #endif
        for (uint j = 0u; j < M; ++j) {
            if (this->is_empty(i, j, false))
                printf_barry("    . ");
            else 
                printf_barry(" %.2f ", static_cast<double>(this->get_cell(i, j, false)));
            
        }

        printf_barry("\n");

    }
    
    
    return;

}

#undef ROW
#undef COL

#undef BARRAY_TYPE
#undef BARRAY_TEMPLATE_ARGS
#undef BARRAY_TEMPLATE

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
    
    for (uint i = 0u; i < nrow(); ++i)
        for (uint j = 0u; j < ncol(); ++j)
            this->operator()(i, j) += rhs.get_cell(i, j);

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator+=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
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
    
    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) -= rhs.get_cell(i, j);
        }
    }

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator-=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) -= rhs;
        }
    }

    return *this;
}

BARRAY_TEMPLATE(BARRAY_TYPE()&, operator*=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {

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

    for (uint i = 0u; i < nrow(); ++i) {

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
    uint N;
    uint M;
    // uint NCells = 0u;
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
     * @param target An unsigned int vector ranging from 0 to M_
     * @param target When `true` tries to add repeated observations.
     * @param value Cell_Type defaul fill-in value (zero, by default.)
     */
    ///@{
    
    /** @brief Zero-size array */
    BArrayDense() : N(0u), M(0u), el(0u), el_rowsums(0u), el_colsums(0u) {};
    
    /** @brief Empty array */
    BArrayDense (uint N_, uint M_, Cell_Type value = static_cast<Cell_Type>(0)) :
        N(N_), M(M_), el(N_ * M_, value),
        el_rowsums(N_, static_cast<Cell_Type>(value * M_)), el_colsums(M_, static_cast<Cell_Type>(value * N_)) {};
    
    /** @brief Edgelist with data */
    BArrayDense (
        uint N_,
        uint M_,
        const std::vector< uint > & source,
        const std::vector< uint > & target,
        const std::vector< Cell_Type > & value,
        bool add = true
    );
    
    /** @brief Edgelist with no data (simpler) */
    BArrayDense (
        uint N_, uint M_,
        const std::vector< uint > & source,
        const std::vector< uint > & target,
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
    void out_of_range(uint i, uint j) const;
    Cell_Type get_cell(uint i, uint j, bool check_bounds = true) const; 
    std::vector< Cell_Type >      get_col_vec(uint i, bool check_bounds = true) const;
    std::vector< Cell_Type >      get_row_vec(uint i, bool check_bounds = true) const;
    void                          get_col_vec(std::vector< Cell_Type > * x, uint i, bool check_bounds = true) const;
    void                          get_row_vec(std::vector< Cell_Type > * x, uint i, bool check_bounds = true) const;
    
    BArrayDenseRow<Cell_Type,Data_Type> & row(uint i, bool check_bounds = true);
    const BArrayDenseRow_const<Cell_Type,Data_Type> row(uint i, bool check_bounds = true) const;

    BArrayDenseCol<Cell_Type,Data_Type> & col(uint j, bool check_bounds = true);
    const BArrayDenseCol_const<Cell_Type,Data_Type> col(uint j, bool check_bounds = true) const;

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
    bool is_empty(uint i, uint j, bool check_bounds = true) const;
    uint nrow() const noexcept;
    uint ncol() const noexcept;
    uint nnozero() const noexcept;
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
    BArrayDense<Cell_Type,Data_Type> & operator+=(const std::pair<uint, uint> & coords);
    BArrayDense<Cell_Type,Data_Type> & operator-=(const std::pair<uint, uint> & coords);
    BArrayDenseCell<Cell_Type,Data_Type> operator()(uint i, uint j, bool check_bounds = true);
    const Cell_Type operator()(uint i, uint j, bool check_bounds = true) const;
    
    void rm_cell(uint i, uint j, bool check_bounds = true, bool check_exists = true);
    
    void insert_cell(uint i, uint j, const Cell< Cell_Type > & v, bool check_bounds, bool check_exists);
    // void insert_cell(uint i, uint j, Cell< Cell_Type > && v, bool check_bounds, bool check_exists);
    void insert_cell(uint i, uint j, Cell_Type v, bool check_bounds, bool check_exists);
    
    void swap_cells(
        uint i0, uint j0, uint i1, uint j1, bool check_bounds = true,
        int check_exists = CHECK::BOTH,
        int * report     = nullptr
        );
    
    void toggle_cell(uint i, uint j, bool check_bounds = true, int check_exists = EXISTS::UKNOWN);
    void toggle_lock(uint i, uint j, bool check_bounds = true);
    ///@}
    
    /**@name Column/row wise interchange*/
    ///@{
    void swap_rows(uint i0, uint i1, bool check_bounds = true);
    void swap_cols(uint j0, uint j1, bool check_bounds = true);
    
    void zero_row(uint i, bool check_bounds = true);
    void zero_col(uint j, bool check_bounds = true);
    ///@}
    
    void transpose();
    void clear(bool hard = true);
    void resize(uint N_, uint M_);
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
    // operator BArrayDense<uint,bool>() const;
    // operator BArrayDense<bool,bool>() const;
    // ///@}
    
    bool is_dense() const noexcept {return dense;};

    const std::vector< Cell_Type > & get_data() const;
    const Cell_Type rowsum(unsigned int i) const;
    const Cell_Type colsum(unsigned int i) const;
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
    uint i;
    uint j;
  
public:
  
    BArrayDenseCell(
        BArrayDense<Cell_Type,Data_Type> * Array_,
        uint i_,
        uint j_,
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
    unsigned int index;
    bool row_filled = false; // after row is filled

    void fill_if_needed()
    {
        if (!row_filled)
        {

            for (unsigned int j = 0u; j < array->M; ++j)
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
        unsigned int i
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

    std::pair<unsigned int,Cell<Cell_Type>> & operator()(unsigned int i)
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
    unsigned int index;

public:
    BArrayDenseRow_const(
        const BArrayDense< Cell_Type,Data_Type > & array_,
        unsigned int i
    ) : array(&array_), index(i)
    {

        for (unsigned int j = 0u; j < array->M; ++j)
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

    const std::pair<unsigned int,Cell<Cell_Type>> operator()(unsigned int i) const
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
    unsigned int index;
    bool col_filled = false;

    void fill_if_needed()
    {
        if (!col_filled)
        {

            for (unsigned int i = 0u; i < array->N; ++i)
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
        unsigned int j
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

    std::pair<unsigned int,Cell_Type*> & operator()(unsigned int i)
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
    unsigned int index;
    Col_type<Cell_Type> col;

public:
    BArrayDenseCol_const(
        const BArrayDense< Cell_Type,Data_Type > & array_,
        unsigned int j
    ) : array(&array_), index(j)
    {

        for (unsigned int i = 0u; i < array->N; ++i)
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

    const std::pair<unsigned int,Cell_Type*> operator()(unsigned int i) const
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


#define BDENSE_TYPE() BArrayDense<Cell_Type, Data_Type>

#define BDENSE_TEMPLATE_ARGS() <typename Cell_Type, typename Data_Type>

#define BDENSE_TEMPLATE(a,b) \
    template BDENSE_TEMPLATE_ARGS() inline a BDENSE_TYPE()::b

#define ROW(a) this->el_ij[a]
#define COL(a) this->el_ji[a]
#define POS(a,b) (b)*N + (a)
#define POS_N(a,b,c) (b)*(c) + (a)

template<typename Cell_Type, typename Data_Type>
Cell_Type BArrayDense<Cell_Type,Data_Type>::Cell_default = static_cast< Cell_Type >(1.0); 

#define ZERO_CELL static_cast<Cell_Type>(0.0)

// Edgelist with data
BDENSE_TEMPLATE(,BArrayDense)(
    uint N_,
    uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
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

    el.resize(N * M, ZERO_CELL);
    el_rowsums.resize(N, ZERO_CELL);
    el_colsums.resize(M, ZERO_CELL);
    
    // Writing the data
    for (uint i = 0u; i < source.size(); ++i)
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
BDENSE_TEMPLATE(, BArrayDense)(
    uint N_, uint M_,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    bool add
) {
  
    std::vector< Cell_Type > value(source.size(), static_cast<Cell_Type>(1.0));

    if (source.size() != target.size())
        throw std::length_error("-source- and -target- don't match on length.");
    if (source.size() != value.size())
        throw std::length_error("-sorce- and -value- don't match on length.");
    
    // Initializing
    N = N_;
    M = M_;

    el.resize(N * M, ZERO_CELL);
    el_rowsums.resize(N, ZERO_CELL);
    el_colsums.resize(M, ZERO_CELL);
    
    // Writing the data
    for (uint i = 0u; i < source.size(); ++i)
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

BDENSE_TEMPLATE(, BArrayDense)(
    const BDENSE_TYPE() & Array_,
    bool copy_data
) : N(Array_.N), M(Array_.M){
  
    // Dimensions
    el.resize(0u);
    el_rowsums.resize(0u);
    el_colsums.resize(0u);
    
    std::copy(Array_.el.begin(), Array_.el.end(), std::back_inserter(el));
    std::copy(Array_.el_rowsums.begin(), Array_.el_rowsums.end(), std::back_inserter(el_rowsums));
    std::copy(Array_.el_colsums.begin(), Array_.el_colsums.end(), std::back_inserter(el_colsums));

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

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator=) (
    const BDENSE_TYPE() & Array_
) {
  
    // Clearing
    if (this != &Array_)
    {
      
        el.resize(0u);
        el_rowsums.resize(0u);
        el_colsums.resize(0u);
        
        // Entries
        std::copy(Array_.el.begin(), Array_.el.end(), std::back_inserter(el));
        std::copy(Array_.el_rowsums.begin(), Array_.el_rowsums.end(), std::back_inserter(el_rowsums));
        std::copy(Array_.el_colsums.begin(), Array_.el_colsums.end(), std::back_inserter(el_colsums));


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

BDENSE_TEMPLATE(, BArrayDense)(
    BDENSE_TYPE() && x
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

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator=)(
    BDENSE_TYPE() && x
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

BDENSE_TEMPLATE(bool, operator==) (
    const BDENSE_TYPE() & Array_
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

BDENSE_TEMPLATE(, ~BArrayDense) () {
    
    if (delete_data && (data != nullptr))
        delete data;
    
    return;
}

BDENSE_TEMPLATE(void, set_data) (
    Data_Type * data_,
    bool delete_data_
) {  

    if ((data != nullptr) && delete_data)
        delete data;
    
    data        = data_;
    delete_data = delete_data_;
    
    return;
    
}

BDENSE_TEMPLATE(Data_Type *, D_ptr) () {
    return this->data;
}

BDENSE_TEMPLATE(const Data_Type *, D_ptr) () const {
    return this->data;
}

BDENSE_TEMPLATE(Data_Type &, D) () {
    return *this->data;
}

BDENSE_TEMPLATE(const Data_Type &, D) () const {
    return *this->data;
}

BDENSE_TEMPLATE(void, out_of_range) (
    uint i,
    uint j
) const {

    if (i >= N)
        throw std::range_error("The row is out of range.");
    else if (j >= M)
        throw std::range_error("The column is out of range.");

    return;

}
    
BDENSE_TEMPLATE(Cell_Type, get_cell) (
    uint i,
    uint j,
    bool check_bounds
) const {
    
    // Checking boundaries  
    if (check_bounds)
        out_of_range(i,j);
    
    return el[POS(i, j)];
    
}

BDENSE_TEMPLATE(std::vector< Cell_Type >, get_row_vec) (
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    std::vector< Cell_Type > ans(ncol(), static_cast< Cell_Type >(false));
    for (uint j = 0u; j < M; ++j) 
        ans[j] = el[POS(i, j)];
    
    return ans;

}

BDENSE_TEMPLATE(void, get_row_vec) (
    std::vector<Cell_Type> * x,
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(i, 0u);

    for (uint j = 0u; j < M; ++j) 
        x->at(j) = el[POS(i, j)];
    
}

BDENSE_TEMPLATE(std::vector< Cell_Type >, get_col_vec)(
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    std::vector< Cell_Type > ans(nrow(), static_cast< Cell_Type >(false));
    for (uint j = 0u; j < N; ++j) 
        ans[j] = el[POS(j, i)];
    
    return ans;

}

BDENSE_TEMPLATE(void, get_col_vec) (
    std::vector<Cell_Type> * x,
    uint i,
    bool check_bounds
) const {

    // Checking boundaries  
    if (check_bounds) 
        out_of_range(0u, i);

    for (uint j = 0u; j < N; ++j) 
        x->at(j) = el[POS(j, i)];//this->get_cell(iter->first, i, false);
    
}
template<typename Cell_Type, typename Data_Type>
inline const BArrayDenseRow_const<Cell_Type,Data_Type> BDENSE_TYPE()::row(
    uint i,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, 0u);

    return BArrayDenseRow_const<Cell_Type,Data_Type>(*this, i);

}

template<typename Cell_Type, typename Data_Type>
inline BArrayDenseRow<Cell_Type,Data_Type> & BDENSE_TYPE()::row(
    uint i,
    bool check_bounds
) {

    if (check_bounds)
        out_of_range(i, 0u);

    return BArrayDenseRow<Cell_Type,Data_Type>(*this, i);

}

template<typename Cell_Type, typename Data_Type>
inline const BArrayDenseCol_const<Cell_Type,Data_Type> 
BArrayDense<Cell_Type,Data_Type>::col(
    uint j,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(0u, j);

    return BArrayDenseCol_const<Cell_Type,Data_Type>(*this, j);

}

template<typename Cell_Type, typename Data_Type>
inline BArrayDenseCol<Cell_Type,Data_Type> & 
BArrayDense<Cell_Type,Data_Type>::col(
    uint j,
    bool check_bounds
) {

    if (check_bounds)
        out_of_range(0u, j);

    return BArrayDenseCol<Cell_Type,Data_Type>(*this, j);

}

BDENSE_TEMPLATE(Entries< Cell_Type >, get_entries)() const {
    
    unsigned int nzero = this->nnozero();

    Entries<Cell_Type> res(nzero);
    
    for (uint i = 0u; i < N; ++i)
    {
        for (uint j = 0u; j < M; ++j)
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

BDENSE_TEMPLATE(bool, is_empty)(
    uint i,
    uint j,
    bool check_bounds
) const {
    
    if (check_bounds)
        out_of_range(i, j);
    
    return el[POS(i, j)] == ZERO_CELL;
    
}

BDENSE_TEMPLATE(unsigned int, nrow)() const noexcept {
    return N;
}

BDENSE_TEMPLATE(unsigned int, ncol)() const noexcept {
    return M;
}

BDENSE_TEMPLATE(unsigned int, nnozero)() const noexcept {

    unsigned int nzero = 0u;
    for (auto & v : el)
        if (v != BARRY_ZERO_DENSE)
            nzero++;

    return nzero;
}

BDENSE_TEMPLATE(Cell< Cell_Type>, default_val)() const {
    return this->Cell_default;
}

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator+=)(
    const std::pair<uint,uint> & coords
) {
    

    unsigned int i = coords.first;
    unsigned int j = coords.second;

    out_of_range(i, j);

    el[POS(i,j)]  += 1;
    el_rowsums[i] += 1;
    el_colsums[j] += 1;
    
    return *this;
    
}

BDENSE_TEMPLATE(BDENSE_TYPE() &, operator-=)(
    const std::pair<uint,uint> & coords
) {
    
    unsigned int i = coords.first;
    unsigned int j = coords.second;

    out_of_range(i, j);

    Cell_Type old = el[POS(i,j)];

    el[POS(i,j)]   = ZERO_CELL;
    el_rowsums[i] -= old;
    el_colsums[j] -= old;
    
    return *this;
    
}

template BDENSE_TEMPLATE_ARGS()
inline BArrayDenseCell<Cell_Type,Data_Type> BDENSE_TYPE()::operator()(  
    uint i,
    uint j,
    bool check_bounds
) {
    
    return BArrayDenseCell<Cell_Type,Data_Type>(this, i, j, check_bounds);
    
}

template BDENSE_TEMPLATE_ARGS()
inline const Cell_Type BDENSE_TYPE()::operator()(  
    uint i,
    uint j,
    bool check_bounds
) const {

    if (check_bounds)
        out_of_range(i, j);
    
    return el[POS(i,j)];
    
}

BDENSE_TEMPLATE(void, rm_cell) (
    uint i,
    uint j,
    bool check_bounds,
    bool check_exists
) {
    
    // Checking the boundaries
    if (check_bounds)
        out_of_range(i,j);

    BARRY_UNUSED(check_exists)
        
    // Remove the pointer first (so it wont point to empty)
    el_rowsums[i] -= el[POS(i, j)];
    el_colsums[j] -= el[POS(i, j)];    
    el[POS(i, j)] = BARRY_ZERO_DENSE;
    
    return;

}

BDENSE_TEMPLATE(void, insert_cell) (
    uint i,
    uint j,
    const Cell< Cell_Type> & v,
    bool check_bounds,
    bool check_exists
) { 
    
    if (check_bounds)
        out_of_range(i,j); 
    
    BARRY_UNUSED(check_exists)

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

BDENSE_TEMPLATE(void, insert_cell)(
    uint i,
    uint j,
    Cell_Type v,
    bool check_bounds,
    bool check_exists
) {
    
    if (check_bounds)
        out_of_range(i,j);

    BARRY_UNUSED(check_exists)
        
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

BDENSE_TEMPLATE(void, swap_cells) (
        uint i0, uint j0,
        uint i1, uint j1,
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

BDENSE_TEMPLATE(void, toggle_cell) (
    uint i,
    uint j,
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

BDENSE_TEMPLATE(void, swap_rows) (
    uint i0,
    uint i1,
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
    for (uint j = 0u; j < M; ++j)
        std::swap(el[POS(i0, j)], el[POS(i1, j)]);

    std::swap(el_rowsums[i0], el_rowsums[i1]);
    
    return;
}

// This swapping is more expensive overall
BDENSE_TEMPLATE(void, swap_cols) (
    uint j0,
    uint j1,
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
    for (uint i = 0u; i < N; ++i)
        std::swap(el[POS(i, j0)], el[POS(i, j1)]);

    std::swap(el_colsums[j0], el_colsums[j1]);
    
    return;
}

BDENSE_TEMPLATE(void, zero_row) (
    uint i,
    bool check_bounds
    ) {
  
    if (check_bounds)
        out_of_range(i, 0u);

    if (el_rowsums[i] == ZERO_CELL)
        return;

    // Else, remove all elements
    for (unsigned int col = 0u; col < M; col++) 
        rm_cell(i, col, false, false);
    
    return;
  
}

BDENSE_TEMPLATE(void, zero_col) (
    uint j,
    bool check_bounds
  ) {
  
    if (check_bounds)
        out_of_range(0u, j);
    
    if (el_colsums[j] == ZERO_CELL)
        return;
    
    // Else, remove all elements
    for (unsigned int row = 0u; row < N; row++) 
        rm_cell(row, j, false, false);
    
    return;
  
}

BDENSE_TEMPLATE(void, transpose) () {
  
    // if (NCells == 0u)
    // {

    //     std::swap(N, M);
    //     return;

    // }

    // Start by flipping the switch 
    visited = !visited;

    // uint N0 = N, M0 = M;
    std::vector< Cell< Cell_Type > > tmp_el(std::move(el));
    el.resize(N * M, ZERO_CELL);
    for (uint i = 0u; i < N; ++i) 
        for (uint j = 0u; j < M; ++j)
            std::swap(tmp_el[POS(i, j)], el[POS_N(j, i, M)]);
    
    // Swapping the values
    std::swap(N, M);
    std::swap(el_rowsums, el_colsums);
    
    return;

}

BDENSE_TEMPLATE(void, clear) (
    bool hard
) {
    
    BARRY_UNUSED(hard)
    
    std::fill(el.begin(), el.end(), ZERO_CELL);
    std::fill(el_rowsums.begin(), el_rowsums.end(), ZERO_CELL);
    std::fill(el_colsums.begin(), el_colsums.end(), ZERO_CELL);
    
    return;
    
}

BDENSE_TEMPLATE(void, resize) (
    uint N_,
    uint M_
) {

    // Moving stuff around
    std::vector< Cell_Type > el_tmp(el);
    el.resize(N_ * M_, ZERO_CELL);
    el_rowsums.resize(N_, ZERO_CELL);
    el_colsums.resize(M_, ZERO_CELL);

    for (unsigned int i = 0u; i < N; ++i)
    {
        // If reached the end
        if (i >= N_)
            break;

        for (unsigned int j = 0u; j < M; ++j)
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

BDENSE_TEMPLATE(void, reserve) () {

    el.reserve(N * M);
    el_rowsums.reserve(N);
    el_colsums.reserve(M);
    return;
  
}

BDENSE_TEMPLATE(void, print) (
    const char * fmt,
    ...
) const
{
  
    std::va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);

    for (uint i = 0u; i < N; ++i)
    {

        printf_barry("[%3i,] ", i);

        for (uint j = 0u; j < M; ++j)
        {

            if (this->is_empty(i, j, false))
                printf_barry("    . ");
            else 
                printf_barry(" %.2f ", static_cast<double>(this->get_cell(i, j, false)));
            
        }

        printf_barry("\n");

    }
    
    return;
    
}

BDENSE_TEMPLATE(const std::vector< Cell_Type > &, get_data)() const
{
    return el;
}

BDENSE_TEMPLATE(const Cell_Type, rowsum)(unsigned int i) const
{
    return el_rowsums[i];
}

BDENSE_TEMPLATE(const Cell_Type, colsum)(unsigned int j) const
{
    return el_colsums[j];
}

#undef ROW
#undef COL
#undef POS
#undef POS_N

#undef BDENSE_TYPE
#undef BDENSE_TEMPLATE_ARGS
#undef BDENSE_TEMPLATE
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
    Cell_Type old      =  dat->el[POS(i,j)];
    dat->el[POS(i,j)]  =  val;
    dat->el_rowsums[i] += (val - old);
    dat->el_colsums[j] += (val - old);

    return *this;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator=(const Cell_Type & val) {

    Cell_Type old      =  dat->el[POS(i,j)];
    dat->el[POS(i,j)]  =  val;
    dat->el_rowsums[i] += (val - old);
    dat->el_colsums[j] += (val - old);
    
}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator+=(const Cell_Type & val) {
    
    dat->el[POS(i,j)]  += val;
    dat->el_rowsums[i] += val;
    dat->el_colsums[j] += val;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator-=(const Cell_Type & val) {
    
    dat->el[POS(i,j)]  -= val;
    dat->el_rowsums[i] -= val;
    dat->el_colsums[j] -= val;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator*=(const Cell_Type & val) {
    
    Cell_Type old = dat->el[POS(i,j)];
    dat->el_colsums[j] += (old * val - old);
    dat->el_rowsums[i] += (old * val - old);
    dat->el[POS(i,j)] *= val;

}

template<typename Cell_Type,typename Data_Type>
inline void BArrayDenseCell<Cell_Type,Data_Type>::operator/=(const Cell_Type & val) {
    
    Cell_Type old = dat->el[POS(i,j)];
    dat->el_rowsums[i] += (old/val - old);
    dat->el_colsums[j] += (old/val - old);
    dat->el[POS(i,j)]  /= val;

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
    
    for (uint i = 0u; i < nrow(); ++i)
        for (uint j = 0u; j < ncol(); ++j)
            this->operator()(i, j) += rhs.get_cell(i, j);

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator+=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
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
    
    for (uint i = 0u; i < nrow(); ++i) {
        for (uint j = 0u; j < ncol(); ++j) {
            this->operator()(i, j) -= rhs.get_cell(i, j);
        }
    }

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator-=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) 
        for (uint j = 0u; j < ncol(); ++j) 
            this->operator()(i, j) -= rhs;
        
    

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator*=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) 
        for (uint j = 0u; j < nrow(); ++j)
            el[POS(i, j)] *= rhs;

    return *this;
}

BDENSE_TEMPLATE(BDENSE_TYPE()&, operator/=) (
    const Cell_Type& rhs
) {

    for (uint i = 0u; i < nrow(); ++i) 
        for (uint j = 0u; j < nrow(); ++j)
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
    double count(Array_Type & Array, uint i, uint j);
    double init(Array_Type & Array, uint i, uint j);
    std::string get_name() const;
    std::string get_description() const;

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
    Counter<Array_Type,Data_Type> & operator[](uint idx);

    /**
     * @brief Number of counters in the set.
     * 
     * @return uint 
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

    void add_hash(
      Hasher_fun_type<Array_Type,Data_Type> fun_
    );
    
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

COUNTER_TEMPLATE(double, count)(Array_Type & Array, uint i, uint j)
{

    if (count_fun == nullptr)
        return 0.0;

    return count_fun(Array, i, j, data);

}

COUNTER_TEMPLATE(double, init)(Array_Type & Array, uint i, uint j)
{

    if (init_fun == nullptr)
        return 0.0;

    return init_fun(Array, i, j, data);

}

COUNTER_TEMPLATE(std::string, get_name)() const {
    return this->name;
}

COUNTER_TEMPLATE(std::string, get_description)() const {
    return this->name;
}

COUNTER_TEMPLATE(void, set_hasher)(Hasher_fun_type<Array_Type,Data_Type> fun) {
    hasher_fun = fun;
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

COUNTERS_TEMPLATE(COUNTER_TYPE() &, operator[])(uint idx) {

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
  
    data.push_back(Counter<Array_Type,Data_Type>(
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

    std::vector< std::string > out(this->size());
    for (unsigned int i = 0u; i < out.size(); ++i)
        out[i] = this->data.at(i).get_name();

    return out;

}

COUNTERS_TEMPLATE(std::vector<std::string>, get_descriptions)() const
{
    
    std::vector< std::string > out(this->size());
    for (unsigned int i = 0u; i < out.size(); ++i)
        out[i] = data.at(i).get_description();

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
    void count_init(uint i, uint j);
    void count_current(uint i, uint j);
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

STATSCOUNTER_TEMPLATE(void, count_init)(uint i,uint j)
{
    
    // Do we have any counter?
    if (counters->size() == 0u)
        throw std::logic_error("No counters added: Cannot count without knowning what to count!");
    
    // Iterating through the functions, and updating the set of
    // statistics.
    current_stats.resize(counters->size(), 0.0);
    // change_stats.resize(counters->size(), 0.0);
    for (uint n = 0u; n < counters->size(); ++n) 
        current_stats[n] = counters->operator[](n).init(EmptyArray, i, j);
    
    return;
}

STATSCOUNTER_TEMPLATE(void, count_current)(uint i, uint j)
{
    
    // Iterating through the functions, and updating the set of
    // statistics.
    for (uint n = 0u; n < counters->size(); ++n) {
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
    for (uint i = 0; i < Array->nrow(); ++i)
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
    for (unsigned int i = 0u; i < Array->nrow(); ++i)
    {

        for (unsigned int j = 0u; j < Array->ncol(); ++j)
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
        uint pos = 0u,
        std::vector< Array_Type > * array_bank = nullptr,
        std::vector< double > * stats_bank = nullptr
    );

    void calc_backend_dense(
        uint pos = 0u,
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
    
public:
    
    uint N, M;
    bool delete_counters  = true;
    bool delete_rules     = true;
    bool delete_rules_dyn = true;
    uint max_num_elements = BARRY_MAX_NUM_ELEMENTS;
    
    // Temp variables to reduce memory allocation
    std::vector< double >                current_stats;
    std::vector< size_t >                coordinates_free;
    std::vector< size_t >                coordinates_locked;
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
    Support(uint N_, uint M_) :
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
    bool eval_rules_dyn(const std::vector<double> & counts, const uint & i, const uint & j);
    // bool eval_rules_dyn(const double * counts, const uint & i, const uint & j);
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
        unsigned int max_num_elements_ = 0u
    );
    
    std::vector< double > get_counts() const;
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

#define SUPPORT_TEMPLATE_ARGS() <typename Array_Type, typename \
    Data_Counter_Type, typename Data_Rule_Type, typename Data_Rule_Dyn_Type>

#define SUPPORT_TYPE() Support<Array_Type,Data_Counter_Type,Data_Rule_Type,\
    Data_Rule_Dyn_Type>

#define SUPPORT_TEMPLATE(a,b) template SUPPORT_TEMPLATE_ARGS() \
    inline a SUPPORT_TYPE()::b

SUPPORT_TEMPLATE(void, init_support)(
    std::vector< Array_Type > * array_bank,
    std::vector< double > * stats_bank
) {
    
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

        for (uint i = 0u; i < coordiantes_n_free; ++i)
            EmptyArray.rm_cell(
                coordinates_free[i * 2u],
                coordinates_free[i * 2u + 1u],
                false, true
                );
                
    }

    // Looked coordinates should still be removed if these are
    // equivalent to zero
    for (unsigned int i = 0u; i < coordiantes_n_locked; ++i)
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
        for (uint n = 0u; n < n_counters; ++n)
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

SUPPORT_TEMPLATE(void, reset_array)() {
    
    data.clear();
    
}

SUPPORT_TEMPLATE(void, reset_array)(const Array_Type & Array_) {
    
    data.clear();
    EmptyArray = Array_;
    N = Array_.nrow();
    M = Array_.ncol();
    
}

SUPPORT_TEMPLATE(void, calc_backend_sparse)(
        uint pos,
        std::vector< Array_Type > * array_bank,
        std::vector< double > * stats_bank
    ) {
    
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
    unsigned int change_stats_different = hashes_initialized[pos] ? 0u : 1u;
    for (uint n = 0u; n < n_counters; ++n)
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
        #ifdef __OPENMP
        #pragma omp simd
        #else
        #pragma GCC ivdep
        #endif
        for (uint n = 0u; n < n_counters; ++n) 
            current_stats[n] -= change_stats[pos * n_counters + n];
    }
        
    
    return;
    
}

SUPPORT_TEMPLATE(void, calc_backend_dense)(
        uint pos,
        std::vector< Array_Type > * array_bank,
        std::vector< double > * stats_bank
    ) {
    
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
    unsigned int change_stats_different = hashes_initialized[pos] ? 0u : 1u;
    for (uint n = 0u; n < n_counters; ++n)
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
        #ifdef __OPENMP
        #pragma omp simd
        #else
        #pragma GCC ivdep
        #endif
        for (uint n = 0u; n < n_counters; ++n) 
            current_stats[n] -= change_stats[pos * n_counters + n];
    }
    
    return;
    
}

SUPPORT_TEMPLATE(void, calc)(
        std::vector< Array_Type > * array_bank,
        std::vector< double > * stats_bank,
        unsigned int max_num_elements_
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

SUPPORT_TEMPLATE(void, add_counter)(
        Counter<Array_Type,Data_Counter_Type> f_
) {
    
    counters->add_counter(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, set_counters)(
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

SUPPORT_TEMPLATE(void, add_rule)(
        Rule<Array_Type, Data_Rule_Type> * f_
) {
    
    rules->add_rule(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, add_rule)(
        Rule<Array_Type,Data_Rule_Type> f_
) {
    
    rules->add_rule(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, set_rules)(
        Rules<Array_Type,Data_Rule_Type> * rules_
) {
    
    // Cleaning up before replacing the memory
    if (delete_rules)
        delete rules;
    delete_rules = false;
    rules = rules_;
    
    return;
    
}

SUPPORT_TEMPLATE(void, add_rule_dyn)(
        Rule<Array_Type, Data_Rule_Dyn_Type> * f_
) {
    
    rules_dyn->add_rule(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, add_rule_dyn)(
        Rule<Array_Type,Data_Rule_Dyn_Type> f_
) {
    
    rules_dyn->add_rule(f_);
    return;
    
}

SUPPORT_TEMPLATE(void, set_rules_dyn)(
        Rules<Array_Type,Data_Rule_Dyn_Type> * rules_
) {
    
    // Cleaning up before replacing the memory
    if (delete_rules_dyn)
        delete rules_dyn;
    delete_rules_dyn = false;
    rules_dyn = rules_;
    
    return;
    
}

SUPPORT_TEMPLATE(bool, eval_rules_dyn)(
    const std::vector< double > & counts,
    const uint & i,
    const uint & j
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

// SUPPORT_TEMPLATE(bool, eval_rules_dyn)(
//     const double * counts,
//     const uint & i,
//     const uint & j
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

SUPPORT_TEMPLATE(std::vector< double >, get_counts)() const {
    
    return data.get_data(); 
    
}

// SUPPORT_TEMPLATE(const MapVec_type<> *, get_counts_ptr)() const {
    
//     return data.get_data_ptr();
      
// }

SUPPORT_TEMPLATE(std::vector< double > *, get_current_stats)() {
    return &this->current_stats;
}

SUPPORT_TEMPLATE(void, print)() const {

    // Starting from the name of the stats
    printf_barry("Position of variables:\n");
    for (uint i = 0u; i < n_counters; ++i) {
        printf_barry("[% 2i] %s\n", i, counters->operator[](i).name.c_str());
    }

    data.print();
}

SUPPORT_TEMPLATE(const FreqTable<double> &, get_data)() const {
    return this->data;
}

template SUPPORT_TEMPLATE_ARGS()
inline Counters<Array_Type,Data_Counter_Type> * SUPPORT_TYPE()::get_counters() {
    return this->counters;
}   
    
template SUPPORT_TEMPLATE_ARGS()
inline Rules<Array_Type,Data_Rule_Type> * SUPPORT_TYPE()::get_rules() {
    return this->rules;
}

template SUPPORT_TEMPLATE_ARGS()
inline Rules<Array_Type,Data_Rule_Dyn_Type> * SUPPORT_TYPE()::get_rules_dyn() {
    return this->rules_dyn;
}

#undef SUPPORT_TEMPLATE_ARGS
#undef SUPPORT_TYPE
#undef SUPPORT_TEMPLATE

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
    void calc_backend_sparse(uint pos = 0u);  
    void calc_backend_dense(uint pos = 0u);  

public:
    Array_Type                         EmptyArray;
    std::vector< Array_Type >          data;
    Rules<Array_Type,Data_Rule_Type> * rules;

    uint N, M;
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
    PowerSet(uint N_, uint M_) :
        EmptyArray(N_, M_), data(0u),
        rules(new Rules<Array_Type,Data_Rule_Type>()), N(N_), M(M_) {};
    PowerSet(const Array_Type & array);

    ~PowerSet();
    ///@}
    
    void init_support();
    void calc();
    void reset(uint N_, uint M_);
    
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
    const Array_Type& operator[](const unsigned int & i) const {return data.at(i);};
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

            for (uint i = 0u; i < n_free; ++i) 
                EmptyArray(
                    coordinates_free[i * 2u],
                    coordinates_free[i * 2u + 1u]
                    ) = 0;

        }
        else
        {

            for (uint i = 0u; i < n_free; ++i) 
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
    uint pos
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
    uint pos
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
        uint N_,
        uint M_
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

private:
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
    std::vector< std::vector< double > > stats_support;          ///< Sufficient statistics of the model (support)
    std::vector< uint >                  stats_support_n_arrays; ///< Number of arrays included per support.
    std::vector< std::vector< double > > stats_target;           ///< Target statistics of the model
    std::vector< uint >                  arrays2support;
    ///@}

    /**
      * @brief Map of types of arrays to support sets
      * @details This is of the same length as the vector `stats_target`.
      */
    MapVec_type< double, uint > keys2support;

    /**
     * @name Container space for the powerset (and its sufficient stats_target)
     * @details This is useful in the case of using simulations or evaluating
     * functions that need to account for the full set of states.
     */
    ///@{
    bool with_pset = false;
    std::vector< std::vector< Array_Type > > pset_arrays; ///< Arrays of the support(s)
    std::vector< std::vector<double> >       pset_stats;  ///< Statistics of the support(s)
    std::vector< std::vector<double> >       pset_probs;  ///< Probabilities of the support(s)
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
     * - `k` unsigned int indicating the number of sufficient statistics
     * 
     * @returns
     * Nothing, but it will modify the model data.
     */
    std::function<std::vector<double>(double *, unsigned int k)>
        transform_model_fun = nullptr;

    std::vector< std::string > transform_model_term_names;
    
public:
    
    void set_rengine(std::mt19937 * rengine_, bool delete_ = false) {

        if (delete_rengine)
            delete rengine;

        rengine        = rengine_;
        delete_rengine = delete_;
        
    };

    void set_seed(unsigned int s) {

        if (rengine == nullptr)
        {
            rengine = new std::mt19937;
            delete_rengine = true;
        }

        rengine->seed(s);

    };
    ///@}
        
    Model();
    Model(uint size_);
    Model(const Model<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> & Model_);
    Model<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> & operator=(
        const Model<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> & Model_
    );

    ~Model() {
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
    uint add_array(const Array_Type & Array_, bool force_new = false);
    
    
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
        const uint & i,
        bool as_log = false
    );
    
    double likelihood(
        const std::vector<double> & params,
        const Array_Type & Array_,
        int i = -1,
        bool as_log = false
    );
    
    double likelihood(
        const std::vector<double> & params,
        const std::vector<double> & target_,
        const uint & i,
        bool as_log = false
    );

    double likelihood(
        const std::vector<double> & params,
        const double * target_,
        const uint & i,
        bool as_log = false
    );
    
    double likelihood_total(
        const std::vector<double> & params,
        bool as_log = false
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
    double get_norm_const(
        const std::vector< double > & params,
        const uint & i,
        bool as_log = false
    );

    const std::vector< Array_Type > * get_pset(
        const uint & i
    );

    const std::vector< double > * get_pset_stats(
        const uint & i
    );
    ///@}
    
    void print_stats(uint i) const;

    /**
     * @brief Prints information about the model
     */
    void print() const;
    
    Array_Type sample(const Array_Type & Array_, const std::vector<double> & params = {});
    Array_Type sample(const uint & i, const std::vector<double> & params);
    
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
        unsigned int i,
        unsigned int j
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
    unsigned int size() const noexcept;
    unsigned int size_unique() const noexcept;
    unsigned int nterms() const noexcept;
    unsigned int support_size() const noexcept;
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
    std::vector< std::vector< double > > * get_stats_support();
    std::vector< unsigned int > * get_arrays2support();
    std::vector< std::vector< Array_Type > > * get_pset_arrays();
    std::vector< std::vector<double> > * get_pset_stats();  ///< Statistics of the support(s)
    std::vector< std::vector<double> > * get_pset_probs(); 
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
        std::function<std::vector<double>(double*,unsigned int)> fun,
        std::vector< std::string > names
        );
    std::vector<double> transform_model(
        double * data,
        unsigned int k
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
    const double * params,
    const double * support,
    size_t k,
    size_t n
)
{
    
    double res = 0.0;
    
    #ifdef __OPENMP
    #pragma omp simd reduction(+:res) 
    #else
    #pragma GCC ivdep
    #endif
    for (unsigned int i = 0u; i < n; ++i)
    {

        double tmp = 0.0;
        const double * support_n = support + i * k + 1u;
        
        for (unsigned int j = 0u; j < (k - 1u); ++j)
            tmp += (*(support_n + j)) * (*(params + j));
        
        res += std::exp(tmp BARRY_SAFE_EXP) * (*(support + i * k));

    }
    
    // This will only evaluate if the option BARRY_CHECK_FINITE
    // is defined
    BARRY_ISFINITE(res)

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
    #ifdef __OPENMP
    #pragma omp simd reduction(+:numerator)
    #else
    #pragma GCC ivdep
    #endif
    for (unsigned int j = 0u; j < params.size(); ++j)
        numerator += *(stats_target + j) * params[j];

    if (!log_)
        numerator = exp(numerator BARRY_SAFE_EXP);
    else
        return numerator BARRY_SAFE_EXP - log(normalizing_constant);

    double ans = numerator/normalizing_constant;

    if (ans > 1.0)
        printf_barry("ooo\n");

    return ans;
    
}

#define MODEL_TYPE() Model<Array_Type,Data_Counter_Type,Data_Rule_Type,\
    Data_Rule_Dyn_Type>

#define MODEL_TEMPLATE_ARGS() <typename Array_Type, typename Data_Counter_Type,\
    typename Data_Rule_Type, typename Data_Rule_Dyn_Type>

#define MODEL_TEMPLATE(a,b) \
    template MODEL_TEMPLATE_ARGS() inline a MODEL_TYPE()::b


MODEL_TEMPLATE(,Model)() :
    stats_support(0u),
    stats_support_n_arrays(0u),
    stats_target(0u), arrays2support(0u),
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

MODEL_TEMPLATE(,Model)(uint size_) :
    stats_support(0u),
    stats_support_n_arrays(0u),
    stats_target(0u), arrays2support(0u), keys2support(0u), 
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

MODEL_TEMPLATE(,Model)(
    const MODEL_TYPE() & Model_) : 
    stats_support(Model_.stats_support),
    stats_support_n_arrays(Model_.stats_support_n_arrays),
    stats_target(Model_.stats_target),
    arrays2support(Model_.arrays2support),
    keys2support(Model_.keys2support),
    pset_arrays(Model_.pset_arrays),
    pset_stats(Model_.pset_stats),
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

MODEL_TEMPLATE(MODEL_TYPE() &, operator)=(
    const MODEL_TYPE() & Model_
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
        stats_support_n_arrays     = Model_.stats_support_n_arrays;
        stats_target               = Model_.stats_target;
        arrays2support             = Model_.arrays2support;
        keys2support               = Model_.keys2support;
        pset_arrays                = Model_.pset_arrays;
        pset_stats                 = Model_.pset_stats;
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

MODEL_TEMPLATE(void, store_psets)() noexcept {
    // if (with_pset)
    //   throw std::logic_error("Powerset storage alreay activated.");
    with_pset = true;
    return;
}

MODEL_TEMPLATE(std::vector< double >, gen_key)(
    const Array_Type & Array_
) {
    return this->counters->gen_hash(Array_);   
}

MODEL_TEMPLATE(void, add_counter)(
        Counter<Array_Type, Data_Counter_Type> & counter
) {
    
    counters->add_counter(counter, Data_Counter_Type());
    return;
}

MODEL_TEMPLATE(void, add_counter)(
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

MODEL_TEMPLATE(void, set_counters)(
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

MODEL_TEMPLATE(void, add_hasher)(
    Hasher_fun_type<Array_Type,Data_Counter_Type> fun_
) {

    counters->add_hash(fun_);

}

////////////////////////////////////////////////////////////////////////////////

MODEL_TEMPLATE(void, add_rule)(
    Rule<Array_Type, Data_Rule_Type> & rules
) {
    
    rules->add_rule(rules, Data_Rule_Type());
    return;
}


MODEL_TEMPLATE(void, set_rules)(
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

MODEL_TEMPLATE(void, add_rule_dyn)(
    Rule<Array_Type, Data_Rule_Dyn_Type> & rules_
) {
    
    rules_dyn->add_rule(rules_, Data_Rule_Dyn_Type());
    return;
}

MODEL_TEMPLATE(void, add_rule_dyn)(
    Rule_fun_type<Array_Type,Data_Rule_Dyn_Type> rule_fun_,
    Data_Rule_Dyn_Type                           data_
) {
    
    rules_dyn->add_rule(
        rule_fun_,
        data_
    );
    
    return;
    
}

MODEL_TEMPLATE(void, set_rules_dyn)(
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

MODEL_TEMPLATE(uint, add_array)(
    const Array_Type & Array_,
    bool force_new
) {
    
    // Array counts (target statistics)
    counter_fun.reset_array(&Array_);
    
    if (transform_model_fun)
    {
        auto tmpcounts = counter_fun.count_all();
        stats_target.push_back(transform_model_fun(&tmpcounts[0u], tmpcounts.size()));
    } else
        stats_target.push_back(counter_fun.count_all());
    
    // If the data hasn't been analyzed earlier, then we need to compute
    // the support
    std::vector< double > key = counters->gen_hash(Array_);
    MapVec_type< double, uint >::const_iterator locator = keys2support.find(key);
    if (force_new | (locator == keys2support.end()))
    {
        
        // Adding to the map
        keys2support[key] = stats_support.size();
        stats_support_n_arrays.push_back(1u);       // How many elements now
        arrays2support.push_back(stats_support.size()); // Map of the array id to the support
        
        // Computing support using the counters included in the model
        support_fun.reset_array(Array_);
        
        /** When computing with the powerset, we need to grow the corresponding
            * vectors on the fly */
        if (with_pset)
        {
            
            // Making space for storing the support
            pset_arrays.resize(pset_arrays.size() + 1u);
            pset_stats.resize(pset_stats.size() + 1u);
            pset_probs.resize(pset_probs.size() + 1u);
            
            try
            {
                
                support_fun.calc(
                    &(pset_arrays[pset_arrays.size() - 1u]),
                    &(pset_stats[pset_stats.size() - 1u])
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

            stats_support.push_back(s_new);

        } else 
            stats_support.push_back(support_fun.get_counts());
        
        // Making room for the previous parameters. This will be used to check if
        // the normalizing constant has been updated or not.
        params_last.push_back(stats_target[0u]);
        normalizing_constants.push_back(0.0);
        first_calc_done.push_back(false);
        
        return arrays2support.size() - 1u;
        
    }
    
    // Increasing the number of arrays in that stat
    ++stats_support_n_arrays[locator->second];
    
    // Adding the corresponding map
    arrays2support.push_back(locator->second);
    
    return arrays2support.size() - 1u;

}

MODEL_TEMPLATE(double, likelihood)(
    const std::vector<double> & params,
    const uint & i,
    bool as_log
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    unsigned int idx = arrays2support[i];

    // Checking if this actually has a change of happening
    if (this->stats_support[idx].size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[idx] || !vec_equal_approx(params, params_last[idx]) )
    {
        
        first_calc_done[idx] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support[idx].size() / k;

        normalizing_constants[idx] = update_normalizing_constant(
            &params[0u], &stats_support[idx][0u], k, n
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

MODEL_TEMPLATE(double, likelihood)(
    const std::vector<double> & params,
    const Array_Type & Array_,
    int i,
    bool as_log
) {
    
    // Key of the support set to use
    int loc;

    if (i < 0)
    {

        std::vector< double > key = counters->gen_hash(Array_);
        MapVec_type< double, uint >::const_iterator locator = keys2support.find(key);
        if (locator == keys2support.end()) 
            throw std::range_error("This type of array has not been included in the model.");

        loc = locator->second;

    }
    else
    {

        if (static_cast<uint>(i) >= arrays2support.size())
            throw std::range_error("This type of array has not been included in the model.");

        loc = arrays2support[i];

    }

    // Checking if this actually has a change of happening
    if (this->stats_support[loc].size() == 0u)
        return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
    
    // Counting stats_target
    StatsCounter< Array_Type, Data_Counter_Type> tmpstats(&Array_);

    tmpstats.set_counters(this->counters);
    
    std::vector< double > target_ = tmpstats.count_all();

    if (transform_model_fun)
        target_ = transform_model_fun(&target_[0u], target_.size());

    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) )
    {
        
        first_calc_done[loc] = true;

        size_t k = params.size() + 1u;
        size_t n = stats_support[loc].size() / k;
        
        normalizing_constants[loc] = update_normalizing_constant(
            &params[0u], &stats_support[loc][0u], k, n
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

MODEL_TEMPLATE(double, likelihood)(
    const std::vector<double> & params,
    const std::vector<double> & target_,
    const uint & i,
    bool as_log
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    uint loc = arrays2support[i];

    // Checking if passes the rules
    if (!support_fun.eval_rules_dyn(target_, 0u, 0u))
    {
        throw std::range_error("The array is not in the support set.");
    }
        

    // Checking if this actually has a change of happening
    if (this->stats_support[loc].size() == 0u)
    {
        // return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
        throw std::logic_error("The support set for this array is empty.");
    }
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) ) {
        
        first_calc_done[loc] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support[loc].size() / k;

        normalizing_constants[loc] = update_normalizing_constant(
            &params[0u], &stats_support[loc][0u], k, n
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

MODEL_TEMPLATE(double, likelihood)(
    const std::vector<double> & params,
    const double * target_,
    const uint & i,
    bool as_log
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    uint loc = arrays2support[i];

    // Checking if passes the rules
    if (support_fun.get_rules_dyn()->size() > 0u)
    {

        std::vector< double > tmp_target(nterms(), 0.0);
        for (size_t t = 0u; t < nterms(); ++t)
            tmp_target[t] = *(target_ + t);

        if (!support_fun.eval_rules_dyn(tmp_target, 0u, 0u))
        {
            throw std::range_error("The array is not in the support set.");
            // return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
        }

    }

    // Checking if this actually has a change of happening
    if (this->stats_support[loc].size() == 0u)
    {
        // return as_log ? -std::numeric_limits<double>::infinity() : 0.0;
        throw std::logic_error("The support set for this array is empty.");
    }
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[loc] || !vec_equal_approx(params, params_last[loc]) ) {
        
        first_calc_done[loc] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support[loc].size() / k;

        normalizing_constants[loc] = update_normalizing_constant(
            &params[0u], &stats_support[loc][0u], k, n
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

MODEL_TEMPLATE(double, likelihood_total)(
    const std::vector<double> & params,
    bool as_log
) {
    
    size_t params_last_size = params_last.size();

    for (uint i = 0u; i < params_last_size; ++i)
    {

        if (!first_calc_done[i] || !vec_equal_approx(params, params_last[i]) )
        {

            size_t k = params.size() + 1u;
            size_t n = stats_support[i].size() / k;
            
            first_calc_done[i] = true;
            normalizing_constants[i] = update_normalizing_constant(
                &params[0u], &stats_support[i][0u], k, n
            );
            
            params_last[i] = params;
            
        }

    }
    
    double res = 0.0;
    if (as_log)
    {

        for (uint i = 0; i < stats_target.size(); ++i) 
            res += vec_inner_prod(
                &stats_target[i][0u],
                &params[0u],
                params.size()
                ) BARRY_SAFE_EXP;
        
        #ifdef __OPENM 
        #pragma omp simd reduction(-:res)
        #else
        #pragma GCC ivdep
        #endif
        for (unsigned int i = 0u; i < params_last_size; ++i)
            res -= (std::log(normalizing_constants[i]) * this->stats_support_n_arrays[i]);

    } else {
        
        res = 1.0;
        size_t stats_target_size = stats_target.size();
        #ifdef __OPENM 
        #pragma omp simd reduction(*:res)
        #else
        #pragma GCC ivdep
        #endif
        for (unsigned int i = 0; i < stats_target_size; ++i)
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

MODEL_TEMPLATE(double, get_norm_const)(
    const std::vector<double> & params,
    const uint & i,
    bool as_log
) {
    
    // Checking if the index exists
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    const auto id = arrays2support[i];
    
    // Checking if we have updated the normalizing constant or not
    if (!first_calc_done[id] || !vec_equal_approx(params, params_last[id]) )
    {
        
        first_calc_done[id] = true;
        
        size_t k = params.size() + 1u;
        size_t n = stats_support[id].size() / k;

        normalizing_constants[id] = update_normalizing_constant(
            &params[0u], &stats_support[id][0u], k, n
        );
        
        params_last[id] = params;
        
    }
    
    return as_log ? 
        std::log(normalizing_constants[id]) :
        normalizing_constants[id]
        ;
    
}

MODEL_TEMPLATE(const std::vector< Array_Type > *, get_pset)(
    const uint & i
) {

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");


    return &pset_arrays[arrays2support[i]];

}

MODEL_TEMPLATE(const std::vector< double > *, get_pset_stats)(
    const uint & i
) {

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    return &pset_stats[arrays2support[i]];

}

MODEL_TEMPLATE(void, print_stats)(uint i) const
{
    
    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    const auto & S = stats_support[arrays2support[i]];

    size_t k       = nterms();
    size_t nunique = S.size() / (k + 1u);

    for (uint l = 0u; l < nunique; ++l)
    {

        printf_barry("% 5i ", l);

        printf_barry("counts: %.0f motif: ", S[l * (k + 1u)]);
        
        for (unsigned int j = 0u; j < k; ++j)
            printf_barry("%.2f, ", S[l * (k + 1) + j + 1]);

        printf_barry("\n");

    }
    
    return;
    
}

MODEL_TEMPLATE(void, print)() const
{

    // Relevant information:
    // - Number of arrays involved
    // - Size of the support
    // - Terms involved

    int min_v = std::numeric_limits<int>::max();
    int max_v = 0;

    for (const auto & stat : this->stats_support)
    {

        if (static_cast<int>(stat.size()) > max_v)
            max_v = static_cast<int>(stat.size());
        
        if (static_cast<int>(stat.size()) < min_v)
            min_v = static_cast<int>(stat.size());

    }  

    // The vectors in the support reflec the size of nterms x entries
    max_v /= static_cast<int>(nterms() + 1);
    min_v /= static_cast<int>(nterms() + 1);

    printf_barry("Num. of Arrays     : %i\n", this->size());
    printf_barry("Support size       : %i\n", this->size_unique());
    printf_barry("Support size range : [%i, %i]\n", min_v, max_v);
    printf_barry("Transform. Fun.    : %s\n", transform_model_fun ? "yes": "no");
    printf_barry("Model terms (%i)   :\n", this->nterms());
    
    for (auto & cn : this->colnames())
        printf_barry(" - %s\n", cn.c_str());

    return;

}

MODEL_TEMPLATE(uint, size)() const noexcept
{
    // INITIALIZED()
    return this->stats_target.size();

}

MODEL_TEMPLATE(uint, size_unique)() const noexcept
{

    // INITIALIZED()
    return this->stats_support.size();

} 

MODEL_TEMPLATE(uint, nterms)() const noexcept
{
 
    if (transform_model_fun)
        return transform_model_term_names.size();
    else
        return this->counters->size();

}

MODEL_TEMPLATE(uint, support_size)() const noexcept
{

    // INITIALIZED()
    uint tot = 0u;
    for (auto& a : stats_support)
        tot += a.size();

    return tot;

}

MODEL_TEMPLATE(std::vector< std::string >, colnames)() const
{
    
    if (transform_model_fun)
        return transform_model_term_names;
    else
        return counters->get_names();

}
    
MODEL_TEMPLATE(Array_Type, sample)(
    const unsigned int & i,
    const std::vector<double> & params
) {

    // Are we recording this?
    if (!this->with_pset)
        throw std::logic_error("Sampling is only available when store_pset() is active.");

    if (i >= arrays2support.size())
        throw std::range_error("The requested support is out of range");

    // Getting the index
    unsigned int a = arrays2support[i];
    
    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    size_t k = params.size();

    // Sampling an array
    unsigned int j = 0u;
    std::vector< double > & probs = pset_probs[a];
    if ((probs.size() > 0u) && (vec_equal_approx(params, params_last[a])))
    // If precomputed, then no need to recalc support
    {

        while (cumprob < r)
            cumprob += probs[j++];

        j--;

    } else { 
       
        probs.resize(pset_arrays[a].size());
        std::vector< double > temp_stats(params.size());
        const std::vector< double > & stats = pset_stats[a];

        int i_matches = -1;
        for (size_t array = 0u; array < probs.size(); ++array)
        {

            // Filling out the parameters
            for (auto p = 0u; p < params.size(); ++p)
                temp_stats[p] = stats[array * k + p];

            probs[array] = this->likelihood(params, temp_stats, i, false);
            cumprob += probs[array];

            if (i_matches == -1 && cumprob >= r)
                i_matches = array;
        }

        j = i_matches;
        
    }
    

    return this->pset_arrays[a][j];   

}

MODEL_TEMPLATE(Array_Type, sample)(
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
    MapVec_type< double, uint >::const_iterator locator = keys2support.find(key);
    if (locator == keys2support.end())
    {
        // throw std::out_of_range("Sampling from an array that has no support in the model.");

        // Adding to the map
        keys2support[key] = stats_support.size();
        stats_support_n_arrays.push_back(1u);       // How many elements now
        arrays2support.push_back(stats_support.size()); // Map of the array id to the support
        
        // Computing support using the counters included in the model
        support_fun.reset_array(Array_);
        
        /** When computing with the powerset, we need to grow the corresponding
            * vectors on the fly */
        if (with_pset)
        {
            
            // Making space for storing the support
            pset_arrays.resize(pset_arrays.size() + 1u);
            pset_stats.resize(pset_stats.size() + 1u);
            pset_probs.resize(pset_probs.size() + 1u);
            
            try
            {
                
                support_fun.calc(
                    &(pset_arrays[pset_arrays.size() - 1u]),
                    &(pset_stats[pset_stats.size() - 1u])
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

            stats_support.push_back(s_new);

        } else 
            stats_support.push_back(support_fun.get_counts());
        
        // Making room for the previous parameters. This will be used to check if
        // the normalizing constant has been updated or not.
        params_last.push_back(stats_target[0u]);
        normalizing_constants.push_back(0.0);
        first_calc_done.push_back(false);
        
        i = arrays2support.size() - 1u;
    } else
        // Retrieving the corresponding position in the support
        i = locator->second;

    // Getting the index
    unsigned int a = arrays2support[i];
    
    // Generating a random
    std::uniform_real_distribution<> urand(0, 1);
    double r = urand(*rengine);
    double cumprob = 0.0;

    size_t k = params.size();

    // Sampling an array
    unsigned int j = 0u;
    std::vector< double > & probs = pset_probs[a];
    if ((probs.size() > 0u) && (vec_equal_approx(params, params_last[a])))
    // If precomputed, then no need to recalc support
    {

        while (cumprob < r)
            cumprob += probs[j++];

        j--;

    } else { 
       
        probs.resize(pset_arrays[a].size());
        std::vector< double > temp_stats(params.size());
        const std::vector< double > & stats = pset_stats[a];

        int i_matches = -1;
        for (size_t array = 0u; array < probs.size(); ++array)
        {

            // Filling out the parameters
            for (auto p = 0u; p < params.size(); ++p)
                temp_stats[p] = stats[array * k + p];

            probs[array] = this->likelihood(params, temp_stats, i, false);
            cumprob += probs[array];

            if (i_matches == -1 && cumprob >= r)
                i_matches = array;
        }

        j = i_matches;
        
    }
    

    return this->pset_arrays[a][j];   

}

MODEL_TEMPLATE(double, conditional_prob)(
    const Array_Type & Array_,
    const std::vector< double > & params,
    unsigned int i,
    unsigned int j
) {

    // Generating a copy of the array so we can update
    Array_Type A(Array_, true);

    // Making sure we add it first
    A.insert_cell(i, j, A.default_val(), true, false);

    // Computing the change stats_target
    std::vector< double > tmp_counts(counters->size());
    for (unsigned int ii = 0u; ii < tmp_counts.size(); ++ii)
        tmp_counts[ii] = counters->operator[](ii).count(A, i, j);

    // If there is a transformation function, it needs to be
    // applied before dealing with the likelihood.
    if (transform_model_fun)
        tmp_counts = transform_model_fun(&tmp_counts[0u], tmp_counts.size());

    return 1.0/
        (1.0 + std::exp(-vec_inner_prod<double>(
            &params[0u], &tmp_counts[0u], params.size()
            )));

    
}

MODEL_TEMPLATE(const std::mt19937 *, get_rengine)() const {
    return this->rengine;
}

template MODEL_TEMPLATE_ARGS()
inline Counters<Array_Type,Data_Counter_Type> * MODEL_TYPE()::get_counters() {
    return this->counters;
}

template MODEL_TEMPLATE_ARGS()
inline Rules<Array_Type,Data_Rule_Type> * MODEL_TYPE()::get_rules() {
    return this->rules;
}

template MODEL_TEMPLATE_ARGS()
inline Rules<Array_Type,Data_Rule_Dyn_Type> * MODEL_TYPE()::get_rules_dyn() {
    return this->rules_dyn;
}

template MODEL_TEMPLATE_ARGS()
inline Support<Array_Type,Data_Counter_Type,Data_Rule_Type,Data_Rule_Dyn_Type> *
MODEL_TYPE()::get_support_fun() {
    return &this->support_fun;
}

MODEL_TEMPLATE(std::vector< std::vector< double > > *, get_stats_target)()
{
    return &stats_target;
}

MODEL_TEMPLATE(std::vector< std::vector< double > > *, get_stats_support)()
{
    return &stats_support;
}

MODEL_TEMPLATE(std::vector< unsigned int > *, get_arrays2support)()
{
    return &arrays2support;
}

MODEL_TEMPLATE(std::vector< std::vector< Array_Type > > *, get_pset_arrays)() {
    return &pset_arrays;
}

MODEL_TEMPLATE(std::vector< std::vector<double> > *, get_pset_stats)() {
    return &pset_stats;
}

MODEL_TEMPLATE(std::vector< std::vector<double> > *, get_pset_probs)() {
    return &pset_probs;
}

MODEL_TEMPLATE(void, set_transform_model)(
    std::function<std::vector<double>(double *,unsigned int)> fun,
    std::vector< std::string > names
    )
{

    if (transform_model_fun)
        throw std::logic_error("A transformation function for the model has already been established.");
    
    transform_model_fun = fun;
    transform_model_term_names = names;

    size_t k = counters->size(); 

    // Applying over the support
    for (auto & s : stats_support)
    {

        // Making room for the new support
        std::vector< double > s_new(0u);
        s_new.reserve(s.size());

        size_t n = s.size() / (k + 1u);

        // Iterating through the unique sets
        for (size_t i = 0; i < n; ++i)
        {

            // Appending size
            s_new.push_back(s[i * (k + 1u)]);

            // Applying transformation and adding to the new set
            auto res = transform_model_fun(&s[i * (k + 1u) + 1u], k);

            if (res.size() != transform_model_term_names.size())
                throw std::length_error("The transform vector from -transform_model_fun- does not match the size of -transform_model_term_names-.");

            std::copy(res.begin(), res.end(), std::back_inserter(s_new));

        }

        // Exchanging with the original
        std::swap(s, s_new);

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
            std::vector< double > new_stats(0u);

            for (auto a = 0u; a < pset_arrays[s].size(); ++a)
            {
                // Computing the transformed version of the data
                auto tmpstats = transform_model_fun(
                    &pset_stats[s][a * k], k
                    );

                // Storing the new values
                for (auto p = 0u; p < k; ++p)
                    new_stats.push_back(tmpstats[p]);
            }

            // Updating the dataset
            std::swap(pset_stats[s], new_stats);

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
bool rule_fun_default(const Array_Type * array, uint i, uint j, Data_Type * dat) {
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
        Data_Type dat_
        ) : fun(fun_), dat(dat_) {};
    ///@}

    ~Rule() {};

    Data_Type & D(); ///< Read/Write access to the data.
    
    bool operator()(const Array_Type & a, uint i, uint j);
    
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

    uint size() const noexcept {
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
        Data_Type                           data_
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

    bool operator()(const Array_Type & a, uint i, uint j);
    
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
inline bool Rule<Array_Type,Data_Type>::operator()(const Array_Type & a, uint i, uint j) {
    return fun(a, i, j, dat);
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
        Data_Type                           data_
) {
       
    data.push_back(Rule<Array_Type,Data_Type>(
        rule_,
        data_
    ));
    
    return;
    
}

template <typename Array_Type, typename Data_Type>
inline bool Rules<Array_Type,Data_Type>::operator()(
    const Array_Type & a, uint i, uint j
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

    
    uint N = a.nrow();
    uint K = a.ncol();
    
    // Reserving some space
    (void) free->empty();
    (void) free->reserve(2u * N * K);
    
    for (uint i = 0u; i < N; ++i)
    {

        for (uint j = 0u; j < K; ++j)
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
  * @brief Data class used to store arbitrary uint or double vectors */
class NetCounterData {
public:
    
    std::vector< uint > indices;
    std::vector< double > numbers;
    
    NetCounterData() : indices(0u), numbers(0u) {};
    NetCounterData(
        const std::vector< uint > indices_,
        const std::vector< double > numbers_
    ): indices(indices_), numbers(numbers_) {};
    
    ~NetCounterData() {};
    
    // const uint get_uint
    
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
inline double (a) (const Tnet & Array, uint i, uint j, NetCounterData & data)

/**Lambda function for definition of a network counter function*/
#define NETWORK_COUNTER_LAMBDA(a) \
Counter_fun_type<Tnet, NetCounterData> a = \
    [](const Tnet & Array, uint i, uint j, NetCounterData & data)

#define NETWORKDENSE_COUNTER_LAMBDA(a) \
Counter_fun_type<NetworkDense, NetCounterData> a = \
    [](const NetworkDense & Array, uint i, uint j, NetCounterData & data)
///@}


/**@name Macros for defining rules
  */
///@{
/**Function for definition of a network counter function*/
#define NETWORK_RULE(a) \
template<typename Tnet = Network>\
inline bool (a) (const Tnet & Array, uint i, uint j, bool & data)

/**Lambda function for definition of a network counter function*/
#define NETWORK_RULE_LAMBDA(a) \
Rule_fun_type<Tnet, bool> a = \
[](const Tnet & Array, uint i, uint j, bool & data)
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
        // for (unsigned int k = 0u; k < Array.nrow(); ++k)
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
        // for (unsigned int k = 0u; k < Array.ncol(); ++k)
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
        unsigned int N = Array.nrow();

        // Self ties do not count
        if (i == j)
            return 0.0;

        // This is the first i sends, so nothing will change
        if (Array.rowsum(i) == BARRY_ZERO_NETWORK_DENSE)
            return 0.0;

        
        double ans = 0.0;
        for (unsigned int k = 0u; k < N; ++k)
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
        #ifdef __OPENM 
        #pragma omp simd reduction(+:ans)
        #else
        #pragma GCC ivdep
        #endif
        for (unsigned int k = 0u; k < Array.nrow(); ++k)
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
        for (unsigned int k = 0u; k < Array.nrow(); ++k)
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
        for (unsigned int k = 0u; k < Array.ncol(); ++k)
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
    uint attr_id,
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
    uint attr_id,
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
    uint attr_id
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
    uint attr_id
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
    uint attr_id
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
    uint attr_id
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
    std::vector< uint > d
) {

    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        uint d = Array.col(j).size();
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
    std::vector< uint > d
) {

    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {
        
        unsigned int indeg = 0u;
        for (unsigned int k = 0u; k < Array.nrow(); ++k)
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
    std::vector<uint> d
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count)
    {
        
        uint d = Array.row(i).size();
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
    std::vector<uint> d
) {
    
    NETWORKDENSE_COUNTER_LAMBDA(tmp_count)
    {
        
        uint d = 0;
        for (unsigned int k = 0u; k < Array.ncol(); ++k)
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
    std::vector<uint> d
) {
    
    NETWORK_COUNTER_LAMBDA(tmp_count) {
        
        uint d = Array.row(i).size();
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
#define CSS_SIZE() \
    uint n = data.indices[0u]; \
    uint s = data.indices[1u]; \
    uint e = data.indices[2u];

// Variables in case that the current cell corresponds to the True
#define CSS_CASE_TRUTH() if ((i < n) && (j < n)) 
#define CSS_TRUE_CELLS() \
    double tji = static_cast<double>(Array(j, i, false)); \
    double pij = static_cast<double>(Array(i + s, j + s, false)); \
    double pji = static_cast<double>(Array(j + s, i + s, false));

// Variables in case that the current cell corresponds to the Perceived
#define CSS_CASE_PERCEIVED() else if (((i >= s) && (i < e)) & ((j >= s) && (j < e)))
#define CSS_PERCEIVED_CELLS() \
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

#define CSS_CHECK_SIZE() for (uint i = 0u; i < end_.size(); ++i) {\
    if (i == 0u) continue; \
    else if (end_[i] < end_[i-1u]) \
        throw std::logic_error("Endpoints should be specified in order.");}

#define CSS_APPEND(name) std::string name_ = (name);\
    for (uint i = 0u; i < end_.size(); ++i) { \
    std::string tmpname = name_ + " (" + std::to_string(i) + ")";\
    counters->add_counter(tmp_count, tmp_init, nullptr, \
            NetCounterData({netsize, i == 0u ? netsize : end_[i-1], end_[i]}, {}),\
            tmpname);}

#define CSS_NET_COUNTER_LAMBDA_INIT() NETWORK_COUNTER_LAMBDA(tmp_init) {\
        CSS_CHECK_SIZE_INIT() \
        return 0.0; \
    };


/** @brief Counts errors of commission 
 * @param netsize Size of the reference (true) network 
 * @param end_ Vector indicating one past the ending index of each network. (see details)
 * @details 
 * The `end_` parameter should be of length `N of networks` - 1. It is
 * assumed that the first network ends at `netsize`.
 */
template<typename Tnet = Network>
inline void counter_css_partially_false_recip_commi(
    NetCounters<Tnet> * counters,
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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

        // At the beginning is all zero
        return n_dbl * (n_dbl - 1.0)/2.0;

    };
    
    // checking sizes
    CSS_CHECK_SIZE()
    CSS_APPEND("(01) Accurate null")
        
    return;

}

template<typename Tnet = Network>
inline void counter_css_census02(
    NetCounters<Tnet> * counters,
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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
    uint netsize,
    const std::vector< uint > & end_
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

#undef CSS_CASE_TRUTH
#undef CSS_TRUE_CELLS
#undef CSS_CASE_PERCEIVED
#undef CSS_PERCEIVED_CELLS
#undef CSS_CASE_ELSE
#undef CSS_CHECK_SIZE_INIT
#undef CSS_CHECK_SIZE
#undef CSS_NET_COUNTER_LAMBDA_INIT
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
    
    rules->add_rule(no_self_tie, false);
    
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
        namespace phylo {
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/counters/phylo.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRAY_PHYLO_H
#define BARRAY_PHYLO_H 1

// Default value that is used for the counters.
#define DEFAULT_DUPLICATION 1u
#define DUPL_SPEC 0u
#define DUPL_DUPL 1u
#define DUPL_EITH 2u


#define MAKE_DUPL_VARS() \
    bool DPL = Array.D_ptr()->duplication; \
    unsigned int DATA_AT = data[0u];

#define IS_EITHER()      (DATA_AT == DUPL_EITH)
#define IS_DUPLICATION() ((DATA_AT == DUPL_DUPL) & (DPL))
#define IS_SPECIATION()  ((DATA_AT == DUPL_SPEC) & (!DPL))

#define IF_MATCHES() MAKE_DUPL_VARS() \
    if (IS_EITHER() | IS_DUPLICATION() | IS_SPECIATION())
#define IF_NOTMATCHES() MAKE_DUPL_VARS() \
    if (!IS_EITHER() & !IS_DUPLICATION() & !IS_SPECIATION())


/**
 * @ingroup counting
 * @details Details about the available counters for `PhyloArray`
 * objects can be found in the \ref counters-phylo section.
 */
///@{

/**
 * @brief Data definition for the `PhyloArray` class.
 * 
 * This holds basic information about a given node.
 * 
 */
class NodeData {
public:
  
    /**
     * Branch length.
     */
    std::vector< double > blengths = {};
    
    /**
     * State of the parent node.
     */
    std::vector< bool > states = {};
    
    /**
     * 
     */
    bool duplication = true;
    
    // NodeData() : blengths(0u), states(0u) {};
    
    NodeData(
        const std::vector< double > & blengths_,
        const std::vector< bool > & states_,
        bool duplication_ = true
    ) : blengths(blengths_), states(states_), duplication(duplication_) {};
    
    // ~NodeData() {};
  
};

// typedef std::vector< uint > PhyloCounterData;
class PhyloCounterData {
private:
    std::vector< uint > data;
    std::vector< double > * counters;

public:
    PhyloCounterData(
        std::vector< uint > data_,
        std::vector< double > * counters_ = nullptr
        ) : data(data_), counters(counters_) {};

    PhyloCounterData() : data(0u) {};

    uint at(uint d) {return data.at(d);};
    uint operator()(uint d) {return data.at(d);};
    uint operator[](uint d) {return data[d];};
    void reserve(uint x) {return data.reserve(x);};
    void push_back(uint x) {return data.push_back(x);};
    void shrink_to_fit()  {return data.shrink_to_fit();};
    uint size() {return data.size();};

    std::vector< uint >::iterator begin() {return data.begin();};
    std::vector< uint >::iterator end() {return data.end();};

    bool empty() {return data.empty();};
    std::vector< double > * get_counters() {return counters;};

};


typedef std::vector< std::pair< uint, uint > > PhyloRuleData;
class PhyloRuleDynData;

/**
 * @name Convenient typedefs for Node objects.
 * */
///@{
typedef BArrayDense<uint, NodeData> PhyloArray;
typedef Counter<PhyloArray, PhyloCounterData > PhyloCounter;
typedef Counters< PhyloArray, PhyloCounterData> PhyloCounters;

typedef Rule<PhyloArray,PhyloRuleData> PhyloRule;
typedef Rules<PhyloArray,PhyloRuleData> PhyloRules;

typedef Rule<PhyloArray,PhyloRuleDynData> PhyloRuleDyn;
typedef Rules<PhyloArray,PhyloRuleDynData> PhyloRulesDyn;

typedef Support<PhyloArray, PhyloCounterData, PhyloRuleData, PhyloRuleDynData > PhyloSupport;
typedef StatsCounter<PhyloArray, PhyloCounterData> PhyloStatsCounter;
typedef Model<PhyloArray, PhyloCounterData, PhyloRuleData, PhyloRuleDynData > PhyloModel;
typedef PowerSet<PhyloArray, PhyloRuleData> PhyloPowerSet;
///@}


/**
 * @brief Extension of a simple counter.
 * 
 * It allows specifying extra arguments, in particular, the corresponding
 * sets of rows to which this statistic may be relevant. This could be important
 * in the case of, for example, counting correlation type statistics between
 * function 1 and 2, and between function 1 and 3.
 * 
 * 
 */
#define PHYLO_COUNTER_LAMBDA(a) Counter_fun_type<PhyloArray, PhyloCounterData> a = \
    [](const PhyloArray & Array, uint i, uint j, PhyloCounterData & data)

#define PHYLO_RULE_DYN_LAMBDA(a) Rule_fun_type<PhyloArray, PhyloRuleDynData> a = \
    [](const PhyloArray & Array, uint i, uint j, PhyloRuleDynData & data)

#define PHYLO_CHECK_MISSING() if (Array.D_ptr() == nullptr) \
    throw std::logic_error("The array data is nullptr."); \
    
inline std::string get_last_name(unsigned int d) {return ((d == 1u)? " at duplication" : ((d == 0u)? " at speciation" : ""));}

/**
 * @weakgroup counters-phylo Phylo counters
 * @brief Counters for phylogenetic modeling
 * @param counters A pointer to a `PhyloCounters` object (`Counters`<`PhyloArray`, `PhyloCounterData`>).
 */
///@{
// -----------------------------------------------------------------------------
/**
 * @brief Overall functional gains
 * @details Total number of gains (irrespective of the function).
 */
inline void counter_overall_gains(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        
        return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {
        IF_NOTMATCHES()
            return 0.0;
      
        return Array.D_ptr()->states[i] ? 0.0 : 1.0;
      
    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Overall gains" + get_last_name(duplication)
    );

    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Functional gains for a specific function (`nfun`).
 */
inline void counter_gains(
    PhyloCounters * counters,
    std::vector<uint> nfun,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        IF_NOTMATCHES()
            return 0.0;

        double ngains = 0.0;
        auto   k = data[1u];
        auto   s = Array.D_ptr()->states[k];

        if (s)
            return 0.0;

        for (auto o = 0u; o < Array.ncol(); ++o)
        {
            if (!s && (Array(k,o) == 1u))
                ngains += 1.0;
        }
        
        return ngains;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is there any gain?
        if (Array.D_ptr()->states[i])
            return 0.0;

        IF_MATCHES()
            return (i == data[1u]) ? 1.0 : 0.0;
        
        return 0.0;

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            PhyloCounterData({duplication, i}),
            "Gains " + std::to_string(i) + get_last_name(duplication)
        );
    
    return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief k genes gain function nfun
 */
inline void counter_gains_k_offspring(
    PhyloCounters * counters,
    std::vector<uint> nfun,
    uint k = 1u,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

      PHYLO_CHECK_MISSING();
      return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this relevant?
        if (i != data[1u])
            return 0.0;

        IF_NOTMATCHES()
            return 0.0;

        // Is there any gain?
        if (Array.D_ptr()->states[i])
            return 0.0;

        // Making the counts
        int counts = 0;
        for (uint k = 0u; k < Array.ncol(); ++k)
            if (k != j)
            {
                if (Array(i, k, false) == 1u)
                    ++counts;
            }

        // Three cases: base on the diff
        int diff = static_cast<int>(data[2u]) - counts + 1;
        // (a) counts were 1 below k, then +1
        if (diff == 1)
            return -1.0;
            // (b) counts were equal to k, then -1
        else if (diff == 0)
        {
            return 1.0;
        } else 
            // (c) Otherwise, nothing happens
            return 0.0;
      

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            PhyloCounterData({duplication, i, k}),
            std::to_string(k) + " genes gain " + std::to_string(i) +
                get_last_name(duplication)
        );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Keeps track of how many genes are changing (either 0, 1, or 2 if dealing
 * with regular trees.)
 */
inline void counter_genes_changing(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        for (auto s : Array.D_ptr()->states)
        {

            if (s) 
                // Yup, we are loosing a function, so break
                return static_cast<double>(Array.ncol());
            
        }

        return 0.0;
      

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;

        // Need to check the other functions
        for (uint k = 0u; k < Array.nrow(); ++k)
        {

            // Nah, this gene was already different.
            if ((k != i) && (Array.D_ptr()->states[k] != (Array(k, j, false) == 1u)))
                return 0.0;
            

        }

        // Nope, this gene is now matching its parent, so we need to 
        // take it out from the count of genes that have changed.
        return Array.D_ptr()->states[i] ? -1.0 : 1.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Num. of genes changing" + get_last_name(duplication)
    );

    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Keeps track of how many pairs of genes preserve pseudostate.
 */
inline void counter_preserve_pseudogene(
    PhyloCounters * counters,
    unsigned int nfunA,
    unsigned int nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        if (Array.D_ptr()->states[data[1u]] | Array.D_ptr()->states[data[2u]])
            return 0.0;

        double n = static_cast<double>(Array.ncol());
        return n * (n - 1.0) / 2.0;
      

    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;

        auto nfunA = data[1u];
        auto nfunB = data[2u];

        if ((i != nfunA) & (i != nfunB))
            return 0.0;

        if (Array.D_ptr()->states[data[1u]] | Array.D_ptr()->states[data[2u]])
            return 0.0;

        unsigned int k = (i == nfunA) ? nfunB : nfunA;

        if (Array(k, j) == 1u)
            return 0.0;

        double res = 0.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {
            if (off == j)
                continue;

            if ((Array(i, off) == 0u) && (Array(k, off) == 0u))
                res -= 1.0;

        }

        return res;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Preserve pseudo gene (" + 
        std::to_string(nfunA) + ", " +
        std::to_string(nfunB) + ")" + get_last_name(duplication)
    );

    
    return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief Keeps track of how many genes are changing (either 0, 1, or 2 if dealing
 * with regular trees.)
 */
inline void counter_prop_genes_changing(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        for (auto s : Array.D_ptr()->states)
        {
            if (s)
                return 1.0;
        }
        
        return 0.0;
      
    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;
        
        // Setup
        bool j_diverges = false;
        const std::vector< bool > & par_state = Array.D_ptr()->states;

        for (unsigned int f = 0u; f < Array.nrow(); ++f)
        {

            // Was the gene annotation different from the parent?
            if (par_state[f] != (Array(f,j) == 1u))
            {
                j_diverges = true;
                break;
            }

        }


        bool j_used_to_diverge = false;
        for (unsigned int f = 0u; f < Array.nrow(); ++f)
        {

            if (f == i)
            {
                if (par_state[f])
                {
                    j_used_to_diverge = true;
                    break;
                }
            }
            else
            {

                if (par_state[f] != (Array(f,j) == 1u))
                {
                    j_used_to_diverge = true;
                    break;
                }

            }

        }

        // Case 1: j hasn't changed
        if ((!j_used_to_diverge & !j_diverges) | (j_used_to_diverge & j_diverges))
            return 0.0;
        // Case 2: j NOW diverges
        else if (j_diverges)
            return 1.0/Array.ncol();
        // Case 3: j USED to diverge, so no more
        else
            return -1.0/Array.ncol();

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Proportion of genes changing" + get_last_name(duplication)
    );

    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Overall functional loss
 */
inline void counter_overall_loss(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
    )
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        if (!Array.D_ptr()->states[i])
            return 0.0;

        IF_MATCHES()
            return -1.0;
        else 
            return 0.0;
        
    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        IF_NOTMATCHES()
            return 0.0;
        
        double res = 0.0;
        for (auto s : Array.D_ptr()->states)
            if (s)
                res += 1.0;

        return res * static_cast<double>(Array.ncol());

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Overall loses" + get_last_name(duplication)
    );
    
    return;

}

// -----------------------------------------------------------------------------
/**
 * @brief Cap the number of functions per gene
 */
inline void counter_maxfuns(
    PhyloCounters * counters,
    uint            lb,
    uint            ub,
    unsigned int duplication = DEFAULT_DUPLICATION
 )
 {

    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();    

        IF_NOTMATCHES()
            return 0.0;

        // At first, all are zero, so we need to check if the lower
        // bound is zero
        if (data[1u] == 0)
            return static_cast<double>(Array.ncol());
        
        return 0.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        int count = Array.colsum(j);
        int ub    = data[2u];
        
        // It now matches
        if (count == static_cast<int>(data[1u]))
            return 1.0;

        // Was within, but now outside
        if (count > ub && ((count - ub) == 1))
            return -1.0;

        // Otherwise nothing happens.
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, lb, ub}),
        "Genes with [" + std::to_string(lb) + ", " + std::to_string(ub) +
            "] funs" + get_last_name(duplication)
    );
    
    return;
  
}
  
// -----------------------------------------------------------------------------
/**
 * @brief Total count of losses for an specific function.
 */
inline void counter_loss(
    PhyloCounters * counters,
    std::vector<uint> nfun,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        if (!Array.D_ptr()->states[i])
            return 0.0;
        
        return (i == data[1u]) ? -1.0 : 0.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;
        
        auto f = data[1u];

        if (!Array.D_ptr()->states[f])
            return 0.0;
        
        return static_cast<double>(Array.ncol());

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            PhyloCounterData({duplication, i}),
            "Loss " + std::to_string(i) + get_last_name(duplication)
        );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Total number of changes. Use this statistic to account for "preservation"
 */
inline void counter_overall_changes(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        if (Array.D_ptr()->states[i])
            return -1.0;
        else 
            return 1.0;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        IF_NOTMATCHES()
            return 0.0;

        PHYLO_CHECK_MISSING();


        // Since we start with all the array at zero,
        // As many chances to change as offspring
        double noff   = static_cast<double> (Array.ncol());
        double counts = 0.0;
        for (uint k = 0u; k < Array.nrow(); ++k)
            if (Array.D_ptr()->states[k])
                counts += noff;

        return counts;

        

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Overall changes" + get_last_name(duplication)
    );
    
    
    return;
  
}


// -----------------------------------------------------------------------------
/**
 * @brief Total count of Sub-functionalization events.
 * @details It requires to specify data = {funA, funB}
 */
inline void counter_subfun(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;

        auto funA = data[1u];
        auto funB = data[2u];
        
        // Are we looking at either of the relevant functions?
        if ((funA != i) && (funB != i))
            return 0.0;
        
        // Are A and B existant? if not, no change
        if (!Array.D_ptr()->states[funA] | !Array.D_ptr()->states[funB])
            return 0.0;
        
        // Figuring out which is the first (reference) function
        uint other = (i == funA)? funB : funA;
        double res = 0.0;
        // There are 4 cases: (first x second) x (had the second function)
        if (Array(other, j, false) == 1u)
        { 
          
            for (uint off = 0u; off < Array.ncol(); ++off)
            {
                
                // Not on self
                if (off == j)
                    continue;
                
                if ((Array(i, off, false) == 1u) && (Array(other, off, false) == 0u))
                    res -= 1.0;
                
            }
          
        } else {
          
            for (uint off = 0u; off < Array.ncol(); ++off)
            {
              
                // Not on self
                if (off == j)
                    continue;
                
                if ((Array(i, off, false) == 0u) && (Array(other, off, false) == 1u))
                    res += 1.0;
              
            }
          
        }
        
        return res;

    };

    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Subfun between " + std::to_string(nfunA) + " and " +
            std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Co-evolution (joint gain or loss)
 * @details Needs to specify pairs of functions (`nfunA`, `nfunB`).
 */
inline void counter_cogain(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        auto d1 = data[1u];
        auto d2 = data[2u];
      
        // Is the function in scope relevant?
        if ((i != d1) && (i != d2))
            return 0.0;
        
        // None should have it
        if (!Array.D_ptr()->states[d1] && !Array.D_ptr()->states[d2])
        {

            uint other = (i == d1)? d2 : d1;

            if (Array(other, j, false) == 1u)
                return 1.0;
            else
                return 0.0;

        } else 
            return 0.0;

    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Co-gains " + std::to_string(nfunA) + " & " + std::to_string(nfunB) +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/** @brief Longest branch mutates (either by gain or by loss) */
inline void counter_longest(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
    )
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // Figuring out which match
        std::vector< bool> is_longest(Array.ncol(), false);
        bool j_mutates = false;
        int nmutate = 0;
        int nmutate_longest = 0;

        auto states  = Array.D_ptr()->states;
        
        for (auto off = 0u; off < Array.ncol(); ++off)
        {

            // On the fly, figuring out if it is longest
            for (auto & l : data)
                if (l == off)
                    is_longest[off] = true;

            for (auto f = 0u; f < Array.nrow(); ++f)
            {
                if ((Array(f, off) == 1u) != states[f])
                {
                    
                    // If it happens that j != off and is not longest
                    // then return 0 (a not longest was mutating prev)
                    if (is_longest[off] && (off != j))
                        return 0.0;

                    if (off == j)
                        j_mutates = true;

                    if (is_longest[j])
                        nmutate_longest++;
                    else
                        nmutate++;

                    break;
                }

            }
        }

        // There was already more than one in difference
        // so nothing to change
        if (std::fabs(nmutate - nmutate_longest) > 1)
            return 0.0;

        // Figuring out previously
        bool j_mutates_prev = false;
        for (auto f = 0u; f < Array.nrow(); ++f)
        {
            // Checking the previous function... was it
            // different before?
            if ((f == i) && states[i])
            {
                j_mutates_prev = true;
                break;
            }
            else if ((Array(f, j) == 1u) != states[f])
            {
                j_mutates_prev = true;
                break;
            }

        }

        // Adjusting the previous count
        auto nmutate_prev         = nmutate;
        auto nmutate_longest_prev = nmutate_longest;
        if (j_mutates & !j_mutates_prev)
        {
            if (is_longest[j])
                nmutate_longest_prev--;
            else
                nmutate_prev--;
        }
        else if (!j_mutates & j_mutates)
        {
            if (is_longest[j])
                nmutate_longest_prev++;
            else
                nmutate_prev++;

        }
        
        // Just compute the change statistic directly
        return
            ( ((nmutate == 0) & (nmutate_longest > 0)) ? 1.0 : 0.0 ) +
            ( ((nmutate_prev == 0) & (nmutate_longest_prev > 0)) ? 1.0 : 0.0 );

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        
        if (Array.D_ptr()->blengths.size() != Array.ncol())
            throw std::logic_error(
                "longest should be initialized with a vec of size Array.ncol()."
            );
          
        // Finding the longest branch (or branches) --
        uint longest_idx = 0u;
        double diff      = 0.0;
        data.reserve(Array.ncol()); 
        data.push_back(0u);
        for (uint ii = 1u; ii < Array.ncol(); ++ii)
        {
            
            diff = Array.D_ptr()->blengths[longest_idx] - Array.D_ptr()->blengths[ii];
            if (diff > 0.0)
                continue;
            else if (diff < 0.0)
            {

                data.empty();
                data.push_back(ii);
                longest_idx = ii;

            }
            else if (diff == 0.0)
                data.push_back(ii);
            
        }

        data.shrink_to_fit();
        
        if (data.size() == 0u)
            throw std::logic_error("The data on the longest branch has size 0.");
        
        // Starting the counter, since all in zero, then this will be equal to
        // the number of functions in 1 x number of longest branches
        for (uint ii = 0u; ii < Array.nrow(); ++ii)
        {
            
            if (Array.D_ptr()->states[ii])
                return (1.0 * static_cast<double>(data.size()));

        }
        
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Longest branch mutates" + get_last_name(duplication)
    );
    
    return;
  
}

//------------------------------------------------------------------------------
/**
 * @brief Total number of neofunctionalization events 
 * @details Needs to specify pairs of function.
 */
inline void counter_neofun(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;
        
        auto funA = data[1u];
        auto funB = data[2u];

        // Is the function in scope relevant?
        if ((i != funA) && (i != funB))
            return 0.0;
        
        // Checking if the parent has both functions
        uint other = (i == funA)? funB : funA;
        bool parent_i     = Array.D_ptr()->states[i];
        bool parent_other = Array.D_ptr()->states[other];
        
        if (!parent_i & !parent_other) 
            return 0.0;
        else if (parent_i & parent_other) 
            return 0.0;
        
        // Figuring out which is the first (reference) function
        double res = 0.0;
        
        if (Array(other, j) == 0u)
        {


            for (auto off = 0u; off < Array.ncol(); ++off)
                if ((Array(i,off) == 0) && (Array(other,off) == 1))
                    res += 1.0;

        }
        else
        {

            for (auto off = 0u; off < Array.ncol(); ++off)
                if ((Array(i,off) == 1) && (Array(other,off) == 0))
                    res -= 1.0;
                
        }
             
        return res;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Neofun between " + std::to_string(nfunA) + " and " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
  
}

//------------------------------------------------------------------------------
/**
 * @brief Total number of neofunctionalization events 
 * sum_u sum_{w < u} [x(u,a)*(1 - x(w,a)) + (1 - x(u,a)) * x(w,a)]
 * change stat: delta{x(u,a): 0->1} = 1 - 2 * x(w,a)
 */
inline void counter_pairwise_neofun_singlefun(
    PhyloCounters * counters,
    uint nfunA,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;
        
        // Is the function in scope relevant?
        if (i != data[1u])
            return 0.0;
        
        // Checking if the parent has the function
        if (Array.D_ptr()->states[i])
            return 0.0;
        
        // Figuring out which is the first (reference) function
        double res = 0.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {

            if (off == j)
                continue;

            if ((Array(i, off) == 0))
                res += 1.0;
            else 
                res -= 1.0;

        }

        return res;

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA}),
        "Pairwise neofun function " + std::to_string(nfunA) +
        get_last_name(duplication)
    );
    
    return;
  
}

//------------------------------------------------------------------------------
/**
 * @brief Total number of neofunctionalization events 
 * @details Needs to specify pairs of function.
 */
inline void counter_neofun_a2b(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Is this node duplication?
        IF_NOTMATCHES()
            return 0.0;
        
        const uint & funA = data[1u];
        const uint & funB = data[2u];

        // Checking scope
        if ((i != funA) && (i != funB))
            return 0.0;
        
        // Checking the parent doesn't have funA or has funB
        if (!Array.D_ptr()->states[funA] | Array.D_ptr()->states[funB]) 
            return 0.0;

        double res = 0.0;

        if (i == funA)
        {

            if (Array(funB, j) == 0u)
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {

                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 0u) && (Array(funB, off) == 1u))
                        res += 1.0;

                }

            }
            else
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {
                    
                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 1u) && (Array(funB, off) == 0u))
                        res -= 1.0;

                }

            }

        }
        else
        {

            if (Array(funA, j) == 0u)
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {

                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 1u) && (Array(funB, off) == 0u))
                        res += 1.0;

                }

            }
            else
            {

                for (auto off = 0u; off < Array.ncol(); ++off)
                {
                    
                    if (off == j)
                        continue;

                    if ((Array(funA, off) == 0u) && (Array(funB, off) == 1u))
                        res -= 1.0;

                }

            }

        }

        return res;
        
    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Neofun from " + std::to_string(nfunA) + " to " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    return;
    
}

// -----------------------------------------------------------------------------
/**
 * @brief Function co-opting
 * @details Function co-opting of functions A and B happens when, for example,
 * function B is gained as a new featured leveraging what function A already does;
 * without losing function A. The sufficient statistic is defined as follows:
 * \f[
 * x_{pa}(1 - x_{pb})\sum_{i<j}\left[x_{ia}^p(1 - x_{ib}^p)x_{ja}^px_{jb}^p + x_{ja}^p(1 - x_{jb}^p)x_{ia}^px_{ib}^p\right]
 * \f]
 * This algorithm implements the change statistic.
 */
inline void counter_co_opt(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB, 
    unsigned int duplication = DEFAULT_DUPLICATION
) {
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    { 

        // Checking whether this is for duplication or not
        IF_NOTMATCHES()
            return 0.0;
        
        const unsigned int funA = data[1u];
        const unsigned int funB = data[2u];

        // If the change is out of scope, then nothing to do
        if ((i != funA) & (i != funB))
            return 0.0;

        // If the parent does not have the initial state, then it makes no sense
        if ((!Array.D_ptr()->states[funA]) | Array.D_ptr()->states[funB])
            return 0.0;

        // Checking whether function A or function B changed
        if (i == funA) {

            // What was the state of the other function? If B is present, then
            // nothing changes.
            if (Array(funB, j, false) == 1u) 
                return 0.0;

            // Iterating through the sibs
            double res = 0.0;
            for (auto c = 0u; c < Array.ncol(); ++c)
                if ((c != j) && (Array(funA, c, false) == 1u) && (Array(funB, c, false) == 1u))
                    res += 1.0;

            return res;

        } else {

            // What was the state of the other function? If A is not present, then
            // nothing changes.
            if (Array(funA, j, false) == 0u) 
                return 0.0;

            // Iterating through the sibs
            double res = 0.0;
            for (auto c = 0u; c < Array.ncol(); ++c)
                if ((c != j) && (Array(funA, c, false) == 1u))
                    res += (Array(funB, c, false) == 0u) ? 1.0 : -1.0;

            return res;

        }

        

    };
    
    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        if (data.size() != 3u)
            throw std::length_error("The counter data should be of length 2.");

        if (data[1u] == data[2u])
            throw std::logic_error("Functions A and B should be different from each other.");

        if (data[1u] >= Array.nrow())
            throw std::length_error("Function A in counter out of range.");

        if (data[2u] >= Array.nrow())
            throw std::length_error("Function B in counter out of range.");

        return 0.0;

    };

    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Coopt of " + std::to_string(nfunA) + " by " +
        std::to_string(nfunB) + get_last_name(duplication)
    );
    
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Indicator function. Equals to one if \f$k\f$ genes changed and zero
 * otherwise.
 */
inline void counter_k_genes_changing(
    PhyloCounters * counters,
    unsigned int k,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        // At the beginning, all offspring are zero, so we need to
        // find at least one state = true.
        for (auto s : Array.D_ptr()->states)
            if (s)
                return Array.ncol() == data[1u] ? 1.0 : 0.0;

        return data[1u] == 0 ? 1.0 : 0.0;
      
    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;
        
        // How many genes diverge the parent
        int              count = 0; 
        bool        j_diverges = false;
        const auto & par_state = Array.D_ptr()->states;

        int k = static_cast<int>(data[1u]);

        for (auto o = 0u; o < Array.ncol(); ++o)
        {

            for (auto f = 0u; f < Array.nrow(); ++f)
            {

                // Was the gene annotation different from the parent?
                if ((Array(f, o) == 1u) != par_state[f])
                {

                    if (o == j)
                        j_diverges = true;

                    count++;
                    break;

                }

            }

        }

        // Counts will only be relevant if (count - k) > 1. Otherwise,
        // having the j gene changed is not relevant
        if (std::abs(count - k) > 1)
            return 0.0;

        // Did it used to diverge?
        bool j_used_to_diverge = false;
        for (auto f = 0u; f < Array.nrow(); ++f)
        {

            if (f == i)
            {
                if (par_state[f]) // Since it is now true, it used to diverge
                {
                    j_used_to_diverge = true;
                    break;
                }
            }
            else
            {

                if (par_state[f] != (Array(f,j) == 1u))
                {
                    j_used_to_diverge = true;
                    break;
                }

            }

        }

        auto count_prev = count;
        // Case 1: j hasn't changed
        if ((!j_used_to_diverge & !j_diverges) | (j_used_to_diverge & j_diverges))
            return 0.0;
        // Case 2: j NOW diverges
        else if (j_diverges)
            count_prev--;
        // Case 3: j USED to diverge
        else
            count_prev++;

        return (count == k ? 1.0 : 0.0) - (count_prev == k ? 1.0 : 0.0);

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, k}),
        std::to_string(k) + " genes changing" + get_last_name(duplication)
    );
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Indicator function. Equals to one if \f$k\f$ genes changed and zero
 * otherwise.
 */
inline void counter_less_than_p_prop_genes_changing(
    PhyloCounters * counters,
    double p,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_init)
    {
        
        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        for (auto s : Array.D_ptr()->states)
            if (s)
                return data[1u] == 100 ? 1.0 : 0.0;

        // Only one if it was specified it was zero
        return 1.0;
      
    };

    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        // Checking the type of event
        IF_NOTMATCHES()
            return 0.0;
        
        // Setup
        double count = 0.0; ///< How many genes diverge the parent

        bool j_diverges = false;
        const std::vector< bool > & par_state = Array.D_ptr()->states;

        for (unsigned int o = 0u; o < Array.ncol(); ++o)
        {

            for (unsigned int f = 0u; f < Array.nrow(); ++f)
            {

                // Was the gene annotation different from the parent?
                if ((Array(f, o) == 1u) != par_state[f])
                {

                    if (o == j)
                        j_diverges = true;

                    count += 1.0;
                    break;

                }

            }

        }


        bool j_used_to_diverge = false;
        for (unsigned int f = 0u; f < Array.nrow(); ++f)
        {

            if (f == i)
            {
                if (par_state[f])
                {
                    j_used_to_diverge = true;
                    break;
                }
            }
            else
            {

                if (par_state[f] != (Array(f,j) == 1u))
                {
                    j_used_to_diverge = true;
                    break;
                }

            }

        }

        auto count_prev = count;
        // Case 1: j hasn't changed
        if ((!j_used_to_diverge & !j_diverges) | (j_used_to_diverge & j_diverges))
            return 0.0;
        // Case 2: j NOW diverges
        else if (j_diverges)
            count_prev -= 1.0;
        // Case 3: j USED to diverge
        else
            count_prev += 1.0;

        double ncol = static_cast<double>(Array.ncol());
        double p    = static_cast<double>(data[1u]) / 100.0;

        return ((count/ncol) <= p ? 1.0 : 0.0) - ((count_prev/ncol) <= p ? 1.0 : 0.0);

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, static_cast<uint>(p * 100)}),
        std::to_string(p) + " prop genes changing" + get_last_name(duplication)
    );
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 */
inline void counter_gains_from_0(
    PhyloCounters * counters,
    std::vector< uint > nfun,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // All must be false
        for (auto s : Array.D_ptr()->states)
        {

            if (s)
                return 0.0;

        }

        // Is this the function?
        if (i != data[1u])
            return 0.0;

        // Now computing the change stats
        double res = static_cast<double>(Array.ncol()) - 1.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {
            if (off  == j)
                continue;

            if (Array(i, off) == 1u)
                res -= 2.0;
        }


        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    for (auto& i : nfun)
        counters->add_counter(
            tmp_count, tmp_init, nullptr,
            PhyloCounterData({duplication, i}),
            "First gain " + std::to_string(i) +
                get_last_name(duplication)
        );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 */
inline void counter_overall_gains_from_0(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // All must be false
        for (auto s : Array.D_ptr()->states)
        {

            if (s)
                return 0.0;

        }

        return 1.0;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Overall first gains" +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 */
inline void counter_pairwise_overall_change(
    PhyloCounters * counters,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        unsigned int funpar = Array.D_ptr()->states[i] == 1u;

        // All must be false
        double res = 0.0;
        for (auto off = 0u; off < Array.ncol(); ++off)
        {
            if (off == j)
                continue;

            if (funpar > Array(i, off))
                res -= 1.0;
            else if (funpar < Array(i, off))
                res += 1.0;
        }
        
        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();

        IF_NOTMATCHES()
            return 0.0;

        double res = 0.0;
        double n   = static_cast<double>(Array.ncol());
        for (auto s : Array.D_ptr()->states)
            if (s)
                res += n * (n - 1.0) / 2.0;

        return res;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication}),
        "Pairs of genes changing" +
            get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 * sum x(a)^3(1-x(b))^3 + x(b)^3(1-x(a))^3 + x(a)^3 * x(b)^3 + (1 - x(a))^3 * (1-x(b))^3
 */
inline void counter_pairwise_preserving(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // Not in the scope
        auto funA = data[1u];
        auto funB = data[2u];
        if ((funA != i) && (funB != i))
            return 0.0;

        unsigned int k = (funA == i) ? funB : funA;

        bool parent_i = Array.D_ptr()->states[i];
        bool parent_k = Array.D_ptr()->states[k];

        // if (!parent_i & !parent_k)
        //     return 0.0;

        double res = 0.0;
        // Case 1: (0,0)
        if (!parent_i & !parent_k)
        {

            if (Array(k, j) == 1u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 0u) && (Array(k, off) == 0u))
                    res -= 1.0;

            }

        }
        else if (parent_i & !parent_k)
        {

            if (Array(k, j) == 1u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 1u) && (Array(k, off) == 0u))
                    res += 1.0;

            }

        }
        else if (!parent_i & parent_k)
        {

            if (Array(k, j) == 0u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 0u) && (Array(k, off) == 1u))
                    res += 1.0;

            }

        }
        else
        {

            if (Array(k, j) == 0u)
                return 0.0; 

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 1u) && (Array(k, off) == 1u))
                    res += 1.0;

            }
        }

        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {


        IF_NOTMATCHES()
            return 0.0;
        
        PHYLO_CHECK_MISSING();
        
        double n = static_cast< double >(Array.ncol());
        if (!Array.D_ptr()->states[data[1u]] && !Array.D_ptr()->states[data[2u]])
            return n * (n - 1.0) / 2.0;

        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "Pariwise preserve (" + std::to_string(nfunA) + ", " +
            std::to_string(nfunB) + ")" +get_last_name(duplication)
    );
    
    return;
  
}

// -----------------------------------------------------------------------------
/**
 * @brief Used when all the functions are in 0 (like the root node prob.)
 * @details Needs to specify function a.
 * sum x(a)^3(1-x(b))^3 + x(b)^3(1-x(a))^3 + x(a)^3 * x(b)^3 + (1 - x(a))^3 * (1-x(b))^3
 */
inline void counter_pairwise_first_gain(
    PhyloCounters * counters,
    uint nfunA,
    uint nfunB,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_COUNTER_LAMBDA(tmp_count)
    {

        IF_NOTMATCHES()
            return 0.0;

        // Not in the scope
        auto funA = data[1u];
        auto funB = data[2u];
        if ((funA != i) && (funB != i))
            return 0.0;

        unsigned int k = (funA == i) ? funB : funA;

        double res = 0.0;
        if (Array(k, j) == 1)
        {

            for (auto off = 0u; off < Array.ncol(); ++off)
            {
                if (off == j)
                    continue;

                if ((Array(i,off) == 0u) && (Array(k,off) == 0u))
                    res -= 1.0;
            }

        }
        else
        {

            for (auto off = 0u; off < Array.ncol(); ++off)
            {

                if (off == j)
                    continue;

                if ((Array(i, off) == 1u))
                {

                    // j: (0,0)\(1,0) -> (1,0)\(1,0), so less 1
                    if (Array(k, off) == 0u)
                        res -= 1.0;

                }
                else
                {

                    if (Array(k, off) == 1u) 
                    // j: (0,0)\(0,1) -> (1,0)\(0,1), so less 1
                        res -= 1.0;
                    else
                    // j: (0,0)\(0,0) -> (1,0)\(0,0), so plus 1
                        res += 1.0;

                }

            }

        }
        

        return res;
        
    };

    PHYLO_COUNTER_LAMBDA(tmp_init) {

        PHYLO_CHECK_MISSING();
        
        return 0.0;

    };
    
    counters->add_counter(
        tmp_count, tmp_init, nullptr,
        PhyloCounterData({duplication, nfunA, nfunB}),
        "First gain (either " + std::to_string(nfunA) + " or " +
            std::to_string(nfunB) + ")" +get_last_name(duplication)
    );
    
    return;
  
}

///@}

/**
 * @weakgroup rules-phylo Phylo rules
 * @brief Rules for phylogenetic modeling
 * @param rules A pointer to a `PhyloRules` object (`Rules`<`PhyloArray`, `PhyloRuleData`>).
 */
///@{

class PhyloRuleDynData {
public:
    const std::vector< double > * counts;
    uint pos;
    uint lb;
    uint ub;
    uint duplication;

    PhyloRuleDynData(
        const std::vector< double > * counts_,
        uint pos_,
        uint lb_,
        uint ub_,
        uint duplication_
        ) :
        counts(counts_), pos(pos_), lb(lb_), ub(ub_), duplication(duplication_) {};
    
    ~PhyloRuleDynData() {};
    
};

/**
 * @brief Overall functional gains
 * @param support Support of a model.
 * @param pos Position of the focal statistic.
 * @param lb Lower bound
 * @param ub Upper bound
 * @details 
 * @return (void) adds a rule limiting the support of the model.
 */
inline void rule_dyn_limit_changes(
    PhyloSupport * support,
    uint pos,
    uint lb,
    uint ub,
    unsigned int duplication = DEFAULT_DUPLICATION
)
{
  
    PHYLO_RULE_DYN_LAMBDA(tmp_rule)
    {

        unsigned int rule_type = data.duplication;
        if (rule_type != DUPL_EITH)
        {

            if (Array.D_ptr()->duplication & (rule_type != DUPL_DUPL))
                return true;
            else if (!Array.D_ptr()->duplication & (rule_type != DUPL_SPEC))
                return true;
                
        }

        if (data.counts->operator[](data.pos) < data.lb)
            return false;
        else if (data.counts->operator[](data.pos) > data.ub)
            return false;
        else
            return true;
      
    };
    
    support->get_rules_dyn()->add_rule(
        tmp_rule,
        PhyloRuleDynData(
            support->get_current_stats(),
            pos, lb, ub, duplication
            )
    );
    
    return;
  
}

///@}

#undef MAKE_DUPL_VARS
#undef IS_EITHER
#undef IS_DUPLICATION
#undef IS_SPECIATION
#undef IF_MATCHES
#undef IF_NOTMATCHES

#undef DEFAULT_DUPLICATION
#undef DUPL_SPEC
#undef DUPL_DUPL
#undef DUPL_EITH


#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/counters/phylo.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


        }
        namespace defm {
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry/counters/defm.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRAY_DEFM_H
#define BARRAY_DEFM_H 1

/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 Start of -include/barry//counters/defm-formula.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


#ifndef BARRY_DEFM_MOTIF_FORMULA_HPP
#define BARRY_DEFM_MOTIF_FORMULA_HPP
/**
 * @brief Parses a motif formula
 * 
 * @details This function will take the formula and generate the corresponding
 * input for defm::counter_transition(). Formulas can be specified in the
 * following ways:
 * 
 * - Intercept effect: {...} No transition, only including the current state.
 * - Transition effect: {...} > {...} Includes current and previous states.
 * 
 * The general notation is `[0]y[column id]_[row id]`. A preceeding zero
 * means that the value of the cell is considered to be zero. The column
 * id goes between 0 and the number of columns in the array - 1 (so it
 * is indexed from 0,) and the row id goes from 0 to m_order.
 * 
 * ## Intercept effects
 * 
 * Intercept effects only involve a single set of curly brackets. Using the
 * 'greater-than' symbol (i.e., '<') is only for transition effects. When
 * specifying intercept effects, users can skip the `row_id`, e.g.,
 * `y0_0` is equivalent to `y0`. If the passed `row id` is different from
 * the Markov order, i.e., `row_id != m_order`, then the function returns
 * with an error. 
 * 
 * Examples:
 * 
 * - `"{y0, 0y1}"` is equivalent to set a motif with the first element equal
 * to one and the second to zero. 
 * 
 * ## Transition effects
 * 
 * Transition effects can be specified using two sets of curly brackets and
 * an greater-than symbol, i.e., `{...} > {...}`. The first set of brackets,
 * which we call LHS, can only hold `row id` that are less than `m_order`.
 * 
 * 
 * 
 * @param formula 
 * @param locations 
 * @param signs 
 * @param m_order 
 * @param y_ncol 
 */
inline void defm_motif_parser(
    std::string formula,
    std::vector< size_t > & locations,
    std::vector< bool > & signs,
    size_t m_order,
    size_t y_ncol
)
{
    // Resetting the results
    locations.clear();
    signs.clear();

    std::regex pattern_intercept(
        "\\{\\s*0?y[0-9]+(_[0-9]+)?(\\s*,\\s*0?y[0-9]+(_[0-9]+)?)*\\s*\\}"
        );
    std::regex pattern_transition(
        std::string("\\{\\s*0?y[0-9]+(_[0-9]+)?(\\s*,\\s*0?y[0-9]+(_[0-9]+)?)*\\}\\s*(>)\\s*") +
        std::string("\\{\\s*0?y[0-9]+(_[0-9]+)?(\\s*,\\s*0?y[0-9]+(_[0-9]+)?)*\\s*\\}")
        );

    auto empty = std::sregex_iterator();

    // This column-major vector indicates true if the variable has already been
    // selected
    std::vector< bool > selected((m_order + 1) * y_ncol, false);

    std::smatch match;
    std::regex_match(formula, match, pattern_transition);
    if (!match.empty())
    {

        if (m_order == 0)
            throw std::logic_error("Transition effects are only valid when the data is a markov process.");

        // Will indicate where the arrow is located at
        size_t arrow_position = match.position(4u);

        // This pattern will match 
        std::regex pattern("(0?)y([0-9]+)(_([0-9]+))?");

        auto iter = std::sregex_iterator(formula.begin(), formula.end(), pattern);

        for (auto i = iter; i != empty; ++i)
        {

            // Baseline position
            size_t current_location = i->position(0u);

            // First value true/false
            bool is_positive;
            if (i->operator[](1u).str() == "")
                is_positive = true;
            else if (i->operator[](1u).str() == "0")
                is_positive = false;
            else
                throw std::logic_error("The number preceding y should be either none or zero.");

            // Variable position
            size_t y_col = std::stoul(i->operator[](2u).str());
            if (y_col >= y_ncol)
                throw std::logic_error("The proposed column is out of range.");

            // Time location
            size_t y_row;
            std::string tmp_str = i->operator[](4u).str();
            if (m_order > 1)
            {
                // If missing, we replace with the location 
                if (tmp_str == "")
                {

                    if (current_location > arrow_position)
                        y_row = m_order;
                    else
                        throw std::logic_error("LHS of transition must specify time when m_order > 1");

                } else
                    y_row = std::stoul(tmp_str);

                if (y_row > m_order)
                    throw std::logic_error("The proposed row is out of range.");


            } else {

                // If missing, we replace with the location 
                if (tmp_str != "")
                    y_row = std::stoul(tmp_str);
                else
                    y_row = (current_location < arrow_position ? 0u: 1u);

            }

            if (selected[y_col * (m_order + 1) + y_row])
                throw std::logic_error(
                    "The term " + i->str() + " shows more than once in the formula.");

            // Only the end of the chain can be located at position after the
            // arrow
            if ((current_location > arrow_position) && (y_row != m_order))
                throw std::logic_error(
                    "Only the row " + std::to_string(m_order) +
                    " can be specified at the RHS of the motif."
                    );

            selected[y_col * (m_order + 1) + y_row] = true;

            locations.push_back(y_col * (m_order + 1) + y_row);
            signs.push_back(is_positive);
            

        }

        return;

    } 
    
    std::regex_match(formula, match, pattern_intercept);
    if (!match.empty()){

        // This pattern will match 
        std::regex pattern("(0?)y([0-9]+)(_([0-9]+))?");

        auto iter = std::sregex_iterator(formula.begin(), formula.end(), pattern);

        for (auto i = iter; i != empty; ++i)
        {
            
            // First value true/false
            bool is_positive;
            if (i->operator[](1u).str() == "")
                is_positive = true;
            else if (i->operator[](1u).str() == "0")
                is_positive = false;
            else
                throw std::logic_error("The number preceding y should be either none or zero.");

            // Variable position
            size_t y_col = std::stoul(i->operator[](2u).str());
            if (y_col >= y_ncol)
                throw std::logic_error("The proposed column is out of range.");

            // Time location
            size_t y_row;
            if (i->operator[](4u).str() == "") // Assume is the last
                y_row = m_order;
            else {

                y_row = std::stoul(i->operator[](4u).str());

                if (y_row != m_order)
                    throw std::logic_error(
                        std::string("Intercept motifs cannot feature past events. ") +
                        std::string("Only transition motifs can: {...} > {...}.")
                        );

            }

            if (selected[y_col * (m_order + 1) + y_row])
                throw std::logic_error(
                    "The term " + i->str() + " shows more than once in the formula.");

            selected[y_col * (m_order + 1) + y_row] = true;

            locations.push_back(y_col * (m_order + 1) + y_row);
            signs.push_back(is_positive);
            

        }

        return;

    } 
    
    throw std::logic_error(
        "The motif specified in the formula: " + formula +
        " has the wrong syntax."
        );
    
}
#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry//counters/defm-formula.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/



/**
 * @ingroup counting 
 * @details Details on the available counters for `DEFMworkData` can be found in
 * the \ref counters-network section.
 * 
 */
///@{

/**
 * @brief Data class for DEFM arrays.
 * 
 * This holds information pointing to the data array, including information
 * regarding the number of observations, the time slices of the observation,
 * and the number of covariates in the data.
 * 
 */

class DEFMData;

typedef BArrayDense<int, DEFMData> DEFMArray;

class DEFMData {
public:
    
    DEFMArray * array; // Pointer to the owner of this data
    const double * covariates; ///< Vector of covariates (complete vector)
    size_t obs_start;    ///< Index of the observation in the data.
    size_t X_ncol; ///< Number of columns in the array of covariates.
    size_t X_nrow; ///< Number of rows in the array of covariates.
    std::vector< size_t > covar_sort; /// Value where the sorting of the covariates is stored.
    std::vector< size_t > covar_used; /// Vector indicating which covariates are included in the model
    
    DEFMData() {};
    
    /**
     * @brief Constructor
     * @param covariates_ Pointer to the attribute data.
     * @param obs_start_ Location of the current observation in the covariates
     *  vector
     * @param X_ncol_ Number of columns (covariates.)
     */
    DEFMData(
        DEFMArray * array_,
        const double * covariates_,
        size_t obs_start_,
        size_t X_ncol_,
        size_t X_nrow_
    ) : array(array_), covariates(covariates_), obs_start(obs_start_),
    X_ncol(X_ncol_), X_nrow(X_nrow_) {}; 

    /**
     * @brief Access to the row (i) colum (j) data
     * 
     * @param i 
     * @param j 
     * @return double 
     */
    double operator()(size_t i, size_t j) const;
    double at(size_t i, size_t j) const;
    size_t ncol() const;
    size_t nrow() const;
    void print() const;
    
    ~DEFMData() {};

};

/**
  * @brief Data class used to store arbitrary uint or double vectors */
class DEFMCounterData {
public:
    
    std::vector< size_t > indices;
    std::vector< double > numbers;
    std::vector< bool >   logical;
    
    DEFMCounterData() : indices(0u), numbers(0u) {};
    DEFMCounterData(
        const std::vector< size_t > indices_,
        const std::vector< double > numbers_,
        const std::vector< bool > logical_
    ): indices(indices_), numbers(numbers_), 
        logical(logical_) {};

    size_t idx(size_t i) const {return indices[i];};
    double num(size_t i) const {return numbers[i];};
    bool is_true(size_t i) const {return logical[i];};
    
    ~DEFMCounterData() {};
    
};

class DEFMRuleData {
public:

    std::vector< double > numbers;
    std::vector< size_t > indices;
    std::vector< bool >   logical;

    bool init = false;

    double num(size_t i) const {return numbers[i];};
    size_t idx(size_t i) const {return indices[i];};
    bool is_true(size_t i) const {return logical[i];};

    DEFMRuleData() {};

    DEFMRuleData(
        std::vector< double > numbers_,
        std::vector< size_t > indices_,
        std::vector< bool > logical_
    ) : numbers(numbers_), indices(indices_), logical(logical_) {};

    DEFMRuleData(
        std::vector< double > numbers_,
        std::vector< size_t > indices_
    ) : numbers(numbers_), indices(indices_), logical(numbers_.size()) {};

};

/**
 * @weakgroup rules-phylo Phylo rules
 * @brief Rules for phylogenetic modeling
 * @param rules A pointer to a `PhyloRules` object (`Rules`<`PhyloArray`, `PhyloRuleData`>).
 */
///@{

class DEFMRuleDynData : public DEFMRuleData {
public:
    const std::vector< double > * counts;
    
    DEFMRuleDynData(
        const std::vector< double > * counts_,
        std::vector< double > numbers_ = {},
        std::vector< size_t > indices_ = {},
        std::vector< bool > logical_ = {}
        ) : DEFMRuleData(numbers_, indices_, logical_), counts(counts_) {};
    
    ~DEFMRuleDynData() {};
    
};

/**
 * @name Convenient typedefs for network objects.
 */
///@{
typedef Counter<DEFMArray, DEFMCounterData > DEFMCounter;
typedef Counters<DEFMArray, DEFMCounterData> DEFMCounters;
typedef Support<DEFMArray, DEFMCounterData, DEFMRuleData,DEFMRuleDynData> DEFMSupport;
typedef StatsCounter<DEFMArray, DEFMCounterData> DEFMStatsCounter;
typedef Model<DEFMArray, DEFMCounterData,DEFMRuleData,DEFMRuleDynData> DEFMModel;


typedef Rule<DEFMArray, DEFMRuleData> DEFMRule;
typedef Rules<DEFMArray, DEFMRuleData> DEFMRules;
typedef Rule<DEFMArray, DEFMRuleDynData> DEFMRuleDyn;
typedef Rules<DEFMArray, DEFMRuleDynData> DEFMRulesDyn;



///@}

inline double DEFMData::operator()(size_t i, size_t j) const
{
    return *(covariates + (obs_start + j * X_nrow + i));
}

inline size_t DEFMData::ncol() const {
    return X_ncol;
}

inline size_t DEFMData::nrow() const {
    return X_nrow;
}

inline void DEFMData::print() const {

    for (size_t i = 0u; i < array->nrow(); ++i)
    {

        printf_barry("row %li (%li): ", i, obs_start + i);
        for (size_t j = 0u; j < X_ncol; ++j)
            printf_barry("% 5.2f, ", operator()(i, j));
        printf_barry("\n");
        
    }

}

#define MAKE_DEFM_HASHER(hasher,a,cov) Hasher_fun_type<DEFMArray,DEFMCounterData> hasher = [cov](const DEFMArray & array, DEFMCounterData * d) { \
            std::vector< double > res; \
            /* Adding the column feature */ \
            for (size_t i = 0u; i < array.nrow(); ++i) \
                res.push_back(array.D()(i, cov)); \
            /* Adding the fixed dims */ \
            for (size_t i = 0u; i < (array.nrow() - 1); ++i) \
                for (size_t j = 0u; j < array.ncol(); ++j) \
                    res.push_back(array(i, j)); \
            return res;\
        };
    

/**@name Macros for defining counters
  */
///@{
/**Function for definition of a network counter function*/
#define DEFM_COUNTER(a) \
inline double (a) (const DEFMArray & Array, uint i, uint j, DEFMCounterData & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_COUNTER_LAMBDA(a) \
Counter_fun_type<DEFMArray, DEFMCounterData> a = \
    [](const DEFMArray & Array, uint i, uint j, DEFMCounterData & data) -> double

///@}

/**@name Macros for defining rules
  */
///@{
/**Function for definition of a network counter function*/
#define DEFM_RULE(a) \
inline bool (a) (const DEFMArray & Array, uint i, uint j, bool & data)

/**Lambda function for definition of a network counter function*/
#define DEFM_RULE_LAMBDA(a) \
Rule_fun_type<DEFMArray, DEFMRuleData> a = \
[](const DEFMArray & Array, uint i, uint j, DEFMRuleData & data) -> bool
///@}

/**Lambda function for definition of a network counter function*/
#define DEFM_RULEDYN_LAMBDA(a) \
Rule_fun_type<DEFMArray, DEFMRuleDynData> a = \
[](const DEFMArray & Array, uint i, uint j, DEFMRuleDynData & data) -> bool
///@}

/**
  * @weakgroup  counters-network DEFMArray counters
  * @brief Counters for network models
  * @param counters A pointer to a `DEFMCounters` object (`Counters`<`DEFMArray`, `DEFMCounterData`>).
  */
///@{
// -----------------------------------------------------------------------------
/**
 * @brief Prevalence of ones
 * 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
inline void counter_ones(
    DEFMCounters * counters,
    int covar_index = -1,
    std::string vname = "",
    const std::vector< std::string > * x_names = nullptr
)
{

    // Weighted by a feature of the array
    if (covar_index >= 0)
    {   

        MAKE_DEFM_HASHER(hasher, array, covar_index)

        DEFM_COUNTER_LAMBDA(counter_tmp)
        {

            // Only count the current
            if (i != (Array.nrow() - 1))
                return 0.0;

            return Array.D()(i, data.idx(0u));

        };


        if (vname == "")
        {
            if (x_names != nullptr)
                vname = x_names->operator[](covar_index);
            else
                vname = std::string("attr")+ std::to_string(covar_index);
        }

        counters->add_counter(
            counter_tmp, nullptr, hasher,
            DEFMCounterData({static_cast<size_t>(covar_index)}, {}, {}), 
            "Num. of ones x " + vname, 
            "Overall number of ones"
        );



    } else {

        DEFM_COUNTER_LAMBDA(count_ones)
        {
            
            // Only count the current
            if (i != (Array.nrow() - 1))
                return 0.0;

            return 1.0;
        };

        counters->add_counter(
            count_ones, nullptr, nullptr,
            DEFMCounterData(),
            "Num. of ones", 
            "Overall number of ones"
        );
    }

    return;

}

inline void counter_logit_intercept(
    DEFMCounters * counters,
    size_t n_y,
    std::vector< size_t > which = {},
    int covar_index = -1,
    std::string vname = "",
    const std::vector< std::string > * x_names = nullptr,
    const std::vector< std::string > * y_names = nullptr
) {


    if (which.size() == 0u)
    {
        which.resize(n_y, 0u);
        std::iota(which.begin(), which.end(), 0u);
    } else {
        for (auto w : which)
            if (w >= n_y)
                throw std::logic_error("Values in `which` are out of range.");
    }

    // Case when no interaction happens, whatsoever.
    if (covar_index < 0)
    {

        DEFM_COUNTER_LAMBDA(tmp_counter)
        {
            if (i != (Array.nrow() - 1))
                return 0.0;

            if (j != data.idx(0u))
                return 0.0;

            return 1.0;
        };

        for (auto i : which)
        {

            if (y_names != nullptr)
                vname = y_names->operator[](i);
            else
                vname = std::to_string(i);

            counters->add_counter(
                tmp_counter, nullptr, nullptr,
                DEFMCounterData({i}, {}, {}), 
                "Logit intercept " + vname, 
                "Equal to one if the outcome " + vname + " is one. Equivalent to the logistic regression intercept."
            );

        }

    } else {

        DEFM_COUNTER_LAMBDA(tmp_counter)
        {
            if (i != Array.nrow() - 1)
                return 0.0;

            if (j != data.idx(0u))
                return 0.0;

            return Array.D()(i, data.idx(1u));
        };

        MAKE_DEFM_HASHER(hasher, array, covar_index)
        bool hasher_added = false;

        std::string yname;
        for (auto i : which)
        {

            if (y_names != nullptr)
                yname = y_names->operator[](i);
            else
                yname = std::to_string(i);

            if (vname == "")
            {
                if (x_names != nullptr)
                    vname = x_names->operator[](covar_index);
                else
                    vname = std::string("attr")+ std::to_string(covar_index);
            }

            if (hasher_added)
                counters->add_counter(
                    tmp_counter, nullptr, nullptr,
                    DEFMCounterData({i, static_cast<size_t>(covar_index)}, {}, {}), 
                    "Logit intercept " + yname + " x " + vname, 
                    "Equal to one if the outcome " + yname + " is one. Equivalent to the logistic regression intercept."
                );
            else {

                hasher_added = true;

                counters->add_counter(
                    tmp_counter, nullptr, hasher,
                    DEFMCounterData({i, static_cast<size_t>(covar_index)}, {}, {}), 
                    "Logit intercept " + yname + " x " + vname, 
                    "Equal to one if the outcome " + yname + " is one. Equivalent to the logistic regression intercept."
                );

            }

        }

    }
    

}

/**
 * @brief Prevalence of ones
 * 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
inline void counter_transition(
    DEFMCounters * counters,
    std::vector< size_t > coords,
    std::vector< bool > signs,
    size_t m_order,
    size_t n_y,
    int covar_index = -1,
    std::string vname = "",
    const std::vector< std::string > * x_names = nullptr,
    const std::vector< std::string > * y_names = nullptr
)
{

    // A vector to store the type of dat
    if (signs.size() == 0u)
        signs.resize(coords.size(), true);
    else if (signs.size() != coords.size())
        throw std::length_error("Size of -coords- and -signs- must match.");

    if (covar_index >= 0)
        coords.push_back(static_cast<size_t>(covar_index));
    else
        coords.push_back(1000u);

    DEFM_COUNTER_LAMBDA(count_init)
    {

        auto indices = data.indices;

        for (size_t i = 0u; i < (indices.size() - 1u); ++i)
        {
            if (
                std::floor(indices[i] / Array.nrow()) >= 
                static_cast<int>(Array.ncol())
                )
                throw std::range_error("The motif includes entries out of range.");
        }
            
        return 0.0;
        
    };

    DEFM_COUNTER_LAMBDA(count_ones)
    {
        
        auto dat = data.indices;
        auto sgn = data.logical;
        int covaridx = dat[dat.size() - 1u];

        // Checking if the observation is in the stat. We 
        const auto & array = Array.get_data();
        size_t loc = i + j * Array.nrow();
        size_t n_cells = dat.size() - 1u;

        // Only one currently needs to be a zero for it
        // to change
        size_t n_now = 0;
        bool baseline_value = false;
        bool i_in_array = false;
        for (size_t e = 0u; e < n_cells; ++e)
        {

            // Is the current cell in the list?
            if (dat[e] == loc)
            {
                i_in_array = true;
                baseline_value = sgn[e];
            }

            if ((sgn[e] & (array[dat[e]] == 1)) | (!sgn[e] & (array[dat[e]] == 0)))
                n_now++;
            
        }

        // If i in array still false, then no change
        if (!i_in_array)
            return 0.0;
        
        size_t n_prev = n_now;
        if (baseline_value)
            n_prev--;
        else
            n_prev++;

        // Computing stats
        if (covaridx < 1000)
        {
            
            double val = Array.D()(Array.nrow() - 1u, covaridx);
            double value_now  = n_now == n_cells ?  val : 0.0;
            double value_prev = n_prev == n_cells ? val : 0.0;

            return value_now - value_prev;

        } 
        else
        {

            double value_now  = n_now == n_cells ? 1.0 : 0.0;
            double value_prev = n_prev == n_cells ? 1.0 : 0.0;

            return value_now - value_prev;

        }

    };

    // Creating name of the structure
    std::string name;
    if (coords.size() == 1u)
        name = "";
    else
        name = "Motif ";

    // Creating an empty motif filled with zeros
    barry::BArrayDense<int> motif(m_order + 1u, n_y, 0);

    // Filling the matrix in, negative values are 0s and 1s are... 1s.
    // Zero are values not used.
    size_t n_cells = coords.size() - 1u;
    for (size_t d = 0u; d < n_cells; ++d)
    {
        size_t c = std::floor(coords[d] / (m_order + 1u));
        size_t r = coords[d] - c * (m_order + 1u);
        motif(r, c) = signs[d] ? 1 : -1;
        
    }

    // Checking if any prior to the event
    bool any_before_event = false;
    
    for (size_t i = 0u; i < m_order; ++i)
    {
        for (size_t j = 0u; j < n_y; ++j)
        {
            if (motif(i,j) != 0)
            {
                any_before_event = true;
                break;
            }

        }
    }
    
    #ifdef BARRY_WITH_LATEX
        name += "$";
    #endif

    if (any_before_event)
        #ifdef BARRY_WITH_LATEX
            name += "(";
        #else
            name += "{";
        #endif

    #ifdef BARRY_WITH_LATEX
        #define UNI_SUB(a) \
            (\
                ((a) == 0) ? "_0" : (\
                ((a) == 1) ? "_1" : (\
                ((a) == 2) ? "_2" : (\
                ((a) == 3) ? "_3" : (\
                ((a) == 4) ? "_4" : (\
                ((a) == 5) ? "_5" : (\
                ((a) == 6) ? "_6" : (\
                ((a) == 7) ? "_7" : (\
                ((a) == 8) ? "_8" : \
                "_9"))))))))\
            )
    #else
        #define UNI_SUB(a) \
            (\
                ((a) == 0) ? "\u2080" : (\
                ((a) == 1) ? "\u2081" : (\
                ((a) == 2) ? "\u2082" : (\
                ((a) == 3) ? "\u2083" : (\
                ((a) == 4) ? "\u2084" : (\
                ((a) == 5) ? "\u2085" : (\
                ((a) == 6) ? "\u2086" : (\
                ((a) == 7) ? "\u2087" : (\
                ((a) == 8) ? "\u2088" : \
                "\u2089"))))))))\
            )
    #endif

    // If order is greater than zero, the starting point of the transtion
    for (size_t i = 0u; i < m_order; ++i)
    {

        bool row_start = true;
        for (size_t j = 0u; j < n_y; ++j)
        {

            // Is it included?
            if (motif(i,j) == 0)
                continue;

            // Is not the first?
            if (row_start)
                row_start = false;
            else
                name += ", ";

            if (y_names != nullptr)
                name += y_names->operator[](j);
            else
                name += (std::string("y") + std::to_string(j));

            #ifdef BARRY_WITH_LATEX
                name += (motif(i,j) < 0 ? "^-" : "^+");
            #else
                name += (motif(i,j) < 0 ? "\u207B" : "\u207A");
            #endif

        }

    }

    // If it has starting point, then need to close.
    if (any_before_event & (m_order > 0u))
        #ifdef BARRY_WITH_LATEX
            name += ") -> (";
        #else
            name += "} \u21E8 {";
        #endif
    else
        #ifdef BARRY_WITH_LATEX
            name += "(";
        #else
            name += "{";
        #endif

    // Looking onto the transtions
    bool row_start = true;
    for (size_t j = 0u; j < n_y; ++j)
    {

        if (motif(m_order, j) == 0)
            continue;

        if (row_start)
            row_start = false;
        else
            name += ", ";

        if (y_names != nullptr)
            name += y_names->operator[](j);
        else
            name += (std::string("y") + std::to_string(j));

        #ifdef BARRY_WITH_LATEX
        name += (motif(m_order, j) < 0 ? "^-" : "^+" );
        #else
        name += (motif(m_order, j) < 0 ? "\u207B" : "\u207A" );
        #endif


    }

    #undef UNI_SUB

    #ifdef BARRY_WITH_LATEX
    name += ")$";
    #else
    name += "}";
    #endif

    if (covar_index >= 0)
    {

        MAKE_DEFM_HASHER(hasher, array, covar_index)

        if (vname == "")
        {
            if (x_names != nullptr)
                vname = x_names->operator[](covar_index);
            else
                vname = std::string("attr")+ std::to_string(covar_index);
        }

        counters->add_counter(
            count_ones, count_init, hasher,
            DEFMCounterData(coords, {}, signs), 
            name + " x " + vname, 
            "Motif weighted by single attribute"
        );

    } else {

        counters->add_counter(
            count_ones, count_init, nullptr,
            DEFMCounterData(coords, {}, signs), 
            name, 
            "Motif"
        );

    }
    

    return;

}

/**
 * @brief Prevalence of ones
 * 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
inline void counter_transition_formula(
    DEFMCounters * counters,
    std::string formula,
    size_t m_order,
    size_t n_y,
    int covar_index = -1,
    std::string vname = "",
    const std::vector< std::string > * x_names = nullptr,
    const std::vector< std::string > * y_names = nullptr
) {

    std::vector< size_t > coords;
    std::vector< bool > signs;

    defm_motif_parser(
        formula, coords, signs, m_order, n_y
    );

    counter_transition(
        counters, coords, signs, m_order, n_y, covar_index, vname,
        x_names, y_names
    );

}

/**
 * @brief Prevalence of ones
 * 
 * @param counters Pointer ot a vector of counters
 * @param covar_index If >= than 0, then the interaction
 */
inline void counter_fixed_effect(
    DEFMCounters * counters,
    int covar_index,
    double k,
    std::string vname = "",
    const std::vector< std::string > * x_names = nullptr
)
{

    DEFM_COUNTER_LAMBDA(count_init)
    {
        return std::pow(Array.D()((size_t) i, data.idx(0u)), data.num(0u));
    };

    DEFM_COUNTER_LAMBDA(count_tmp)
    {
        return 0.0;
    };

    MAKE_DEFM_HASHER(hasher, array, covar_index)

    if (x_names != nullptr)
        vname = x_names->operator[](covar_index);
    else
        vname = std::string("attr")+ std::to_string(covar_index);

    counters->add_counter(
        count_tmp, count_init, hasher,
        DEFMCounterData({static_cast<size_t>(covar_index)}, {k}, {}), 
        "Fixed effect feature (" + vname + ")^" + std::to_string(k)
    );

    return;

}

/**
 * @name Returns true if the cell is free
 * @param rules A pointer to a `DEFMRules` object (`Rules`<`DEFMArray`, `bool`>).
 */
///@{
// -----------------------------------------------------------------------------
/**@brief Number of edges */
inline void rules_markov_fixed(
    DEFMRules * rules,
    size_t markov_order
    ) {
    
    DEFM_RULE_LAMBDA(no_self_tie) {
        return i >= data.idx(0u);
    };
    
    rules->add_rule(
        no_self_tie,
        DEFMRuleData({},{markov_order})
        );
    
    return;
}

/**
 * @brief Blocks switching a one to zero.
 * 
 * @param rules 
 * @param ids Ids of the variables that will follow this rule.
 */
inline void rules_dont_become_zero(
    DEFMSupport * support,
    std::vector<size_t> ids
    ) {
    
    
    DEFM_RULE_LAMBDA(rule) {

        if (!data.init)
        {
            std::vector< size_t > tmp(Array.ncol(), 0u);

            for (auto v : data.indices)
            {
                if (v >= Array.ncol())
                    throw std::range_error("The specified id for `dont_become_zero` is out of range.");

                tmp[v] = 1u;
            }

            data.indices.resize(Array.ncol());
            for (size_t v = 0u; v < tmp.size(); ++v)
                data.indices[v] = tmp[v];

            data.init = true;
        }

        // If not considered, then continue
        if (data.indices[j] == 0u)
            return true;

        // The last observation is always included
        if (i == (Array.nrow() - 1))
            return true;

        // This is now one, is the next different zero? If so,
        // we can include it (1->1)
        return (Array(i + 1, j) != 0); // |
            // (Array(i, j) != 1);

    };
    
    support->get_rules()->add_rule(
        rule,
        DEFMRuleData({}, {ids})
        );
    
    return;
}

/**
 * @brief Blocks switching a one to zero.
 * 
 * @param rules 
 * @param ids Ids of the variables that will follow this rule.
 */
inline void rules_exclude_all_ones(
    DEFMSupport * support
    ) {
    

    DEFM_RULEDYN_LAMBDA(rule) {

        if (!data.init)
        {
            data.init = true;
            data.indices[0u] = (Array.nrow() * Array.ncol());
        }

        return Array.nnozero() != data.idx(0u);

    };
    
    support->get_rules_dyn()->add_rule(
        rule,
        DEFMRuleDynData(nullptr, {}, {0u})
        );
    
    return;
}

///@}

///@}

#endif
/*//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

 End of -include/barry/counters/defm.hpp-

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/


        }
    }
    
}

namespace netcounters = barry::counters::network;
namespace phylocounters = barry::counters::phylo;
namespace defmcounters = barry::counters::defm;

#define COUNTER_FUNCTION(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
    inline double (a) (const Array_Type & Array, uint i, uint j, Data_Type & data)\

#define COUNTER_LAMBDA(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
    Counter_fun_type<Array_Type, Data_Type> a = \
    [](const Array_Type & Array, uint i, uint j, Data_Type & data)

#define RULE_FUNCTION(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
    inline bool (a) (const Array_Type & Array, uint i, uint j, Data_Type & data)\

#define RULE_LAMBDA(a) template <typename Array_Type = barry::BArray<>, typename Data_Type = bool> \
    Rule_fun_type<Array_Type, Data_Type> a = \
    [](const Array_Type & Array, uint i, uint j, Data_Type & data)

#endif
