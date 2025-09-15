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
<<<<<<< HEAD
#define BARRY_VERSION_MINOR 1
#define BARRY_VERSION_PATCH 0

static const size_t barry_version_mayor = BARRY_VERSION_MAYOR;
static const size_t barry_version_minor = BARRY_VERSION_MINOR;
static const size_t barry_version_patch = BARRY_VERSION_PATCH;

#define BARRY_VERSION BARRY_VERSION_MAYOR ## . ## BARRY_VERSION_MINOR ## . ## BARRY_VERSION_PATCH
=======
#define BARRY_VERSION_MINOR 2
#define BARRY_VERSION_PATCH 0
#define BARRY_VERSION BARRY_VERSION_MAYOR ## . ## BARRY_VERSION_MINOR ## . ## BARRY_VERSION_PATCH

static const int barry_version_major = BARRY_VERSION_MAYOR;
static const int barry_version_minor = BARRY_VERSION_MINOR;
static const int barry_version_patch = BARRY_VERSION_PATCH;
>>>>>>> master

/**
  * @brief barry: Your go-to motif accountant
  */
namespace barry {
    
    //! Tree class and TreeIterator class
    #include "typedefs.hpp"
    #include "barry-macros.hpp"
    #include "freqtable.hpp"

    #include "cell-bones.hpp"
    #include "cell-meat.hpp"

    #include "barray-bones.hpp"
    #include "barraycell-bones.hpp"
    #include "barray-meat.hpp"
    #include "barraycell-meat.hpp"
    #include "barray-meat-operators.hpp"

    #include "barraydense-bones.hpp"
    #include "barraydensecell-bones.hpp"

    #include "barraydenserow-bones.hpp"
    #include "barraydensecol-bones.hpp"

    #include "barraydense-meat.hpp"
    #include "barraydensecell-meat.hpp"
    #include "barraydense-meat-operators.hpp"
    
    #include "counters-bones.hpp"
    #include "counters-meat.hpp"

    #include "statscounter-bones.hpp"
    #include "statscounter-meat.hpp"

    #include "support-bones.hpp"
    #include "support-meat.hpp"

    #include "powerset-bones.hpp"
    #include "powerset-meat.hpp"

    #include "model-bones.hpp"
    #include "model-meat.hpp"
    
    #include "rules-bones.hpp"
    #include "rules-meat.hpp"
    
    namespace counters {
        namespace network {
            #include "counters/network.hpp"
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