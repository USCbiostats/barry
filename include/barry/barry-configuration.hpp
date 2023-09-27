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

#ifdef __OPENMP
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