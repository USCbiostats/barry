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
    #define BARRY_MAX_NUM_ELEMENTS static_cast< unsigned int >(UINT_MAX/2u)
#endif

#ifndef BARRY_DEBUG_LEVEL
    #define BARRY_DEBUG_LEVEL 0
#else
    #define BARRY_DEBUG_HEADER "[barry]"
    #define BARRY_DEBUG_MSG(a) \
        printf_barry("%s %s\n", BARRY_DEBUG_HEADER, (a));
    template <typename T>
    void BARRY_DEBUG_VEC_PRINT(const std::vector<T> & a) {
        printf_barry("  [");
        for(const auto & iter : (a)) 
            printf_barry("%.4f ", static_cast<double>(iter));
        printf_barry("]\n");
        return;
    };
#endif

///@}

#endif