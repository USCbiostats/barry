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
 * - `BARRY_CHECK_FINITE` When specified, it will introduce a macro
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

#ifdef BARRY_CHECK_FINITE
  #define BARRY_ISFINITE(a) {if (!std::isfinite( (a) )) \
    throw std::overflow_error("The likelihood function has overflowed.");}
#else
  #define BARRY_ISFINITE(a) 
#endif
///@}

#endif