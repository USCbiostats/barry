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
#define BARRY_NCORES_ARG(default) size_t 
#endif


#endif