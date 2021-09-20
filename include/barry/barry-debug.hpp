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