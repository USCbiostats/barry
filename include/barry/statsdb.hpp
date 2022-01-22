// #include <vector>
// #include <functional>
// #include <unordered_map>
#include "typedefs.hpp"

#ifndef BARRY_STATSDB_HPP 
#define BARRY_STATSDB_HPP 1
  
/**
 * @brief Database of statistics.
 * 
 * This is mostly used in `Support`.
 * 
 */
template<typename T = double> 
class FreqTable {
private:

    MapVec_type<T, int> index;
    std::vector< double > data;
    size_t k = 0u;
    size_t n = 0u;

    typename MapVec_type<T,int>::iterator iter;
        
public:
    // uint ncols;
    FreqTable() {};
    ~FreqTable() {};
    
    void add(const std::vector< T > & x);
    
    Counts_type                 as_vector() const;
    const std::vector< double > & get_data() const {return data;};
    const MapVec_type<T,int> & get_index() const {return index;};
    // const MapVec_type<T,uint> * get_data_ptr() const;
    
    void clear();
    void reserve(unsigned int n);
    void print() const;
    size_t size() const noexcept;
    // void rehash();
    
    
};

template<typename T>  
inline void FreqTable<T>::add(const std::vector< T > & x) { 
    
    // The term exists, then we add it to the list and we initialize it
    // with a single count
    if (k == 0u)
    {

        index[x] = 0u;

        data.push_back(1.0);
        data.insert(data.end(), x.begin(), x.end());

        k = x.size();

        n++;

        return;

    }
    else
    {

        if (x.size() != k)
            throw std::length_error(
                "The value you are trying to add doesn't have the same lenght used in the database."
                );
            
        iter = index.find(x);

        if (iter == index.end())
        {

            index[x] = data.size();

            data.push_back(1.0);

            data.insert(data.end(), x.begin(), x.end());

            n++;
            
            return;

        }

        data[(*iter).second] += 1.0;

    }
    
    return; 

}

template<typename T>
inline Counts_type FreqTable<T>::as_vector() const { 
    
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

// template<typename T>
// inline MapVec_type<T,uint> FreqTable<T>::get_data() const {
//     return data;
// }

// template<typename T>
// inline const MapVec_type<T,uint> * FreqTable<T>::get_data_ptr() const {
//     return &data;
// }

template<typename T>
inline void FreqTable<T>::clear() {
    index.clear();
    data.clear();
    n = 0u;
    k = 0u;
    return;
}

template<typename T>
inline void FreqTable<T>::reserve(
    unsigned int n
) {
    data.reserve(n);
    return;
}

// inline void StatsDB::rehash() {
//   stats.rehash();
//   return;
// }

template<typename T>
inline void FreqTable<T>::print() const {

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
inline size_t FreqTable<T>::size() const noexcept {

    return n;

}

#endif
