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
    MapVec_type<T, uint> data;
    
    /**
    * A nasty way to declare when using :: on template members
    * https://stackoverflow.com/questions/6571381/dependent-scope-and-nested-templates/6571836
    */
    typename MapVec_type<T,uint>::iterator iter;
        
public:
    // uint ncols;
    FreqTable() {};
    ~FreqTable() {};
    
    void add(const std::vector< T > & x);
    
    Counts_type                 as_vector() const;
    MapVec_type<T,uint>         get_data() const;
    const MapVec_type<T,uint> * get_data_ptr() const;
    
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
    iter = data.find(x);
    if (iter == data.end()) {
        data.insert(std::make_pair(x, 1u));
    } else // We increment the counter
        ++(iter->second);
    
    return; 
}

template<typename T>
inline Counts_type FreqTable<T>::as_vector() const { 
    
    Counts_type ans;
    ans.reserve(data.size());
    for (auto iter = data.begin(); iter != data.end(); ++iter)
        ans.push_back(*iter);
    
    
    return ans;
}

template<typename T>
inline MapVec_type<T,uint> FreqTable<T>::get_data() const {
    return data;
}

template<typename T>
inline const MapVec_type<T,uint> * FreqTable<T>::get_data_ptr() const {
    return &data;
}

template<typename T>
inline void FreqTable<T>::clear() {
    data.clear();
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

    uint grand_total = 0u;
    printf("%7s | %s\n", "Counts", "Stats");
    for (auto i = data.begin(); i != data.end(); ++i)
    {

        printf("%7i | ", i->second);

        for (const auto& j : i->first)
            printf(" %.2f", j);
        printf("\n");

        grand_total += i->second;

    }

    printf("Grand total: %i\n", grand_total);

    return;

}

template<typename T>
inline size_t FreqTable<T>::size() const noexcept {

    return this->data.size();

}

#endif
