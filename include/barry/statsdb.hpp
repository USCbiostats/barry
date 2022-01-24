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

    std::unordered_map<size_t, size_t> index;
    std::vector< double > data;
    size_t k = 0u;
    size_t n = 0u;

    typename std::unordered_map<size_t, size_t>::iterator iter;
        
public:
    // uint ncols;
    FreqTable() {};
    ~FreqTable() {};
    
    void add(const std::vector< T > & x);
    
    Counts_type                 as_vector() const;
    const std::vector< double > & get_data() const {return data;};
    const std::unordered_map<size_t,size_t> & get_index() const {return index;};
    
    void clear();
    void reserve(unsigned int n);
    void print() const;

    /**
     * @brief Number of unique elements in the table.
     * (
     * @return size_t 
     */
    size_t size() const noexcept;

    size_t make_hash(const std::vector< double > & x) const;
    
};

template<typename T>  
inline void FreqTable<T>::add(const std::vector< T > & x) { 
    
    // The term exists, then we add it to the list and we initialize it
    // with a single count
    if (k == 0u)
    {

        index.insert({make_hash(x), 0u});

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
        
        size_t h = make_hash(x);
        iter = index.find(h);

        if (iter == index.end())
        {

            index.insert({h, data.size()});
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
    unsigned int n
)
{

    data.reserve(n);

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
inline size_t FreqTable<T>::make_hash(const std::vector< double > & x) const
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
