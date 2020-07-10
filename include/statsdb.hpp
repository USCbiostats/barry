// #include <vector>
// #include <functional>
// #include <unordered_map>
#include "typedefs.hpp"

#ifndef LBARRAY_STATSDB_HPP 
#define LBARRAY_STATSDB_HPP 1
 
/**@brief Database of statistics.
 * 
 * This is mostly used in `Support`.
 * 
 */
template<typename T = double> 
class FreqTable {
private:
  MapVec_type<T, uint> data;
  
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
  // void rehash();
  
  
};

template<typename T>  
inline void FreqTable<T>::add(const std::vector< T > & x) { 
  
  // The term exists, then we add it to the list and we initialize it
  // with a single count
  if (data.find(x) == data.end()) {
    data[x] = 1u;
  } else // We increment the counter
    data.at(x)++;
  
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

#endif
