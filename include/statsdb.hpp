// #include <vector>
// #include <functional>
// #include <unordered_map>
#include "typedefs.hpp"

#ifndef LBARRAY_STATSDB_HPP 
#define LBARRAY_STATSDB_HPP 1

template <typename T>
struct vecHasher {
  std::size_t operator()(std::vector< T > const&  dat) const noexcept {
    
    std::hash< T > hasher;
    std::size_t hash = hasher(dat[0u]);
    
    // ^ makes bitwise XOR
    // 0x9e3779b9 is a 32 bit constant (comes from the golden ratio)
    // << is a shift operator, something like lhs * 2^(rhs)
    if (dat.size() > 1u)
      for (unsigned int i = 1u; i < dat.size(); ++i)
        hash ^= hasher(dat[i]) + 0x9e3779b9 + (hash<<6) + (hash>>2);
    
    return hash;
    
  }
};

/**@brief Database of statistics.
 * 
 * This is mostly used in `Support`.
 * 
 */
class StatsDB {
private:
  std::unordered_map< std::vector< double >, uint, vecHasher< double > > stats;
  
public:
  // uint ncols;
  StatsDB() {};
  ~StatsDB() {};
  
  void add(const std::vector< double > & x) {
    
    // The term exists, then we add it to the list and we initialize it
    // with a single count
    if (stats.find(x) == stats.end()) {
      stats[x] = 1u;
    } else // We increment the counter
      stats.at(x)++;
   
    return; 
  };
  
  Counts_type get_entries() const { 
    
    Counts_type ans;
    ans.reserve(stats.size());
    for (auto iter = stats.begin(); iter != stats.end(); ++iter)
      ans.push_back(*iter);
    
    
    return ans;
  };
  
  void clear();
  
  
};

inline void StatsDB::clear() {
  stats.clear();
  return;
}


#endif
