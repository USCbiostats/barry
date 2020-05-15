#include <vector>
#include <unordered_map>

#ifndef BARRAY_META_HPP
#define BARRAY_META_HPP 1

class Meta {
public:
  Meta() : data(0u) {};
  ~Meta() {};
  std::unordered_map< std::string, bool > data;
  bool is(const std::string & attr) const;
  void set(const std::string & attr, bool value);
};

inline bool Meta::is(const std::string & attr) const {
  
  auto iter = data.find(attr);
  if (iter != data.end())
    return iter->second; 
  
  return false;
}

inline void Meta::set(const std::string & attr, bool value) {
  data[attr] = value;
  return;
}

#endif