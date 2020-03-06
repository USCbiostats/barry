#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double hash_it(double x) {
  std::size_t H = std::hash< double >{}(x);
  std::cout << "This is the hash: " << H << "\n";
  return (double) H;
}


// [[Rcpp::export]]
double vector2hash(std::vector< double > x) {
  
  std::hash< double > hasher;
  std::size_t hash = hasher(x[0]);

  if (x.size() > 1u)
    for (int i = 1; i < x.size(); ++i)
      hash ^= hasher(x[i]) + 0x9e3779b9 + (hash<<6) + (hash>>2);
  
  std::cout << "This is the hash: " << hash << "\n";  
  return (double) hash;
}

// Here I'm creating a hash function to be used with unordered_map

// A simple pair structure of integers
struct coords{
  int i;
  int j;
};

// Needs to have a way to be compared
bool operator==(const coords& lhs, const coords& rhs) {
  return lhs.i == rhs.i && lhs.j == rhs.j;
}

namespace std {
  template<> struct hash<coords> {
    std::size_t operator()(coords const&  dat) const noexcept {
      
      std::hash<int> hasher;
      std::size_t hash = hasher(dat.i);
      
      return hash ^ (hasher(dat.j) + 0x9e3779b9 + (hash << 6) + (hash >> 2));
      
    }
  };
}

struct CoordsHasher {
  std::size_t operator()(coords const&  dat) const noexcept {
    
    std::hash<int> hasher;
    std::size_t hash = hasher(dat.i);
    
    // ^ makes bitwise XOR
    // 0x9e3779b9 is a 32 bit constant (comes from the golden ratio)
    // << is a shift operator, something like lhs * 2^(rhs)
    return hash ^ (hasher(dat.j) + 0x9e3779b9 + (hash << 6) + (hash >> 2));
    
  }
};


// [[Rcpp::export]]
int test(int i, int j) {
  
  // Alt 1
  std::unordered_map< coords, double > M1;
  
  // Alt 2
  std::unordered_map< coords, double, CoordsHasher > M2;
  
  coords c0 = {i, j};
  coords c1 = {i - 10, j + 10};
  
  // Writing
  M1[c0] = (double) (pow(c0.i, 2) + pow(c0.j, 2));
  M1[c1] = (double) (pow(c1.i, 2) + pow(c1.j, 2));
  
  M2[c0] = (double) (pow(c0.i, 2) + pow(c0.j, 2));
  M2[c1] = (double) (pow(c1.i, 2) + pow(c1.j, 2));
  
  // Recovering
  std::cout << "For case 1:\n" ;
  std::cout << "For values " << c0.i << ", " << c0.j << " the value is " << M1.at(c0) <<"\n";
  std::cout << "For values " << c1.i << ", " << c1.j << " the value is " << M1.at(c1) <<"\n";
  
  std::cout << "For case 2:\n" ;
  std::cout << "For values " << c0.i << ", " << c0.j << " the value is " << M2.at(c0) <<"\n";
  std::cout << "For values " << c1.i << ", " << c1.j << " the value is " << M2.at(c1) <<"\n";
  
  return 0;
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
z <- hash_it(500.11)

h1 <- vector2hash(c(1,2,1,.9))
h2 <- vector2hash(sort(c(1,2,1,.9)))
h3 <- vector2hash(c(1,2,1,.9) + 1)
h4 <- vector2hash(c(1,2,1,.9))

test(1, 5)
test(1, 3)
*/
