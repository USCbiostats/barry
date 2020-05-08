#include <Rcpp.h>
using namespace Rcpp;




// Definition of the class structure

// Edgelist
typedef unsigned int uint;
typedef std::unordered_map< uint, double > umap_int_dbl;

class EdgeList {
public:
  uint N;
  uint M;
  std::vector< umap_int_dbl > EL;
  
  // Empty datum
  EdgeList (uint N_, uint M_) : N(N_), M(M_), EL(N_) {};
  
  // Edgelist with data
  EdgeList (
      uint N_, uint M_,
      const std::vector< uint > & source,
      const std::vector< uint > & target,
      const std::vector< double > & value,
      bool add = true
    ) {
    
    if (source.size() != target.size())
      Rcpp::stop("Must match the size.");
    if (source.size() != value.size())
      Rcpp::stop("Must match the size.");
    
    // Initializing
    N = N_;
    M = M_;
    EL.resize(N);
    
    // Writing the data
    for (uint i = 0u; i < source.size(); ++i) {
      
      // Checking range
      if (source.at(i) >= N_ | target.at(i) >= M_)
        Rcpp::stop("Out of range.");
      
      // Checking if it exists
      auto search = EL.at(source.at(i)).find(target.at(i));
      if (search != EL.at(source.at(i)).end()) {
        if (!add)
          Rcpp::stop("The value already exists");
        
        EL.at(source.at(i))[target.at(i)] = search->second + value.at(i);
        continue;
      }
        
      // Adding the value
      EL.at(source.at(i))[target.at(i)] = value.at(i);
    }
    
    return;
    
  }
  
  // Function to access the elements
  double get_cell(uint i, uint j) const {
    
    if (this->EL.at(i).size() == 0u)
      return 0.0;
    
    // If it is not empty, then find and return
    auto search = EL.at(i).find(j);
    if (search != EL.at(i).end())
      return search->second;
    
    // This is if it is empty
    return 0.0;
    
  }
  
  umap_int_dbl get_row(uint i) const {
    return EL.at(i);
  }
};


// [[Rcpp::export]]
SEXP new_EdgeList(
    int N, int M,
    const std::vector< uint > & source,
    const std::vector< uint > & target,
    const std::vector< double > & value
) {
  
  Rcpp::XPtr< EdgeList > ptr(
    new EdgeList((uint) N, (uint) M, source, target, value),
    true
  );
  
  return ptr;
  
}

// [[Rcpp::export]]
double get_cell(SEXP x, int i, int j) {
  Rcpp::XPtr< EdgeList > xptr(x);
  return xptr->get_cell(i, j);
}

// Returning the i-th row
// [[Rcpp::export]]
NumericVector get_row(SEXP x, int i) {
  
  Rcpp::XPtr< EdgeList > xptr(x);
  NumericVector ans(xptr->M, 0);
  umap_int_dbl m = xptr->get_row(i);
  
  for (auto row = m.begin(); row != m.end(); ++row)
    ans[row->first] = row->second;
  
  return ans;
  
}

/***R

set.seed(123)
N <- 10000
M <- 20000

nedges  <- 1e4
source <- sample.int(N, nedges, replace = TRUE)
target <- sample.int(M, nedges, replace = TRUE)
values <- runif(nedges)

el <- new_EdgeList(N, M, source - 1L, target - 1L, values)
cbind(cbind(source, target) - 1, values)

ans <- sapply(1:nedges, function(i) get_cell(el, source[i]-1, target[i]-1))
range(values - ans)

rowz <- get_row(el, 1)
mat <- matrix(nrow = N, ncol = M)

*/
