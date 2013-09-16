#ifndef SAGE_CHOMP_MATRIXCOMPLEX_H
#define SAGE_CHOMP_MATRIXCOMPLEX_H

#include <stdint.h>
#include <vector>
#include <string>

#include "boost/unordered_map.hpp"

#include "chomp/Complex.h"
#include "chomp/Chain.h"

using chomp::Index;

class SageMatrixComplex : public chomp::Complex 
{
 public:
  typedef int64_t Integer;

  CHOMP_COMPLEX(Integer);
  
  virtual void boundary(chomp::Chain* output, const chomp::Index input, int dim) const;
  virtual void coboundary(chomp::Chain * output, const chomp::Index input, int dim) const;
  
 private:
  std::vector< boost::unordered_map<chomp::Index, chomp::Chain> > boundary_data;

 public:
  SageMatrixComplex(const std::vector<std::vector<Integer> >& data);
  
  //! Return a human-readable string representation for the boundary in dimension dim
  std::string to_string(const int dim) const;
};


#endif

