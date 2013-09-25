#include "sage_matrix_complex.h"

#include <sstream>

using chomp::Complex;
using chomp::Chain;
using chomp::Index;
using chomp::Term;
using std::vector;
using std::string;
using std::stringstream;


SageMatrixComplex::SageMatrixComplex
(const vector<vector<SageMatrixComplex::Cell> >& data)
{
  int current_dimension = 0;
  boundary_data  .resize(current_dimension + 2);
  coboundary_data.resize(current_dimension + 2);
  for (vector<vector<Cell> >::const_iterator ri = data.begin(); ri != data.end(); ri++) {
    const vector<Cell>& row = *ri;
    if (row.size() == 0) {
    } else if (row.size() == 1) {
      current_dimension = row[0];
      if ( current_dimension >= boundary_data.size()) {
        boundary_data  .resize(current_dimension + 2);
        coboundary_data.resize(current_dimension + 2);
      }
    } else if (row.size() == 3) {
      insertCell(row[0], current_dimension);
      insertCell(row[1], current_dimension + 1);
      const Index bdindex  = cellToIndex(row[0], current_dimension);
      const Index cbdindex = cellToIndex(row[1], current_dimension + 1);
      Chain& bd = boundary_data[current_dimension+1][cbdindex];
      bd += Term(bdindex, row[2]);
      bd.dimension()  = current_dimension;
      Chain& cbd = coboundary_data[current_dimension][bdindex];
      cbd += Term(cbdindex, row[2]);
      cbd.dimension() = current_dimension + 1;
    }
  }
}


inline void SageMatrixComplex::boundary(Chain* output, const Index input, int dim) const 
{
  const boost::unordered_map<Index, Chain>& bdry = boundary_data[dim];
  const boost::unordered_map<Index, Chain>::const_iterator 
    found = bdry.find(input);
  if (found != bdry.end())
    *output = found->second;
}


inline void SageMatrixComplex::coboundary(Chain* output, const Index input, int dim) const 
{
  const boost::unordered_map<Index, Chain>& cbdry = coboundary_data[dim];
  const boost::unordered_map<Index, Chain>::const_iterator 
    found = cbdry.find(input);
  if (found != cbdry.end())
    *output = found->second;
}
  

string SageMatrixComplex::to_string(const int dim) const
{
  stringstream ss("");
  for (Index i = 0; i < size(dim); i++) {
    ss << "bd(" << i << ", " << dim << ") = " << boundary(i, dim) << "\n";
  }
  return ss.str();
}


bool SageMatrixComplex::is_chain_complex() const
{
  for (int dim = 0; dim <= dimension(); dim++)
    for (Index i = 0; i < size(dim); i++) {
      Chain dd = simplify(boundary(boundary(i, dim)));
      // std::cout << "bd(bd(" << i << ", " << dim << ")) = " << dd << "\n";
      if (dd().size() != 0)
        return false;
    }
  return true;
}
