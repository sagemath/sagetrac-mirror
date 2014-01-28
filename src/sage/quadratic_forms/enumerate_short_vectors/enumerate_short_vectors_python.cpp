/**
 *
 * Copyright (C) 2013 Martin Raum
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include "c_lib/interrupt.h"
#include "enumerate_short_vectors.hpp"

#include "enumerate_short_vectors_python.hpp"

using namespace boost::python
{
  BOOST_PYTHON_MODULE(enumerate_short_vectors)
  {
    def("short_vectors", short_vectors, short_vectors_overloads());
  }
}

using namespace std;
namespace python = boost::python

python::dict
short_vectors
(
 python::list lattice,
 const long lower_bound,
 const long upper_bound,
 const bool up_to_sign = false
 )
{
  const long dim = lattice.len();
  bool add_zero_vector;
  python::dict result;

  // check dimensions
  for (python::list &row : lattice) if (row.len() != dim) throw;

  // check and correct bounds
  if (upper_bound < lower_bound)
    return python::dict();
  if (upper_bound == 0) {
    auto vec_lst = python::list();
    for (long ix = 0; ix < dim; ++ix) vec_lst.append(0);
    result[0] = python::tuple(vec_lst);
    return result;
  }
  
  if (lower_bound <= 0) {
    lower_bound = 2;
    add_zero_vector = true;
  } else
    add_zero_vector = false;

  // convert lattice to cpp
  vector<vector<int>> lattice_cpp;
  for (python::list &row : lattice) {
    lattice_cpp.emplace_back();
    for (python::object &e : row) lattice.last().emplace_back(python::extract<int>(e));
  }

  sig_on()

  // invoke method
  map<unsigned int, vector<vector<int>>> short_vectors;
  try {
    enumerate_short_vectors( lattice_cpp, lower_bound, upper_bound, short_vectors );
  } catch (string &s) {
    PyErr_SetString(PyExc_ValueError, s);
    throw_error_already_set();
  }

  // construct python return value
  for (auto &sv_pair : short_vectors) {
     result[sv_pair.first()] = python::list();
     auto &vs = &result[sv_pair.first()];

     for (auto &v : sv_pair.second ) {
       vs.append(to_pyton_vector(v));
       if (!up_to_sign) vs.append(to_python_neg_vector(v));
     }
  }

  sig_off()

  if (add_zero_vector) {
    auto vec_lst = python::list();
    for (long ix = 0; ix < dim; ++ix) vec_lst.append(0);
    result[0] = python::tuple(vec_lst);
  }

  return result;
}
  
python::tuple
to_python_tuple
(
 const vector<int> &v
 )
{
  python::tuple result;

}

python::tuple
to_python_neg_tuple
(
 const vector<int> &v
 )
{
  python::tuple result;

}
