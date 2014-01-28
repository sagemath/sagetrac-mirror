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

#include <boost/python/errors.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/module.hpp>
#include <boost/python/stl_iterator.hpp>
#include "interrupt.h"
#include "enumerate_short_vectors.hpp"

#include "enumerate_short_vectors_boost_python.hpp"

using namespace std;
namespace python = boost::python;

BOOST_PYTHON_MODULE(enumerate_short_vectors_boost_python)
{
  python::def("enumerate_short_vectors_boost_python", enumerate_short_vectors_boost_python);
}

python::dict
enumerate_short_vectors_boost_python
(
 python::list lattice,
 long lower_bound,
 long upper_bound,
 bool up_to_sign = false
 )
{
  const long dim = python::len(lattice);
  bool add_zero_vector;
  python::dict result;

  // check dimensions
  python::stl_input_iterator<python::list> lattice_end;
  for (auto it = python::stl_input_iterator<python::list>(lattice);
       it != lattice_end; ++it)
    if (python::len(*it) != dim) {
      PyErr_SetString(PyExc_ValueError, "Lattice must represent a square matrix");
      python::throw_error_already_set();
    }

  // check and correct bounds
  if (upper_bound < lower_bound)
    return python::dict();
  if (upper_bound == 0) {
    result[0] = to_python_tuple(vector<int>(dim, 0));
    return result;
  }

  if (lower_bound <= 0) {
    lower_bound = 2;
    add_zero_vector = true;
  } else
    add_zero_vector = false;

  // convert lattice to cpp
  vector<vector<int>> lattice_cpp;
  for (auto it = python::stl_input_iterator<python::list>(lattice);
       it != lattice_end; ++it) {
    python::stl_input_iterator<int> row_begin(*it), row_end;
    lattice_cpp.emplace_back();
    for (auto row_it = row_begin; row_it != row_end; ++row_it)
      lattice_cpp.back().emplace_back(*row_it);
  }


  sig_on();

  // invoke method
  map<unsigned int, vector<vector<int>>> short_vectors;
  try {
    enumerate_short_vectors( lattice_cpp, lower_bound, upper_bound, short_vectors );
  } catch (string &s) {
    PyErr_SetString(PyExc_ValueError, s.c_str());
    python::throw_error_already_set();
  }

  // construct python return value
  for (auto &sv_pair : short_vectors) {
    python::list vs;
     for (auto &v : sv_pair.second ) {
       vs.append(to_python_tuple(v));
       if (!up_to_sign) vs.append(to_python_neg_tuple(v));
     }
     result[sv_pair.first] = vs;
  }
  if (add_zero_vector) {
    python::list vs;
    vs.append( to_python_tuple(vector<int>(dim, 0)) );
    result[0] = vs;
  }

  sig_off();


  return result;
}
  
python::tuple
to_python_tuple
(
 const vector<int> &v
 )
{
  python::list result;
  for (auto e : v) result.append(e);
  return python::tuple(result);
}

python::tuple
to_python_neg_tuple
(
 const vector<int> &v
 )
{
  python::list result;
  for (auto e : v) result.append(-e);
  return python::tuple(result);
}
