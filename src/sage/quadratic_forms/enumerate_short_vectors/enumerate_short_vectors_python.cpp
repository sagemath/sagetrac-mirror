// enumerate_short_vectors__python.cc --- 

// Copyright (C) 2014 Martin Raum

// Author: Martin Raum

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

using namespace boost::python
{
  BOOST_PYTHON_MODULE(enumerate_short_vectors)
  {
    def("short_vectors", enumerate_short_vectors__python);
  }
}

using namespace std;
namespace python = boost::python;

// C++ declaration
// void
// enumerate_short_vectors( const std::vector<std::vector<int>>&, unsigned int, unsigned int, std::map<unsigned int, std::vector<std::vector<int>>>& );

// Cython declaration
// cpdef object enumerate_short_vectors__python( object lattice, lower_bound, upper_bound, up_to_sign = False ) :
