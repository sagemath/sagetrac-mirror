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

boost::python::dict
short_vectors( boost::python::list, const long, const long, const bool );

BOOST_PYTHON_FUNCTION_OVERLOADS(short_vectors_overloads, short_vectors, 3, 4)

boost::python::tuple
to_python_vector( const std::vector<int> & );

boost::python::tuple
to_python_neg_vector( const std::vector<int> & );
