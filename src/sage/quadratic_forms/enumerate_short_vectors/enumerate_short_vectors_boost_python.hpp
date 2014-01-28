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

#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/tuple.hpp>


boost::python::dict
enumerate_short_vectors_boost_python( boost::python::list, long, long, bool );

boost::python::tuple
to_python_tuple( const std::vector<int> & );

boost::python::tuple
to_python_neg_tuple( const std::vector<int> & );
