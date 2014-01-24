#clang C++

#===============================================================================
# 
# Copyright (C) 2014 Martin Raum
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

include "interrupt.pxi"

from enumerate_short_vectors_python cimport short_vectors as svs

cpdef object short_vectors( object lattice, lower_bound, upper_bound, up_to_sign = False ) :
    r"""
    Enumerate vectors of minimal norm ``lower_bound`` and maximal norm ``upper_bound``, either up to sign or not.

    INPUT:
    
    - ``lattice`` - a list of list of integers which corresponds to the Gram matrix of a binary quadratic form.

    - ``lower_bound`` - a positive integer.

    - ``upper_bound`` - a positive integer.

    - ``up_to_sign`` - a boolean (default: ``False``).

    OUPUT:

    - A dictionary mapping integers to a list of tuples.  Each tuple corresponds to a vector.

    """
    sig_on()
    result = svs( lattice, lower_bound, upper_bound, up_to_sign )
    sig_off()
    return result
