"""
CHomP C-Library Interface


EXAMPLES::

    sage: from sage.libs.chomp.matrix_complex import MatrixComplex
    sage: circle = [[0], [0, 0, 1], [1, 0, -1], [0, 1, 1], [1, 1, -1]]
    sage: circle = MatrixComplex(circle)
    sage: circle
    bd(index, dim) = boundary chain
    bd(0, 0) = 0
    bd(1, 0) = 0
    bd(0, 1) = 1[0] + 2[1]
    bd(1, 1) = 1[0] + 2[1]
"""

###############################################################################
#       Copyright (C) 2013, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 3 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################

include 'sage/ext/interrupt.pxi'

from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport int64_t

from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer



                        
cdef extern from "sage_matrix_complex.h":
 
    cdef cppclass SageMatrixComplex:
        SageMatrixComplex(vector[vector[int64_t]] data)
        string to_string(int dim)
        int dimension()
        int size()
        int size(int dim)




cdef class MatrixComplex(SageObject):
    
    cdef SageMatrixComplex* thisptr

    def __cinit__(self, data):
        """
        Construct a chain complex from boundary matrix data

        INPUT:

        - ``data`` -- list of list of integers. 

        EXAMPLES::

            sage: from sage.libs.chomp.matrix_complex import MatrixComplex
            sage: cplx = MatrixComplex([[0], [0, 0, 1], [1, 0, -1], [0, 1, 1], [1, 1, -1]])
            sage: type(cplx)
            <type 'sage.libs.chomp.matrix_complex.MatrixComplex'>
        """
        cdef vector[vector[int64_t]] data_c
        cdef vector[int64_t] row_c
        for row in data:
            row_c.clear()
            for item in row:
                row_c.push_back(item)
            data_c.push_back(row_c)
        self.thisptr = new SageMatrixComplex(data_c)

    cpdef dim(self):
        """
        Return the dimension

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.chomp.examples import matrix_complex_circle
            sage: circle = matrix_complex_circle()
            sage: circle.dim()
            1
        """
        return self.thisptr.dimension()

    cpdef size(self, dim):
        """
        Return the number of cells in the given dimension.

        INPUT:

        - ``dim`` -- integer. The dimension.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.chomp.examples import matrix_complex_circle
            sage: circle = matrix_complex_circle()
            sage: circle.size(0)
            2
            sage: circle.size(1)
            2
            sage: circle.size(2)
            0
        """
        return self.thisptr.size(dim)

    cpdef total_size(self):
        """
        Return the total number of cells.

        INPUT:

        - ``dim`` -- integer. The dimension.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.libs.chomp.examples import matrix_complex_circle
            sage: circle = matrix_complex_circle()
            sage: circle.total_size()
            4
        """
        return self.thisptr.size()

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.
        
        EXAMPLES::

            sage: from sage.libs.chomp.examples import matrix_complex_circle
            sage: matrix_complex_circle()._repr_()
            'bd(index, dim) = boundary chain\nbd(0, 0) = 0\nbd(1, 0) = 0\nbd(0, 1) = 1[0] + 2[1]\nbd(1, 1) = 1[0] + 2[1]'
        """
        result = 'bd(index, dim) = boundary chain\n'
        for dim in range(self.dim()+1):
            result += self.thisptr.to_string(dim).c_str()
        return result.rstrip()
