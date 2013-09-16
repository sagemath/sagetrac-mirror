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
    bd(0, 1) = 1[0] + -1[1]
    bd(1, 1) = 1[0] + -1[1]
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
from libcpp.pair cimport pair
from libcpp.string cimport string
from libc.stdint cimport int64_t

from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer



cdef extern from "sage_matrix_complex.h":
 
    cdef cppclass Complex:
        int dimension()

    cdef cppclass SageMatrixComplex(Complex):
        SageMatrixComplex(vector[vector[int64_t]] data)
        string to_string(int dim)
        bint is_chain_complex()
        int size()
        int size(int dim)



cdef extern from "chomp/Ring.h" namespace "chomp":

    cdef cppclass Long:
        Long(int64_t)
        bint equals "operator ==" (Long)

    ctypedef Long Ring



cdef extern from "chomp/Chain.h" namespace "chomp":

    cdef cppclass Chain:
        int dimension()
        



cdef extern from "chomp/Generators.h" namespace "chomp":

    cdef cppclass Generators_t(vector[vector[pair[Chain, Ring]]]):
        Generators_t(Generators_t src)

    Generators_t SmithGenerators(Complex)
    Generators_t MorseGenerators(Complex)



cdef class Generators(SageObject):

    cdef Generators_t* thisptr

    def __cinit__(self):
        self.thisptr = NULL

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        assert self.thisptr!=NULL, 'Do not construct Generators manually, only via make_Generators!'
        del self.thisptr

    def betti(self, dim):
        """
        Return the Betti number

        INPUT:

        - ``dim`` -- integer. The dimension (degree).

        OUTPUT:

        An integer. The Betti number in the given dimension.

        EXAMPLES::

            sage: from sage.libs.chomp.examples import matrix_complex_segment
            sage: segment = matrix_complex_segment()
            sage: gens = segment.Morse_generators()
            sage: [gens.betti(d) for d in range(0,3)]

            sage: from sage.libs.chomp.examples import matrix_complex_circle
            sage: circle = matrix_complex_circle()
            sage: gens = circle.Smith_generators()
            sage: [gens.betti(d) for d in range(0,3)]

            sage: gens = circle.Morse_generators()
            sage: [gens.betti(d) for d in range(0,3)]
        """
        print 'betti', self.thisptr[0].size()
        if dim < 0 or dim >= self.thisptr[0].size():
            return 0
        cdef int betti = 0
        cdef int i
        for i in range(self.thisptr[0][dim].size()):
            print dim, i, self.thisptr[0][dim][i].first.dimension()
            if self.thisptr[0][dim][i].second.equals(Ring(0)):
                betti += 1
        return betti


    
cdef make_Generators(Generators_t gens_c):
    """
    Construct a :class:`Generators` Cython class

    INPUT:

    - ``gens_c`` -- A ``chomp::Generators_t`` C++ object.
    """
    cdef Generators gens = Generators()
    gens.thisptr = new Generators_t(gens_c)
    return gens



cdef class MatrixComplex(SageObject):
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
    
    cdef SageMatrixComplex* thisptr

    def __cinit__(self, data):
        """
        the Cython constructor
        """
        cdef vector[vector[int64_t]] data_c
        cdef vector[int64_t] row_c
        for row in data:
            row_c.clear()
            for item in row:
                row_c.push_back(item)
            data_c.push_back(row_c)
        self.thisptr = new SageMatrixComplex(data_c)

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

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

    def is_chain_complex(self):
        """
        Verify the chain complex condition

        OUTPUT:

        Boolean. Whether `d\circ d=0` for all differentials.

        EXAMPLES::
        
            sage: from sage.libs.chomp.examples import matrix_complex_circle
            sage: circle = matrix_complex_circle()
            sage: circle.is_chain_complex()
            True
        """
        return self.thisptr.is_chain_complex()

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.
        
        EXAMPLES::

            sage: from sage.libs.chomp.examples import matrix_complex_circle
            sage: matrix_complex_circle()._repr_()
            'bd(index, dim) = boundary chain\nbd(0, 0) = 0\nbd(1, 0) = 0\nbd(0, 1) = 1[0] + -1[1]\nbd(1, 1) = 1[0] + -1[1]'
        """
        result = 'bd(index, dim) = boundary chain\n'
        for dim in range(self.dim()+1):
            result += self.thisptr.to_string(dim).c_str()
        return result.rstrip()

    def Smith_generators(self):
        """
        Return the generators computed via the Smith normal form

        EXAMPLES::

            sage: from sage.libs.chomp.examples import matrix_complex_circle
            sage: circle = matrix_complex_circle()
            sage: circle.Smith_generators()
            <type 'sage.libs.chomp.matrix_complex.Generators'>
        """
        return make_Generators(SmithGenerators(self.thisptr[0]))

    def Morse_generators(self):
        """
        Return the generators computed via the discrete Morse complex

        EXAMPLES::

            sage: from sage.libs.chomp.examples import matrix_complex_circle
            sage: circle = matrix_complex_circle()
            sage: circle.Morse_generators()
            <type 'sage.libs.chomp.matrix_complex.Generators'>
        """
        return make_Generators(MorseGenerators(self.thisptr[0]))
