r"""
Inline conversions between LinBox and Sage

Each LinBox type has a corresponding Sage types and we use the following
conventions for conversion functions

- ``new_linbox_XXX`` : create a new linbox object
- ``new_sage_XXX``   : create a new Sage object
- ``set_linbox_XXX`` : set the entries of the linbox object
- ``set_sage_XXX``   : set the entries of the Sage object

For matrices that uses a flint datastructure, see the lower level conversions
in the module ``linbox_flint_interface``.
"""
#*****************************************************************************
#       Copyright (C) 2007 Martin Albrecht
#       Copyright (C) 2008 Clement Pernet
#       Copyright (C) 2018 Vincent Delecroix
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .givaro cimport Modular_uint64
from .linbox cimport SparseMatrix_Modular_uint64, DenseVector_Modular_uint64

from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse
from sage.modules.vector_modn_dense cimport Vector_modn_dense

from sage.modules.vector_modn_sparse cimport c_vector_modint

################################################
# matrix_modn_sparse (sparse matrix over Z/nZ) #
################################################

# set the entries of a LinBox matrix from a Sage matrix
# INPUT:
#   A - LinBox matrix
#   m - Sage matrix
cdef inline void set_linbox_matrix_modn_sparse(SparseMatrix_Modular_uint64& A, Matrix_modn_sparse m):
    cdef c_vector_modint * row
    cdef size_t i, j
    for i in range(m._nrows):
        row = m.rows + i
        for j in range(row.num_nonzero):
            A.setEntry(i, row.positions[j], row.entries[j])

# return a new LinBox matrix from a Sage matrix
# (such matrix has to be deallocated with a "del" statement)
# INPUT:
#   F - LinBox field
#   m - Sage matrix
cdef inline SparseMatrix_Modular_uint64 * new_linbox_matrix_modn_sparse(Modular_uint64 &F, Matrix_modn_sparse m):
    cdef SparseMatrix_Modular_uint64 * A = new SparseMatrix_Modular_uint64(F, m._nrows, m._ncols)
    set_linbox_matrix_modn_sparse(A[0], m)
    return A

##############################################
# vector_modn_dense (dense vector over Z/nZ) #
##############################################

# set the entries of a LinBox vector from a Sage vector
# INPUT:
#   res - LinBox vector
#   v - Sage vector
cdef inline void set_linbox_vector_modn_dense(DenseVector_Modular_uint64& res, Vector_modn_dense v):
    cdef size_t i
    for i in range(v._degree):
        res.setEntry(i, v._entries[i])

# return a new LinBox vector from a Sage vector
# (such vector has to be deallocated with a "del" statement)
# INPUT:
#   F - LinBox field
#   v - Sage vector
cdef inline DenseVector_Modular_uint64 * new_linbox_vector_modn_dense(Modular_uint64& F, Vector_modn_dense v):
    cdef DenseVector_Modular_uint64 * V = new DenseVector_Modular_uint64(F, v._degree)
    set_linbox_vector_modn_dense(V[0], v)
    return V

# return a new Sage vector from a LinBox one
# INPUT:
#   P - parent for the Sage vector
#   v - LinBox vector
cdef inline Vector_modn_dense new_sage_vector_modn_dense(P, DenseVector_Modular_uint64& v):
    cdef Vector_modn_dense res = <Vector_modn_dense?> P()

    vec = &v.refRep()
    cdef size_t i
    for i in range(res._degree):
        res._entries[i] = vec[0][i]
    return res
