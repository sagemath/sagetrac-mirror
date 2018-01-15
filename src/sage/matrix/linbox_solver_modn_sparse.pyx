r"""
Solver for sparse matrices over `\ZZ/n\ZZ` using LinBox.
"""
#*****************************************************************************
#       Copyright (C) 2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream& write(const char*, int) except +

    ostream cout

from sage.libs.linbox.conversion cimport (
    set_linbox_matrix_modn_sparse,
    set_linbox_vector_modn_dense,
    new_sage_vector_modn_dense)

from sage.libs.linbox.linbox cimport solve, Method

from .matrix2 cimport Matrix
from .matrix_modn_sparse cimport Matrix_modn_sparse

from sage.modules.vector_modn_sparse cimport c_vector_modint, set_entry

cdef int METHOD_DEFAULT = 0
cdef int METHOD_BLAS_ELIMINATION = 1
cdef int METHOD_SPARSE_ELIMINATION = 2
cdef int METHOD_BLACKBOX = 3
cdef int METHOD_WIEDEMANN = 4

cdef int get_method(str algo):
    if algo is None or algo == "default":
        return METHOD_DEFAULT
    elif algo == "blas_elimination" or \
         algo == "linbox_blas_elimination" or \
         algo == "LinBox::BlasElimination":
        return METHOD_BLAS_ELIMINATION
    elif algo == "sparse_elimination" or \
         algo == "linbox_sparse_elimination" or \
         algo == "LinBox::SparseElimination":
        return METHOD_SPARSE_ELIMINATION
    elif algo == "blackbox" or \
         algo == "linbox_blackbox" or \
         algo == "LinBox::Blackbox":
        return METHOD_BLACKBOX
    elif algo == 'wiedemann' or \
         algo == "linbox_wiedemann" or \
         algo == "LinBox::Wiedeman":
        return METHOD_WIEDEMANN
    else:
        raise ValueError("unknown algo")


# TODO: Cython bug
# https://groups.google.com/forum/#!topic/cython-users/FYIh0OZ1qeE
def __dummy():
    raise RuntimeError()

cdef class LinBoxSolver_modn_sparse(object):
    r"""
    Class to solve linear equations in Sage using LinBox.

    This class handle all allocation/deallocation as well as setting
    properly the underlying LinBox objects. The whole interface is
    accessed with Sage vectors and matrices.

    The two main methods are :meth:`LinBoxSolver_modn_sparse.solve_vector`
    and :meth:`LinBoxSolver_modn_sparse.solve_matrix`.

    EXAMPLES::

        sage: from sage.matrix.linbox_solver_modn_sparse import LinBoxSolver_modn_sparse
        sage: m = matrix(GF(7), [[1,3],[2,4]], sparse=True)
        sage: L = LinBoxSolver_modn_sparse(m)
        sage: v = vector(GF(3), [1,1])
        sage: L.solve_vector(v)
    """
    def __cinit__(self, Matrix_modn_sparse mat):
        self.F = new Modular_uint64(mat.base_ring().characteristic())
        self.A = new SparseMatrix_Modular_uint64(self.F[0], mat._nrows, mat._ncols)
        self.res = new DenseVector_Modular_uint64(self.F[0], mat._ncols)
        self.b = new DenseVector_Modular_uint64(self.F[0], mat._nrows)

    def __dealloc__(self):
        del self.F
        del self.A
        del self.b
        del self.res

    def __init__(self, Matrix_modn_sparse mat):
        r"""
        INPUT:

        - ``mat`` - a sparse matrix over a prime finite field

        EXAMPLES::

            sage: from sage.matrix.linbox_solver_modn_sparse import LinBoxSolver_modn_sparse
            sage: m = matrix(GF(3), [[1,2],[1,0]], sparse=True)
            sage: L = LinBoxSolver_modn_sparse(m)

            sage: v = vector(GF(3), [0,1])
            sage: L.solve_vector(v)
            (1, 1)
        """
        if mat._nrows == 0 or mat._ncols == 0:
            raise ValueError("not implemented for nrows=0 or ncols=0")
        if mat.base_ring().characteristic() == 2:
            raise NotImplementedError("not implemented in characteristic 2")

        # LinBox "solve" is mostly broken for singular matrices. The
        # conditions below could be removed once all LinBox issues
        # have been solved.
        if mat._nrows != mat._ncols or mat.rank() != mat._nrows:
            raise ValueError("only available for full rank square matrices")

        set_linbox_matrix_modn_sparse(self.A[0], mat)
        self.Vin = mat._column_ambient_module().dense_module()
        self.Vout = mat._row_ambient_module().dense_module()

    def solve_vector(self, v, algorithm=None):
        r"""
        Solve the equation ``A x = v`` where ``A`` is the matrix stored in this solver.

        INPUT:

        - ``v`` - vector

        - ``algorithm`` - a string to specify an algorithm to be used by LinBox

        EXAMPLES::

            sage: from sage.matrix.linbox_solver_modn_sparse import LinBoxSolver_modn_sparse
            sage: m = matrix(GF(3), [[1,2],[1,0]], sparse=True)
            sage: L = LinBoxSolver_modn_sparse(m)

            sage: v = vector(GF(3), [0,1])
            sage: u = L.solve_vector(v)
            sage: u
            (1, 1)
            sage: m * u == v
            True

            sage: v = vector(GF(3), [1,2])
            sage: u = L.solve_vector(v)
            sage: u
            (2, 1)
            sage: m * u == v
            True

        You can specify the algorithm to be used by LinBox::

            sage: v = vector(GF(3), [1,1])
            sage: L.solve_vector(v, algorithm='default')
            (1, 0)
            sage: L.solve_vector(v, algorithm='blas_elimination')
            (1, 0)
            sage: L.solve_vector(v, algorithm='sparse_elimination')
            (1, 0)
            sage: L.solve_vector(v, algorithm='blackbox')
            (1, 0)
            sage: L.solve_vector(v, algorithm='wiedemann')
            (1, 0)

        TESTS::

            sage: from sage.matrix.linbox_solver_modn_sparse import LinBoxSolver_modn_sparse
            sage: algos = ["default", "blas_elimination", "sparse_elimination",
            ....:          "blackbox", "wiedemann"]

            sage: m = matrix(GF(3), [[1,1],[0,1]], sparse=True)
            sage: L = LinBoxSolver_modn_sparse(m)
            sage: v = vector(GF(3), [1,1])
            sage: for algo in algos:
            ....:     u = L.solve_vector(v, algorithm=algo)
            ....:     assert m * u == v

        Random testing::

            sage: for i in range(200):
            ....:     dim = randint(1, 30)
            ....:     p = random_prime(10000)
            ....:     if p == 2: p = 3
            ....:     M = MatrixSpace(GF(p), dim, sparse=True)
            ....:     for density in [0.1, 0.25, 0.5]:
            ....:         m = M.random_element(density=density)
            ....:         while m.rank() != dim:
            ....:             m = M.random_element(density=density)
            ....:         L = LinBoxSolver_modn_sparse(m)
            ....:         U = m.column_space()
            ....:         v = U.random_element()
            ....:         w = U.complement().random_element()
            ....:         for algo in algos:
            ....:             u = L.solve_vector(U.zero(), algorithm=algo)
            ....:             assert u.is_zero()
            ....:             u = L.solve_vector(v, algorithm=algo)
            ....:             assert m * u == v
        """
        set_linbox_vector_modn_dense(self.b[0], self.Vin(v))

        cdef int algo = get_method(algorithm)

        if algo == METHOD_DEFAULT:
            solve(self.res[0], self.A[0], self.b[0])
        elif algo == METHOD_BLAS_ELIMINATION:
            solve(self.res[0], self.A[0], self.b[0], Method.BlasElimination())
        elif algo == METHOD_SPARSE_ELIMINATION:
            solve(self.res[0], self.A[0], self.b[0], Method.SparseElimination())
        elif algo == METHOD_BLACKBOX:
            solve(self.res[0], self.A[0], self.b[0], Method.Blackbox())
        elif algo == METHOD_WIEDEMANN:
            solve(self.res[0], self.A[0], self.b[0], Method.Wiedemann())

        return new_sage_vector_modn_dense(self.Vout, self.res[0])

    def solve_matrix(self, mat, algorithm=None):
        r"""
        Solve the equation ``A x = mat`` where ``A`` is the matrix stored in this solver.

        EXAMPLES::

            sage: from sage.matrix.linbox_solver_modn_sparse import LinBoxSolver_modn_sparse
            sage: m = matrix(GF(3), [[1,2],[1,0]], sparse=True)
            sage: b = matrix(GF(3), 2, 4, [1,0,2,0,1,1,2,0])
            sage: L = LinBoxSolver_modn_sparse(m)
            sage: u = L.solve_matrix(b)
            sage: u
            [1 1 2 0]
            [0 1 0 0]
            sage: m * u == b
            True

        TESTS::

            sage: algos = ["default", "blas_elimination", "sparse_elimination",
            ....:          "blackbox", "wiedemann"]

            sage: for _ in range(100):
            ....:     p = random_prime(10000)
            ....:     if p == 2: p = 3
            ....:     dim = randint(1, 100)
            ....:     m = random_matrix(GF(p), dim)
            ....:     while m.rank() != dim:
            ....:         m = random_matrix(GF(p), dim)
            ....:     L = LinBoxSolver_modn_sparse(m)
            ....:     b = random_matrix(GF(3), 30, 10)
            ....:     for algo in algos:
            ....:         u = L.solve_matrix(b, algo)
            ....:         assert m * u == b
        """
        # TODO: for now we just copied the code that used to be in
        # Matrix_modn_sparse. We might want to use dense matrices
        # and call more direct function from LinBox.
        cdef Matrix_modn_sparse B
        cdef Matrix_modn_sparse X
        cdef c_vector_modint * row
        R = self.Vin.base_ring()

        from .constructor import matrix

        if not isinstance(mat, Matrix):
            B = <Matrix_modn_sparse?> matrix(R, mat, sparse=True)
        else:
            B = <Matrix_modn_sparse?> mat.change_ring(R).sparse_matrix()

        cdef int algo = get_method(algorithm)

        X = matrix(R, B.ncols(), self.A.coldim(), sparse=True)
        B = B.transpose()
        cdef size_t i, j
        for i in range(X.nrows()):
            # set self.b to the i-th row of B
            row = B.rows + i
            # TODO: no method in LinBox to set a vector to zero?
            for j in range(self.A.coldim()):
                self.b.setEntry(j, 0)
            for j in range(row.num_nonzero):
                self.b.setEntry(row.positions[j], row.entries[j])

            # solve the current row
            if algo == METHOD_DEFAULT:
                solve(self.res[0], self.A[0], self.b[0])
            elif algo == METHOD_BLAS_ELIMINATION:
                solve(self.res[0], self.A[0], self.b[0], Method.BlasElimination())
            elif algo == METHOD_SPARSE_ELIMINATION:
                solve(self.res[0], self.A[0], self.b[0], Method.SparseElimination())
            elif algo == METHOD_BLACKBOX:
                solve(self.res[0], self.A[0], self.b[0], Method.Blackbox())
            elif algo == METHOD_WIEDEMANN:
                solve(self.res[0], self.A[0], self.b[0], Method.Wiedemann())

            # set i-th row of X to be self.res
            for j in range(self.A.coldim()):
                set_entry(X.rows + i, j, self.res[0].getEntry(j))

        return X.transpose()
