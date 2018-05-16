#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

##############################################################################
#
# ntl_mat_ZZ_p: Matrices over ZZ_p via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version of mat_GF2E (based on code by William Stein)
#  - Alex J. Best <alex.j.best@gmail.com>
#    2018-02: transplanted to mat_ZZ_p
#
##############################################################################
from __future__ import absolute_import

from cysignals.signals cimport sig_on, sig_off

include 'misc.pxi'
include 'decl.pxi'

from cpython.object cimport Py_EQ, Py_NE
from .ntl_ZZ_p cimport ntl_ZZ_p
from .ntl_ZZ_pContext import ntl_ZZ_pContext
from .ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.rings.integer cimport Integer
from sage.misc.randstate cimport randstate, current_randstate
from cpython.object cimport PyObject_RichCompare
from sage.ext.cplusplus cimport ccrepr

from sage.libs.ntl.ntl_ZZ import unpickle_class_args


cdef class ntl_mat_ZZ_p(object):
    r"""
    The \class{mat_ZZ_p} class implements arithmetic with matrices over $\ZZ/p\ZZ$.
    """
    def __init__(self, modulus = None, nrows=0, ncols=0, v=None):
        """
        Constructs a matrix over ntl.ZZ_p.

        INPUT:
            modulus -- ZZ_p context
            nrows -- number of rows
            ncols -- number of columns
            v     -- either a list or a matrix over Z/pZ

        EXAMPLES::

            sage: k = Integers(20)
            sage: ctx = ntl.ZZ_pContext(k.characteristic())
            sage: ntl.mat_ZZ_p(ctx, 5,5, [0..24])
            [[0 1 2 3 4]
            [5 6 7 8 9]
            [10 11 12 13 14]
            [15 16 17 18 19]
            [0 1 2 3 4]
            ]
            sage: ntl.mat_ZZ_p(ctx, 5,5)
            [[0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            ]
            sage: A= matrix(k,5,5,[k(_) for _ in range(25)])
            sage: ntl.mat_ZZ_p(ctx, A)
            [[0 1 2 3 4]
            [5 6 7 8 9]
            [10 11 12 13 14]
            [15 16 17 18 19]
            [0 1 2 3 4]
            ]
        """
        if modulus is None:
            raise ValueError("You must specify a modulus when creating a ZZ_p.")

        cdef unsigned long _nrows, _ncols
        cdef unsigned long i, j

        from sage.structure.element import is_Matrix
        if is_Matrix(nrows):
            _nrows = nrows.nrows()
            _ncols = nrows.ncols()
            v     = nrows.list()
        else:
            _nrows = nrows
            _ncols = ncols

        self.x.SetDims(_nrows, _ncols)
        if v is not None:
            sig_on()
            for i from 0 <= i < _nrows:
                for j from 0 <= j < _ncols:
                    elem = v[i*_ncols+j]
                    if not isinstance(elem, ntl_ZZ_p):
                        elem = ntl_ZZ_p(elem, modulus)
                    mat_ZZ_p_setitem(&self.x, i, j, &(<ntl_ZZ_p>elem).x)
            sig_off()

    def __cinit__(self, modulus=None, nrows=0, ncols=0, v=None):
        #################### WARNING ###################
        ## Before creating a ZZ_p, you must create a  ##
        ## ZZ_pContext, and restore it.  In Python,   ##
        ## the error checking in __init__ will prevent##
        ## you from constructing an ntl_ZZ_p          ##
        ## inappropriately.  However, from Cython, you##
        ## could do r = ntl_ZZ_p.__new__(ntl_ZZ_p) without
        ## first restoring a ZZ_pContext, which could ##
        ## have unfortunate consequences.  See _new  ##
        ## defined below for an example of the right  ##
        ## way to short-circuit __init__ (or just call##
        ## _new in your own code).                    ##
        ################################################
        if modulus is None:
            return
        if isinstance(modulus, ntl_ZZ_pContext_class):
            self.c = <ntl_ZZ_pContext_class>modulus
            self.c.restore_c()
        else:
            self.c = <ntl_ZZ_pContext_class>ntl_ZZ_pContext(modulus)
            self.c.restore_c()

    cdef ntl_ZZ_p _new_element(self):
        cdef ntl_ZZ_p r
        self.c.restore_c()
        r = ntl_ZZ_p.__new__(ntl_ZZ_p)
        r.c = self.c
        return r

    cdef ntl_mat_ZZ_p _new(self):
        cdef ntl_mat_ZZ_p r
        self.c.restore_c()
        r = ntl_mat_ZZ_p.__new__(ntl_mat_ZZ_p)
        r.x.SetDims(self.x.NumRows(),self.x.NumCols())
        r.c = self.c
        return r

    def modulus_context(self):
        """
        Returns the structure that holds the underlying NTL ZZ_p modulus.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(135)
            sage: a = ntl.ZZ_p(200, ctx)
            sage: A= ntl.mat_ZZ_p(ctx, 1, 1, [a])
            sage: cty = A.modulus_context(); cty
            NTL modulus 135
            sage: ctx == cty
            True
        """
        return self.c

    def __reduce__(self):
        """
        EXAMPLES::

            sage: k = Integers(20)
            sage: ctx = ntl.ZZ_pContext(k.characteristic())
            sage: A = ntl.mat_ZZ_p(ctx, 5,5, [0..24])
            sage: A == loads(dumps(A))
            True
        """
        return unpickle_class_args, (ntl_mat_ZZ_p, (self.modulus_context(), self.x.NumRows(), self.x.NumCols(), self.list()))

    def __repr__(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(200)
            sage: ntl.mat_ZZ_p(ctx, 2,2,range(4)).__repr__()
            '[[0 1]\n[2 3]\n]'
        """
        self.c.restore_c()
        return ccrepr(self.x)

    def __mul__(ntl_mat_ZZ_p self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: n = ntl.mat_ZZ_p(ctx, 5,5,[3..27])
            sage: m*n ## indirect doctest
            [[0 10 0 10 0]
            [5 0 15 10 5]
            [10 10 10 10 10]
            [15 0 5 10 15]
            [0 10 0 10 0]
            ]
        """
        cdef ntl_mat_ZZ_p r = self._new()
        if not isinstance(other, ntl_mat_ZZ_p):
            other = ntl_mat_ZZ_p(other, self.c)
        if not self.c is (<ntl_mat_ZZ_p>other).c:
            raise ValueError("You can not perform arithmetic with matrices over different rings.")
        sig_on()
        mat_ZZ_p_mul(r.x, self.x, (<ntl_mat_ZZ_p>other).x)
        sig_off()
        return r

    def __sub__(ntl_mat_ZZ_p self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: n = ntl.mat_ZZ_p(ctx, 5,5,[3..27])
            sage: m-n ## indirect doctest
            [[17 17 17 17 17]
            [17 17 17 17 17]
            [17 17 17 17 17]
            [17 17 17 17 17]
            [17 17 17 17 17]
            ]
        """
        cdef ntl_mat_ZZ_p r = self._new()
        if not isinstance(other, ntl_mat_ZZ_p):
            other = ntl_mat_ZZ_p(other, self.c)
        if not self.c is (<ntl_mat_ZZ_p>other).c:
            raise ValueError("You can not perform arithmetic with matrices over different rings.")
        sig_on()
        mat_ZZ_p_sub(r.x, self.x, (<ntl_mat_ZZ_p>other).x)
        sig_off()
        return r

    def __add__(ntl_mat_ZZ_p self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: n = ntl.mat_ZZ_p(ctx, 5,5,[3..27])
            sage: m+n ## indirect doctest
            [[3 5 7 9 11]
            [13 15 17 19 1]
            [3 5 7 9 11]
            [13 15 17 19 1]
            [3 5 7 9 11]
            ]
        """
        cdef ntl_mat_ZZ_p r = self._new()
        if not isinstance(other, ntl_mat_ZZ_p):
            other = ntl_mat_ZZ_p(other, self.c)
        if not self.c is (<ntl_mat_ZZ_p>other).c:
            raise ValueError("You can not perform arithmetic with matrices over different rings.")
        sig_on()
        mat_ZZ_p_add(r.x, self.x, (<ntl_mat_ZZ_p>other).x)
        sig_off()
        return r

    def __neg__(ntl_mat_ZZ_p self):
        """
        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: n = -m; n ## indirect doctest
            [[0 19 18 17 16]
            [15 14 13 12 11]
            [10 9 8 7 6]
            [5 4 3 2 1]
            [0 19 18 17 16]
            ]
            sage: m + n
            [[0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            ]
        """
        cdef ntl_mat_ZZ_p r = self._new()
        sig_on()
        mat_ZZ_p_negate(r.x, self.x)
        sig_off()
        return r

    def __pow__(ntl_mat_ZZ_p self, long e, ignored):
        """
        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: m**2 == m*m ## indirect doctest
            True
        """
        cdef ntl_mat_ZZ_p r = self._new()
        sig_on()
        mat_ZZ_p_power(r.x, self.x, e)
        sig_off()
        return r

    def __richcmp__(ntl_mat_ZZ_p self, other, int op):
        """
        Compare self to other.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: n = ntl.mat_ZZ_p(ctx, 5,5,[3..27])
            sage: m == n
            False
            sage: m == m
            True
            sage: m == []
            False
        """
        self.c.restore_c()

        cdef ntl_mat_ZZ_p b
        try:
            b = <ntl_mat_ZZ_p?>other
        except TypeError:
            return NotImplemented
        return PyObject_RichCompare(self.list(), other.list(), op)

    def NumRows(self):
        """
        Return the number of rows in self.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24]) ; m.NumRows()
            5
        """
        return int(self.x.NumRows())

    def NumCols(self):
        """
        Return the number of columns in self.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24]) ; m.NumCols()
            5
        """
        return int(self.x.NumCols())

    def __setitem__(self, ij, x):
        """
        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: m[0,0]
            0
            sage: m[0,0] = 1
            sage: m[0,0]
            1
        """
        cdef int i, j
        if not isinstance(x, ntl_ZZ_p):
            x = ntl_ZZ_p(x, self.c)

        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.x.NumCols()==1 and (isinstance(ij, Integer) or isinstance(ij, int)):
            i = ij
            j = 0
        elif self.x.NumRows()==1 and (isinstance(ij, Integer) or isinstance(ij, int)):
            i = 0
            j = ij
        else:
            raise TypeError('ij is not a matrix index')

        if i < 0 or i >= self.x.NumRows() or j < 0 or j >= self.x.NumCols():
            raise IndexError("array index out of range")

        if not (<ntl_ZZ_p>x).c is self.c:
            raise ValueError("You can not assign elements from different rings.")

        self.c.restore_c()

        mat_ZZ_p_setitem(&self.x, i, j, &(<ntl_ZZ_p>x).x)

    def __getitem__(self, ij):
        """
        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: m[0,1]
            1
            sage: m[0,0] = 0
            sage: m[0,0]
            0
        """
        cdef int i, j
        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.x.NumCols() == 1 and (isinstance(ij, Integer) or isinstance(ij, int)):
            i = ij
            j = 0
        elif self.x.NumRows() == 1 and (isinstance(ij, Integer) or isinstance(ij, int)):
            i = 0
            j = ij
        else:
            raise TypeError('ij is not a matrix index')

        if i < 0 or i >= self.x.NumRows() or j < 0 or j >= self.x.NumCols():
            raise IndexError("array index out of range")

        cdef ntl_ZZ_p e = self._new_element()
        e.x = self.x.get( i+1, j+1 )
        return e

    def determinant(self):
        """
        Returns the determinant.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(37)
            sage: m = ntl.mat_ZZ_p(ctx, 3,3,[0..8])
            sage: m[0,0] = 1
            sage: m.determinant()
            34
        """
        cdef ntl_ZZ_p r = self._new_element()
        sig_on()
        r.x = mat_ZZ_p_determinant(self.x)
        sig_off()
        return r

    def gauss(self,ncols=-1):
        """
        Performs unitary row operations so as to bring this matrix
        into row echelon form.  If the optional argument \code{ncols}
        is supplied, stops when first ncols columns are in echelon
        form.  The return value is the rank (or the rank of the first
        ncols columns).

        INPUT:

        - ``ncols`` - number of columns to process (default: all)

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(37)
            sage: n = ntl.mat_ZZ_p(ctx, 5,5,[3..27])
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[3..27])
            sage: m.gauss()
            2
            sage: n == m
            False
            sage: ntl.mat_ZZ_p(ctx, 5,5).gauss()
            0
        """
        if ncols == -1:
            ncols = self.x.NumCols()
        return int(mat_ZZ_p_gauss(self.x, int(ncols)))

    def list(self):
        """
        Returns a list of the entries in this matrix

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(10)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[ntl.ZZ_p_random(10) for x in range(5*5)])
            sage: m.list()
            [3, 7, 9, 7, 4, 1, 3, 0, 5, 0, 7, 0, 2, 3, 8, 3, 0, 4, 7, 7, 4, 9, 2, 5, 7]
        """
        return [self[i,j] for i in range(self.NumRows()) for j in range(self.x.NumCols())]

    def IsZero(self):
        """
        Return True if self is zero, and false otherwise.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: n = ntl.mat_ZZ_p(ctx, 5,5)
            sage: m.IsZero()
            False
            sage: n.IsZero()
            True
        """
        cdef long isZero
        sig_on()
        isZero = mat_ZZ_p_IsZero(self.x)
        sig_off()
        return bool(isZero)

    def _sage_(ntl_mat_ZZ_p self, k=None):
        """
        Returns a ``Matrix`` over a ``IntegerMod`` representation
        of this element.

        INPUT:

        - ``k`` - optional IntegerModRing(N)

        OUTPUT:
        - Sage matrix over ``\ZZ/p\ZZ`` or ``k``, if supplied

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 2,3,[9..16])
            sage: m
            [[9 10 11]
            [12 13 14]
            ]
            sage: m._sage_()
            [ 9 10 11]
            [12 13 14]
            sage: m._sage_(Integers(10))
            [9 0 1]
            [2 3 4]
        """
        cdef ZZ_c p
        if k is None:
            from sage.rings.finite_rings.integer_mod_ring import Integers
            k = Integers(self.c.modulus())

        l = [k(e._sage_()) for e in self.list()] # we actually can do faster than this

        from sage.matrix.constructor import Matrix
        return Matrix(k,self.x.NumRows(),self.x.NumCols(),l)

    def transpose(ntl_mat_ZZ_p self):
        """
        Returns the transposed matrix of self.

        OUTPUT:
            transposed Matrix

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 5,5,[0..24])
            sage: n = m.transpose()
            sage: n == m
            False
            sage: n.transpose() == m
            True
        """
        cdef ntl_mat_ZZ_p r = self._new()
        sig_on()
        mat_ZZ_p_transpose(r.x, self.x)
        sig_off()
        return r

    def __invert__(self):
        """
        Return $X = A^{-1}$; an error is raised if A is singular.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(5)
            sage: m = ntl.mat_ZZ_p(ctx, 2,2,[0..3])
            sage: n = ~m
            sage: o = n*m
            sage: o.IsIdent()
            True
        """
        cdef ntl_mat_ZZ_p r = self._new()
        sig_on()
        mat_ZZ_p_inv(r.x, self.x)
        sig_off()
        return r

    def IsIdent(self, n = -1):
        """
        test if A is the n x n identity matrix

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(5)
            sage: m = ntl.mat_ZZ_p(ctx, 2,2,[0..3])
            sage: n = ~m
            sage: o = n*m
            sage: o.IsIdent()
            True
        """
        if n < 0:
            n = self.NumRows()
        return bool(mat_ZZ_p_IsIdent(self.x, n))

    def IsDiag(self, long n, ntl_ZZ_p d):
        """
        Test if X is an  n x n diagonal matrix with d on diagonal.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(20)
            sage: m = ntl.mat_ZZ_p(ctx, 3,3,[4,0,0,0,4,0,0,0,4])
            sage: m.IsDiag(2, ntl.ZZ_p(4,ctx))
            False
            sage: m.IsDiag(3, ntl.ZZ_p(3,ctx))
            False
            sage: m.IsDiag(3, ntl.ZZ_p(4,ctx))
            True
        """
        return bool(mat_ZZ_p_IsDiag(self.x, n, d.x))

    def image(self):
        """
        The rows of X are computed as basis of A's row space.  X is
        row echelon form.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(37)
            sage: m = ntl.mat_ZZ_p(ctx, 3,3,[0..8])
            sage: m.image()
            [[3 4 5]
            [0 1 2]
            ]
        """
        cdef ntl_mat_ZZ_p X = self._new()
        sig_on()
        mat_ZZ_p_image(X.x, self.x)
        sig_off()
        return X

    def kernel(self):
        """
        Computes a basis for the kernel of the map ``x -> x*A``, where
        ``x`` is a row vector.

        EXAMPLES::

            sage: ctx = ntl.ZZ_pContext(37)
            sage: m = ntl.mat_ZZ_p(ctx, 3,3,[0..8])
            sage: m.kernel()
            [[1 35 1]
            ]
        """
        cdef ntl_mat_ZZ_p X = self._new()
        sig_on()
        mat_ZZ_p_kernel(X.x, self.x)
        sig_off()
        return X

    def randomize(self, density=1, nonzero=False):
        """
        Randomize ``density`` proportion of the entries of this matrix,
        leaving the rest unchanged.

        INPUT:

        -  ``density`` - float; proportion (roughly) to be considered for
           changes
        -  ``nonzero`` - Bool (default: ``False``); whether the new entries
           are forced to be non-zero

        EXAMPLES::

            sage: k = Integers(20)
            sage: ctx = ntl.ZZ_pContext(k.characteristic())
            sage: A = ntl.mat_ZZ_p(ctx, 100,100)
            sage: A.randomize()
            sage: len([e for e in A.list() if e!=0])
            9486

            sage: A = ntl.mat_ZZ_p(ctx, 100,100)
            sage: A.randomize(nonzero=True)
            sage: len([e for e in A.list() if e!=0])
            10000

            sage: A = ntl.mat_ZZ_p(ctx, 100,100)
            sage: A.randomize(nonzero=True, density=0.1)
            sage: len([e for e in A.list() if e!=0])
            994

        """
        cdef long i,j
        cdef ZZ_p_c tmp

        cdef float _density = density
        cdef randstate rstate = current_randstate()

        if _density <= 0:
            return
        if _density > 1:
            _density = 1.0

        if not nonzero:
            if _density == 1.0:
                for i in xrange(self.x.NumRows()):
                    for j in xrange(self.x.NumCols()):
                        ZZ_p_random(tmp)
                        mat_ZZ_p_setitem(&self.x, i, j, &tmp)
            else:
                for i in xrange(self.x.NumRows()):
                    for j in xrange(self.x.NumCols()):
                        if rstate.c_rand_double() <= _density:
                            ZZ_p_random(tmp)
                            mat_ZZ_p_setitem(&self.x, i, j, &tmp)
        else:
            if _density == 1.0:
                for i in xrange(self.x.NumRows()):
                    for j in xrange(self.x.NumCols()):
                        ZZ_p_random(tmp)
                        while ZZ_p_IsZero(tmp):
                            ZZ_p_random(tmp)
                        mat_ZZ_p_setitem(&self.x, i, j, &tmp)
            else:
                for i in xrange(self.x.NumRows()):
                    for j in xrange(self.x.NumCols()):
                        if rstate.c_rand_double() <= _density:
                            ZZ_p_random(tmp)
                            while ZZ_p_IsZero(tmp):
                                ZZ_p_random(tmp)
                            mat_ZZ_p_setitem(&self.x, i, j, &tmp)
