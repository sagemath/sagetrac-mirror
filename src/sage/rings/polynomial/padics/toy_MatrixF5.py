"""
This module implements Groebner bases computation via Matrix-F5 algorithm.

Computations over finite-precision complete discrete valuation field are available with the option "weak" in toy_matrixF5.
Computations of tropical Groebner bases are available via the tropical option of toy_matrixF5.

AUTHOR:

- Tristan Vaccon (2013--2016)
"""

#*****************************************************************************
#  Copyright (C) 2013 Tristan Vaccon <tristanvaccon@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ


from sage.modules.free_module_element import FreeModuleElement_generic_dense
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.matrix.matrix_generic_dense import Matrix_generic_dense
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix, Matrix, zero_matrix, random_matrix
from sage.rings.arith import binomial
from copy import copy

def list_monom(P,d):
    r"""

    Return an ordered list of the monomials of P of degree d
 
    INPUT:

    - ``P`` -- a polynomial ring
    - ``d`` -- an integer

    OUTPUT:

    The list of monomials in P of degree d, ordered from the biggest to the smallest according to the monomial ordering on P

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_MatrixF5 import list_monom
        sage: list_monom(S,2)
        [x^2, x*y, y^2, x*z, y*z, z^2]
    """
    if d == 0:
        return [P(1)]
    l1 = [ a for a in P.gens()]
    if d == 1:
        return l1
    ld1 = list_monom(P,d-1)
    ld = [[x*u for x in l1] for u in ld1]
    ld = sum(ld,[])
    s = []
    for i in ld:
        if i not in s:
            s.append(i)
    s.sort()
    s.reverse()
    return s

def firstvar(mon):
    r"""

    INPUT:

    - ``mon`` -- a monomial

    OUTPUT:

    The index of the first variable appearing in the monomial mon, or -1 for a constant monomial

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_MatrixF5 import firstvar
        sage: firstvar(x*y^2*z)
        0
    """
    if mon.degree() == 0:
        return (-1)
    else:
        P = mon.parent()
        variables = list_monom(P,1)
        variables.reverse()
        n = len(variables)
        l = mon.degrees()
        i = 0
        while l[i] == 0:
            i += 1
        return i





class MacaulayMatrix:
    r"""
    The class of Macaulay matrices.
    """

    def __init__(self, base, monomials, entries, labels=None, signatures=None, s=None, indexes=None):
        r"""

        INPUT:

        - ``base`` -- the polynomial base ring.
        - ``monomials`` -- the list of monomials (in decreasing order) of the Macaulay matrix.
        - ``entries`` -- the coefficients of the matrix.
        - ``labels`` -- the labels of the rows of the matrix.
        - ``signatures`` -- the signatures of the rows of the matrix.
        - ``s`` -- the length of the labels (the number of generating polynomials)
        - ``indexes`` -- the list of the indexes (column of the first non-zero entry) of the rows of the matrix, or the number of columns if the rows has no non-zero entry.

        OUTPUT:

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [1,0,0,0,2,3]), Matrix(T, 2, 2, [2,0,3,4]))
        """
        self._monomials = monomials
        self._ncols = len(monomials)
        self._entries = entries
        self._labels = labels
        self._signatures = signatures
        self._s = s
        self._indexes = indexes
        self._base = base
        self._zero_poly = base.zero()
        self._zero_coeff = base.base_ring().zero()

    def nrows(self):
        r"""

        INPUT:

        OUTPUT:

        The number of rows of the Macaulay matrix self

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [1,0,0,0,2,3]), Matrix(T, 2, 2, [2,0,3,4]))
            sage: M.nrows()
            2
        """
        return self._entries.nrows()

    def ncols(self):
        r"""

        INPUT:

        OUTPUT:

        The number of columns of the Macaulay matrix self

        EXAMPLES::
            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [1,0,0,0,2,3]), Matrix(T, 2, 2, [2,0,3,4]))
            sage: M.ncols()
            3
        """
        return self._ncols

    def labels(self):
        r"""
        INPUT:

        OUTPUT:

        The matrix of labels of this matrix.

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2], Matrix(QQ, 2, 3, [1,0,0,0,2,3]), Matrix(T, 2, 2, [2,0,3,4]))
            sage: M.labels()
            [2 0]
            [3 4]
        """
        return self._labels

    def label(self,i):
        r"""

        INPUT:

        - ``i`` -- an integer

        OUTPUT:

        The label of the i-th line of this matrix.

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [1,0,0,0,2,3]), Matrix(T, 2, 2, [2,0,3,4]))
            sage: M.label(1)
            (3, 4)
        """
        return self._labels.row(i)

    def signature_via_label(self,i):
        """

        INPUT:

        - ``i`` -- an integer

        OUTPUT:

        The signature of the i-th line of this matrix, as read on the label of this line.

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2], Matrix(QQ, 2, 3, [1,0,0,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]) , s = 2)
            sage: M.signature(1)
            [y, 1]
        """
        if self.labels() is None:
            return None
        else:
            s = self._s
            lab = self.label(i)
            i = 0
            while i<s and (lab[s-1-i].is_zero()):
                i += 1
            return [lab[s-1-i].lm(),s-1-i]

    def signature(self,i):
        """
        INPUT:

        - ``i`` -- an integer

        OUTPUT:

        The signature of the i-th line of this matrix

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2], Matrix(QQ, 2, 3, [1,0,0,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]) , signatures = [[x,0],[y,1]], s = 2)
            sage: M.signature(1)
            [y, 1]
        """
        if self._signatures is not None:
            return self._signatures[i]
        if self._labels is not None:
            return self.signature_via_label(i)

    def signature_tropical(self,i,tropical_weight, tropical_order_ring):
        """
        INPUT:

        - ``i`` -- an integer

        OUTPUT:

        The tropical signature of the i-th line of this matrix.

        EXAMPLES::

            sage: T.<x,y,z> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: U.<a,b,c> = PolynomialRing(Qp(7), order='lex')
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2,x*z,y*z,z^2],Matrix(QQ, 2, 6, [1,0,0,0,0,0,0,2,3,0,0,0]), Matrix(T, 2, 2, [x,0,x+y,y]), s=2)
            sage: M.signature_tropical(1,[1,2,4],U)
            [y, 1]
        """
        if self.labels() is None:
            return None
        else:
            s = self._s
            lab = self.label(i)
            i = 0
            while i<s and (lab[s-1-i].is_zero()):
                i += 1
            listmonom = lab[s-1-i].monomials()
            l = len(listmonom)
            nvar = self._base.ngens()
            T = tropical_order_ring.base_ring()
            max_monomial = listmonom[0]
            max_index = 0
            max_trop_weight = sum([listmonom[0].degrees()[v]*tropical_weight[v] for v in range(nvar)])
            for j0 in range(1,l):
                term_trop_weight = sum([listmonom[j0].degrees()[w]*tropical_weight[w] for w in range(nvar)])
                if term_trop_weight < max_trop_weight or (term_trop_weight == max_trop_weight and tropical_order_ring(listmonom[j0])> tropical_order_ring(max_monomial)):
                    max_monomial = listmonom[j0]
                    max_index = j0
                    max_trop_weight = term_trop_weight
            return [max_monomial,s-1-i]

    def indexes(self):
        """
        INPUT:

        OUTPUT:

        The indexes of the lines of this matrix.

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [1,0,0,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]), None, 2,[0,1])
            sage: M.indexes()
            [0, 1]
        """
        return self._indexes

    def index(self,i):
        """
        INPUT:

        - ``i`` -- an integer

        OUTPUT:

        The index of the i-th line of this matrix.

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]), None, 2,[2,1])
            sage: M.index(0)
            2
        """
        return self._indexes[i]

    def monomials(self):
        """
        INPUT:

        OUTPUT:

        The monomials attached to this matrix

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]), None, 2,[2,1])
            sage: M.monomials()
            [x^2, x*y, y^2]
        """
        return self._monomials

    def polynomial(self, i):
        """
        INPUT:
        - ``i`` -- an integer

        OUTPUT:

        The polynomial corresponding to line of index i of the Macaulay matrix

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]), None, 2,[2,1])
            sage: M.polynomial(1)
            2*x*y + 3*y^2
        """
        row = self._entries.row(i)
        monomials = self._monomials
        return sum([ row[j]*monomials[j] for j in range(self._indexes[i], self._ncols) ])

    def copy(self):
        """
        Return a copy of self

        TESTS::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]), None, 2,[2,1])
            sage: M.polynomial(1) == M.copy().polynomial(1)
            True
        """
        if self._labels is not None:
            labels = Matrix(self._base, self.nrows(), self._s ,self._labels.list())
        else:
            labels = None
        indexes = self._indexes
        if indexes is not None:
            indexes = [copy(x) for x in self._indexes]
        l = Matrix(self._base.base_ring(), self.nrows(), self.ncols() ,self._entries.list())
        monom = [ x for x in self._monomials]
        signs = None
        if self._signatures is not None:
            signs = [y for y in self._signatures]
        return self.__class__(self._base, monom, l, labels, signs, self._s, indexes)

    def addf(self, f, flab=None, fsign=None):
        """
        Add the polynomial f with label flab and signature fsign at the bottom of self

        Note::

            We assume that f has the correct degree and labels and signatures are correct

        INPUT:

        - ``f`` -- a polynomial
        - ``flab`` -- its label
        - ``fsign`` -- its signature

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]), signatures = [[x,0],[y,1]], s=2,indexes= [2,1])
            sage: M.addf(x*y-y^2, [y,x], [x,1])
            sage: M.nrows()
            3
            sage: M.index(2)
            1
            sage: M.polynomial(2)
            x*y - y^2
        """
        list = f.monomials()
        monom = self._monomials
        l = len(monom)
        flist = [ self._zero_coeff ] * l
        ind = l+1
        for mon in list:
            if mon in monom:
                i = monom.index(mon)
                flist[i] = f.monomial_coefficient(mon)
                ind = min(i,ind)
        fmat = Matrix(self._base.base_ring(), 1, l, flist)
        self._entries = self._entries.stack(fmat)
        self._indexes.append(ind)
        if self._labels is not None:
            self._labels = self._labels.stack(Matrix(self._base, 1, self._s, flab))
        if self._signatures is not None:
            self._signatures.append(fsign)

    def add_multiple_of_row(self, i, rpiv, kpiv, eliminate = False):
        """
        INPUT:

        - ``i`` -- a natural number
        - ``rpiv`` -- a natural number
        - ``kpiv`` -- an element of self.base_ring()
        - ``eliminate`` -- a boolean

        OUTPUT:

        Performs self.row(i)=self.row(i)+kpiv*self.row(rpiv),
        with the index of the row i being computed, and with
        the same performance on the matrix of labels.
        If eliminate is True, we are in the case of an inexact ring, where a real zero has to be put by hand.

        Note::

            When signatures are involved, we assume that this operation is consistent.

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]), None, 2,[2,1])
            sage: M.add_multiple_of_row(0,1,3)
            sage: M.polynomial(0)
            6*x*y + 13*y^2
            sage: M.index(0)
            1
        """
        ncols = self._ncols
        index = self._indexes[rpiv]
        self._entries.add_multiple_of_row(i,rpiv,kpiv)
        if eliminate:
            self._entries[i,index]=self._zero_coeff
        if (self._indexes[i] >= index) :
            j = index
            while (j < ncols and self._entries[i,j].is_zero()):
                j += 1
            self._indexes[i] = j
        if self._labels is not None:
            self._labels.add_multiple_of_row(i,rpiv,kpiv)

    def insert_multiple(self, Mac, i, ind, mon):
        """
        INPUT:

        - ``Mac`` -- a matrix with label
        - ``i`` -- a natural number
        - ``ind`` -- a natural number
        - ``mon`` -- a monomial

        Note::

            This operation is usually not compatible with signatures.

        OUTPUT:

        insert mon times the ind-th row of Mac to the i-th row of self

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]),None,2,[2,1])
            sage: Mac = MacaulayMatrix(T,[x,y],Matrix(QQ, 1, 2, [0,1]), Matrix(T, 1, 2, [1,0]),None,2,[1])
            sage: M.insert_multiple(Mac,1,0,x)
            sage: M.polynomial(1)
            x*y
            sage: M.index(1)
            1
            sage: M.label(1)
            (x, 0)
        """
        f=mon*Mac.polynomial(ind)
        list = f.monomials()
        monom=self._monomials
        l = len(monom)
        flist = [ self._zero_coeff ] * l
        index = l+1
        for mono in list:
            if mono in monom:
                j = monom.index(mono)
                flist[j] = f.monomial_coefficient(mono)
                index = min(j,index)
        self._entries.set_row(i,flist)

        self._indexes[i]=index
        if self._labels is not None:
            self._labels.set_row(i, [mon * lab for lab in Mac._labels.row(ind)])

    def column_transposition(self,j1,j2):
        """
        performs the transposition of the columns j1 and j2 of the Macaulay matrix self.

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]),None,2,[2,1])
            sage: M.column_transposition(1,2)
            sage: M.monomials()
            [x^2, y^2, x*y]
            sage: M.indexes()
            [1, 1]
            sage: M.polynomial(0)
            4*y^2
        """
        if j1 != j2:
            nrows = self.nrows()
            ncols = self._ncols
            temp = self._monomials[j1]
            self._monomials[j1] = self._monomials[j2]
            self._monomials[j2] = temp
            self._entries.swap_columns(j1,j2)
            for i in range(nrows):
                if self._indexes[i] != ncols:
                    if self._indexes[i]>= min(j1,j2) and self._indexes[i] <=  max(j1,j2):
                        j = min(self._indexes[i],j1,j2)
                        while (j < ncols and self._entries[i,j].is_zero()):
                            j += 1
                        self._indexes[i] = j 





    def row_echelon(self, pivot_choice = False, Macdd = None, l = None, compute_pivots = False, tropical = False, tropical_weight = None, tropical_order_ring = None, precMF5loc = None ):
        """
        Returns the row-echelon form of self

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,2,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]),None,2,[1,1])
            sage: M2 = M.row_echelon(pivot_choice = True)[0]
            sage: M2.polynomial(1)
            -y^2
            sage: M2.polynomial(0)
            2*x*y + 4*y^2
        """
        if Macdd is not None:
            weak = True
        else:
            weak = False
        if precMF5loc is not None:
            precMF5 = 0
        else:
            precMF5 = None
        P = self._base.base_ring()
        pbm_entry_non_zero = False
        to_fill = False
        M = self.copy()
        monom = self.monomials()
        entries = M._entries
        m = self.nrows()
        n = self._ncols
        nmin = min(m,n)
        rows = range(m)
        nvar = self._base.ngens()
        u = 0
        if compute_pivots: rows_pivot = [ ]
        j = 0
        if l is None:
            l = n
        while (j<l) and len(rows)>0:
            if not(tropical):
                #looking for a first real non-zero term in the i column
                maxabs = -1
                for k0 in rows:
                    entry = entries[k0,j]
                    if entry <> 0:
                        if not(pivot_choice):
                            piv = entry
                            rpiv = k0
                            maxabs = P(entry).abs()
                            break
                        else:
                            curabs = P(entry).abs()
                            if curabs > maxabs:
                                piv = entry
                                rpiv = k0
                                maxabs = curabs
                #performing the elimination
                if maxabs == -1:
                    #print("Watch out, a problem appeared with false zeroes")
                    if weak:
                        to_fill = True
                        break
                else:
                    if precMF5 is not None:
                        precMF5 = precMF5 + P(piv).valuation()
                    rows.remove(rpiv)
                    if compute_pivots: 
                        rows_pivot.append(rpiv)
                    if not(pivot_choice):
                        for k1 in range(rpiv+1,m):
                            if entries[k1,j] <> 0:
                                M.add_multiple_of_row(k1, rpiv, -entries[k1,j]/piv)
                    else:
                        for k1 in rows:
                            #here might happen the loss of precision...
                            eliminate = not(self._base.is_exact())
                            if eliminate:
                                M.add_multiple_of_row(k1, rpiv, -entries[k1,j]/piv, eliminate)
                            else:
                                if entries[k1,j] <> 0:
                                    M.add_multiple_of_row(k1, rpiv, -entries[k1,j]/piv)
                j += 1


            else:
                if pivot_choice:
                    #looking for the greatest term of the remaining rows
                    T = tropical_order_ring.base_ring()
                    max_term = None
                    max_col = None
                    max_row = None
                    entry_non_zero = None
                    pbm_entry_non_zero = False
                    for k0 in rows:
                        for j0 in range(M._indexes[k0],n):
                                if entries[k0,j0] != M._entries[k0,j0]:
                                    print("some issue has occured when handling the coefficients")

                                if max_term is None:
                                    max_term = entries[k0,j0]
                                    max_col = j0
                                    max_row = k0
                                    max_monomial = M._monomials[j0]
                                    if not(entries[k0,j0].parent().is_exact()):
                                        pbm_entry_non_zero = not(entries[k0,j0].precision_relative()>0)
                                    max_trop_weight = T(entries[k0][j0]).valuation()+sum([M._monomials[j0].degrees()[i]*tropical_weight[i] for i in range(nvar)])
                                elif entries[k0,j0].parent().is_exact() and entries[k0,j0] != 0:
                                    term_trop_weight = T(entries[k0,j0]).valuation()+sum([M._monomials[j0].degrees()[i]*tropical_weight[i] for i in range(nvar)])
                                    if term_trop_weight < max_trop_weight or (term_trop_weight == max_trop_weight and tropical_order_ring(M._monomials[j0])> tropical_order_ring(max_monomial)):
                                        max_term = entries[k0,j0]
                                        max_col = j0
                                        max_row = k0
                                        max_trop_weight = term_trop_weight
                                        max_monomial = M._monomials[j0]
                                elif not(entries[k0,j0].parent().is_exact()):
                                    term_trop_weight = T(entries[k0,j0]).valuation()+sum([M._monomials[j0].degrees()[i]*tropical_weight[i] for i in range(nvar)])
                                    if term_trop_weight < max_trop_weight or (term_trop_weight == max_trop_weight and tropical_order_ring(M._monomials[j0])> tropical_order_ring(max_monomial)):
                                        max_term = entries[k0,j0]
                                        max_col = j0
                                        max_row = k0
                                        pbm_entry_non_zero = not(entries[k0,j0].precision_relative()>0)
                                        max_trop_weight = term_trop_weight
                                        max_monomial = M._monomials[j0]
                    if pbm_entry_non_zero:
                        print("Error, the leading term could not be determined")
                        if compute_pivots:
                            return M, rows_pivot, pbm_entry_non_zero, precMF5
                        else:
                            return M, pbm_entry_non_zero, precMF5
                    elif max_term == 0 or max_term is None:
                        print("Beware of the result, there are rows reduced to zero")
                        if compute_pivots:
                            return M, rows_pivot, pbm_entry_non_zero, precMF5
                        else:
                            return M, pbm_entry_non_zero, precMF5
                    else:
                        M.column_transposition(j, max_col)
                        rows.remove(max_row)
                        if compute_pivots:
                            rows_pivot.append(max_row)
                        for k1 in rows:
                            #here might happen the loss of precision...
                            piv = max_term
                            if M._entries[k1,j] <> 0:
                                M.add_multiple_of_row(k1, max_row, -M._entries[k1,j]/piv)
                else:
                    #looking for the pivot on the u-th row
                    T = tropical_order_ring.base_ring()
                    max_term = None
                    max_col = None
                    entry_non_zero = None
                    pbm_entry_non_zero = False
                    k0=u
                    for j0 in range(M._indexes[k0],n):
                            if entries[k0,j0] != M._entries[k0,j0]:
                                print("Some issue has occured when handling the coefficients")

                            if max_term is None:
                                max_term = entries[k0,j0]
                                max_col = j0
                                max_monomial = M._monomials[j0]
                                if not(entries[k0,j0].parent().is_exact()):
                                    pbm_entry_non_zero = not(entries[k0,j0].precision_relative()>0)
                                max_trop_weight = T(entries[k0,j0]).valuation()+sum([M._monomials[j0].degrees()[i]*tropical_weight[i] for i in range(nvar)])
                            elif entries[k0,j0].parent().is_exact() and entries[k0,j0] != 0:
                                term_trop_weight = T(entries[k0,j0]).valuation()+sum([M._monomials[j0].degrees()[i]*tropical_weight[i] for i in range(nvar)])
                                if term_trop_weight < max_trop_weight or (term_trop_weight == max_trop_weight and tropical_order_ring(M._monomials[j0])> tropical_order_ring(max_monomial)):
                                    max_term = entries[k0,j0]
                                    max_col = j0
                                    max_trop_weight = term_trop_weight
                                    max_monomial = M._monomials[j0]
                            elif not(entries[k0,j0].parent().is_exact()):
                                term_trop_weight = T(entries[k0,j0]).valuation()+sum([M._monomials[j0].degrees()[i]*tropical_weight[i] for i in range(nvar)])
                                if term_trop_weight < max_trop_weight or (term_trop_weight == max_trop_weight and tropical_order_ring(M._monomials[j0])> tropical_order_ring(max_monomial)):
                                    max_term = entries[k0,j0]
                                    max_col = j0
                                    max_trop_weight = term_trop_weight
                                    max_monomial = M._monomials[j0]                                    
                                    pbm_entry_non_zero = not(entries[k0,j0].precision_relative()>0)
                    if pbm_entry_non_zero:
                        print("Error, the leading term could not be determined")
                        if compute_pivots:
                            return M, rows_pivot, pbm_entry_non_zero, precMF5
                        else:
                            return M, pbm_entry_non_zero, precMF5
                    else:
                        if max_term == 0 or max_term is None:
                            print("Beware of the result, there is a row that has reduced to zero")
                            rows.remove(u)
                        else:
                            M.column_transposition(j, max_col)
                            rows.remove(u)
                            if compute_pivots: 
                                rows_pivot.append(u)   
                            for k1 in range(u+1,m):
                                piv = max_term
                                if M._entries[k1,j] <> 0:                            
                                    M.add_multiple_of_row(k1, u, -M._entries[k1,j]/piv)
                    u += 1
                j += 1

        
        #performing a filling with Mac_{d-1,i}
        if to_fill:
            P = self._base
            ndd =  Macdd.nrows()
            ncoldd = Macdd.ncols()
            list_monomials_dd = [ ]
            for k in range(ndd):
                if Macdd.index(k) != ncoldd:
                    for x in P.gens():
                        list_monomials_dd.append([k,x])
            
            monomdd = Macdd._monomials
            dict_monom_dd = {}
            for ll in list_monomials_dd:
                dict_monom_dd[(ll[1] * monomdd[Macdd.index(ll[0])])] = ll
            
            #performing the filling
            while j < n and len(rows)>0:
                if dict_monom_dd.has_key(monom[j]):
                    a = dict_monom_dd[monom[j]]
                    i = rows.pop()
                    M.insert_multiple(Macdd, i, a[0], a[1])
                j += 1
            if len(rows)>0:
                print("We tried a weak-filling, and it was not enough")
                weak_was_enough = False
            else:
                weak_was_enough = True
        else:
            weak_was_enough = True
        if compute_pivots:
            if weak:
                return M, rows_pivot, weak_was_enough, precMF5
            elif tropical:
                return M, rows_pivot, pbm_entry_non_zero, precMF5
            else:
                return M, rows_pivot, precMF5
        else:
            if weak:
                return M, weak_was_enough, precMF5
            elif tropical:
                return M, pbm_entry_non_zero, precMF5
            else:
                return M, precMF5

    def F5crit(self,mon):
        """
        INPUT:

        - ``mon`` -- a monomial

        OUTPUT: 

        True if a row in self has mon as leading monomial, False elsewhere.

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]),None,2,[2,1])
            sage: M.F5crit(y^2)
            True
            sage: M.F5crit(x^2)
            False
        """
        return self._monomials.index(mon) in self._indexes

    def update_mac(self, d, i, fi, di, Macdi, allmonomials, Macdd=None, use_criteria=True, tropical = False, tropical_weight = None, tropical_order_ring = None, grevlex = False):
        """
        INPUT:

        - ``d`` -- a natural number
        - ``i`` -- a natural number
        - ``fi`` -- a polynomial
        - ``di`` -- a natural number
        - ``Macdi`` -- a labeled matrix
        - ``allmonomials`` -- the list of the monomials of the ambient multi-polynomial ring of all degree up to d.
        - ``Macdd`` -- the Macaulay matrix Mac[d-1][i], in the case of fast MatrixF5 algorithm
        - ``grev`` -- a boolean, to decide whether we want to use the fastest trick for graded reverse lexicographical orderings

        OUTPUT:

        If Maccdd is None, the function adds to the bottom of self all the `x^\alpha fi`, with
            `\vert \alpha \vert + deg(fi) = d`
        and `x^\alpha` does not verify the F5 criteria (unless use_criteria is False),
        ordered by label.
        Else, the rows added are created by multiples of rows of Mac[d-1][i] that pass the F5 criteria.

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), Matrix(T, 2, 3, [x,0,x-y,x+y,y,x-2*y]),None,3,[2,1])
            sage: allmonomials = [[],[x,y],[x^2,x*y,y^2]]
            sage: Macdi = MacaulayMatrix(T,[x,y],Matrix(QQ, 2, 2, [0,1,0,1]), Matrix(T, 2, 3,[1,0,1,-1,0,2]),None,3,[2,2])
            sage: M.update_mac(2,3,x,1,Macdi,allmonomials)
            sage: M.nrows()
            4
            sage: M.polynomial(2)
            x*y
        """
        if Macdd is None:
            monom = allmonomials[d-di]
            P = self._base
            l = len(monom)
            labels = self._labels
            s = self._s
            for k in range(l):
                mon = monom[l-1-k]
                if not (use_criteria and Macdi.F5crit(mon)):
                    if labels is None:
                        self.addf(mon*fi)
                    else:
                        lab = [ self._zero_poly ] * s
                        lab[i-1] = mon
                        self.addf(mon*fi, lab)
        else:
            variables = allmonomials[1]
            nvar = len(variables)
            n = Macdd.nrows()
            if grevlex:
                for k in range(n):
                    sign = Macdd.signature(k)
                    if sign[1] == i-1 and Macdd.index(k) != Macdd.ncols():
                        indmax = firstvar(sign[0])
                        if indmax < 0:
                            indmax = nvar-1 
                        for l in range(indmax+1):
                            x = variables[indmax-l]
                            if not( Macdi.F5crit(sign[0]*x)):
                                if self._labels is not None:
                                    if self._signatures is not None:
                                        nsign = [aa for aa in Macdd.signature(k)]
                                        nsign[0] = x*nsign[0]
                                        self.addf( x*Macdd.polynomial(k), [x*lab for lab in Macdd.label(k)], nsign)
                                    else:
                                        self.addf( x*Macdd.polynomial(k), [x*lab for lab in Macdd.label(k)])
                                else: #If we are here with no label, there is necessarily usage of signatures
                                    nsign = [aa for aa in Macdd.signature(k)]
                                    nsign[0] = x*nsign[0]
                                    self.addf( x*Macdd.polynomial(k), fsign = nsign)
            else:
                listsign = [ ]
                dict_sign ={ }
                #We establish here a list of all the admissible (for the F5 criterion) signatures in degree d, along with a way to obtain them from the Macaulay matrix of degree d-1
                for k in range(n):
                    if tropical:
                        sign = Macdd.signature_tropical(k, tropical_weight, tropical_order_ring)
                    else:
                        sign = Macdd.signature(k)
                    if sign[1] == i-1 and Macdd.index(k) != Macdd.ncols():
                        for x in variables:
                            if not( sign[0]*x in listsign) and not( Macdi.F5crit(sign[0]*x)):         
                                lensign = len(listsign)
                                mon = sign[0]*x
                                dict_sign.update({mon:(x,k)})
                                if lensign == 0:
                                    listsign.append(mon)
                                else:
                                    found = False
                                    if tropical:
                                        trop_weight0 = sum([mon.degrees()[vv]*tropical_weight[vv] for vv in range(nvar)])
                                        for uu in range(lensign):
                                            trop_weight = sum([listsign[uu].degrees()[ww]*tropical_weight[ww] for ww in range(nvar)])
                                            if trop_weight < trop_weight0 or (trop_weight == trop_weight0 and tropical_order_ring(listsign[uu])> tropical_order_ring(mon)):
                                                listsign.insert(uu,mon)
                                                found = True
                                                break
                                    else:
                                        for uu in range(lensign):
                                            if mon < listsign[uu]:
                                                listsign.insert(uu,mon)
                                                found = True
                                                break
                                    if not(found):
                                        listsign.append(mon)
                for mono in listsign:
                    y,row = dict_sign[mono]
                    if self._labels is not None:
                        if self._signatures is not None:
                            nsign = [aa for aa in Macdd.signature(row)]
                            nsign[0] = y*nsign[0]
                            self.addf( y*Macdd.polynomial(row), [y*lab for lab in Macdd.label(row)], nsign)
                        else:
                            self.addf( y*Macdd.polynomial(row), [y*lab for lab in Macdd.label(row)])
                    else: #If we are here with no label, there is necessarily usage of signatures
                        nsign = [aa for aa in Macdd.signature(k)]
                        nsign[0] = y*nsign[0]
                        self.addf( y*Macdd.polynomial(row), fsign = nsign)






    def macaulay_inter_reduction(self,f, flab = None, tropical = False):
        """
        INPUT:

        - ``self`` --  is a Matrix With Label, which is supposed to be under row-echelon form
        - ``f`` -- a polynomial with the same degree as self, which is supposed to come from a row of self

        OUTPUT:

        - the inter-reduction of f by self (the leading monomial of f is not reduced)

        EXAMPLES::

            sage: T.<x,y> = QQ[]
            sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
            sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,0,4,0,2,3]), indexes = [2,1])
            sage: M.macaulay_inter_reduction(x*y+3*y^2)
            x*y
        """
        fmonom = f. monomials()
        monom = self._monomials
        nmonom = self.ncols()
        dict_monom ={}
        for i in range(self.nrows()):
            if self._indexes[i] < nmonom:
                dict_monom[monom[self._indexes[i]]]=i
        g = f
        if flab is not None:
           glab = [copy(lab) for lab in flab]
        if tropical:
            found_proposed_lm = False
            for mon in monom:
                if mon in g.monomials():
                    if dict_monom.has_key(mon):
                        if found_proposed_lm:
                            i = dict_monom[u]
                            coef = g.monomial_coefficient(mon)/self._entries[i, self.index(i)]
                            g = g - coef* self.polynomial(i)
                            if flab is not None:
                                glab =  [glab[j] + coef * self._labels[j] for j in range(s)]
                        else:
                            glm = mon
                            glc = g.monomial_coefficient(mon)
            g = (1/glc)*g
            if flab is not None:
                glab = [(1/glc)*lab for lab in glab]
                return g,glab
            return g
          
        else:
            g = copy(f)
            r = g.lt()
            g = g - r
            if flab is not None:
                s = len(flab)
                rlab = [copy(lab) for lab in flab]
            while g.degree() >= 0:
                u = g.lm()
                if dict_monom.has_key(u):
                    i = dict_monom[u]
                    coef = (g.lc()/self._entries[i, self.index(i)] )
                    g = g - coef* self.polynomial(i)
                    if flab is not None:
                        rlab =  [rlab[j] + coef * self._labels[j] for j in range(s)] 
                else:
                    r = r + g.lt()
                    g = g - g.lt()

            if flab is not None:
                rlab = [(1/r.lc())*lab for lab in rlab]
            r = (1/r.lc())*r

            if flab is not None:
                return r, rlab
            else:
                return r


def update_G(Mlab, Mactilda, G, label_of_output = None, tropical = False):
    r"""
    INPUT:

    - ``Mlab`` -- a matrix with label
    - ``Mactilda`` -- required to be the row echelon form of Mlab
    - ``G`` -- a list of polynomials
    - ``label_of_output`` -- a boolean to decide if we want the labels of the elements of G

    OUTPUT:

    Add to G the polynomials corresponding to rows of mactilda such that they
    have a leading monomial not in Mlab.

    EXAMPLES::

        sage: T.<x,y> = QQ[]
        sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
        sage: G = [ ]
        sage: M = MacaulayMatrix(T,[x^2,x*y,y^2],Matrix(QQ, 2, 3, [0,2,4,0,2,3]), Matrix(T, 2, 2, [x,0,x+y,y]),None,2,[1,1])
        sage: M2 = M.row_echelon(pivot_choice = True)[0]
        sage: update_G(M, M2, G, False)
        sage: G
        [-y^2]
    """
    monom = Mlab.monomials()
    n = Mlab.nrows()
    m = Mlab.ncols()
    for i in range(n):
        if tropical:
            if Mactilda.index(i) < m and Mlab.monomials()[Mlab.index(i)] <> Mactilda.monomials()[Mactilda.index(i)]  :
                if label_of_output:
                    G.append([Mactilda.polynomial(i),Mactilda.label(i)])
                else:
                    G.append(Mactilda.polynomial(i))
        else:
            if Mlab.index(i) <> Mactilda.index(i) and Mactilda.index(i) < m:
                if label_of_output:
                    G.append([Mactilda.polynomial(i),Mactilda.label(i)])
                else:
                    G.append(Mactilda.polynomial(i))

def toy_matrixF5(F, D, P, label = True, fast = True, weak_filling = False, matrices = False, use_criteria=True, minimal = True, reduced = False, label_of_output = False, tropical = False, tropical_weight = None, tropical_order_ring = None, prec_theo = None):
    r"""
    INPUT:

    - ``F`` -- the list of the polynomials [f_1, ..., f_s]
    - ``D`` -- a natural number
    - ``P`` -- the multipolynomial ring in which f_1,...,f_s live. 
    - ``label`` -- a boolean to decide whether the algorithm works with labels or not
    - ``fast`` -- a boolean to decide whether the algorithm builds Mac_{d,i} with Mac_{d-1,i} (fast = True in the non-tropical setting does not work with an ordering which is not graded-reverse lexicographical) or with f_i (fast = False). Using the option fast triggers the usage of signatures for the Macaulay matrices that will be built
    - ``weak_filling`` -- a boolean to decide whether we want to try to fill Mac_{d,i} with multiples of rows from Mac_{d-1,i}, when a column with no non-zero pivot is encountered
    - ``matrices`` -- a boolean to decide whether we return the matrices Mac_{d,i} and Mactilda_{d,i} or not.
    - ``minimal`` -- a boolean to decide whether the algorithm returns a minimal basis or not
    - ``reduced`` -- a boolean to decide whether the algorithm returns a reduced basis or not
    - ``tropical`` -- a boolean to decide whether we work with tropical bases or Groebner bases
    - ``label_of_output`` -- a boolean to decide whether we return the expresion of the polynomials returned by the algorithm in function of the initial polynomials

    OUTPUT:

    G, a Groebner basis of `Id(F)`, if `D` is large enough, and there is no issue with precision

    EXAMPLES::

        sage: T.<x,y,z> = QQ[]
        sage: from sage.rings.polynomial.padics.toy_MatrixF5 import*
        sage: list = [x^3 - 2*x*y*z, x^2*y - 2*y^2*z + x*z^2]
        sage: D = 5
        sage: toy_matrixF5(list, 5,T, reduced = True)
        [x^3 - 2*x*y*z, x^2*y - 2*y^2*z + x*z^2, x^2*z^2, y^2*z^3 - 1/2*x*z^4, x*y*z^3]
        sage: toy_matrixF5(list, 5,T, label = False, reduced = True)
        [x^3 - 2*x*y*z, x^2*y - 2*y^2*z + x*z^2, x^2*z^2, y^2*z^3 - 1/2*x*z^4, x*y*z^3]
        sage: S.<a,b,c> = PolynomialRing(Qp(7), order='lex')
        sage: Mac, Mactilda, G = toy_matrixF5(list,5,T, matrices = True, use_criteria = True, minimal = True, label = False, tropical = True, tropical_weight = [1,2,4], tropical_order_ring = S, fast = False)
        sage: G
        [x^3 - 2*x*y*z, x^2*y - 2*y^2*z + x*z^2, x^2*z^2, 2*y^2*z^3 - x*z^4, 2*x*y*z^3]
        sage: S.<a,b,c> = PolynomialRing(Qp(2))
        sage: list = [2*x+z, x^2+y^2-2*z^2, 4*y^2+y*z+8*z^2]
        sage: Mac, Mactilda, G = toy_matrixF5(list ,4,T, matrices = True, minimal = True, use_criteria = True, tropical = True, tropical_weight = [0,0,0], tropical_order_ring = S)
        sage: G
        [2*x + z, x^2 + y^2 - 2*z^2, -2*x*y + 60/7*y^2, 1786/49*y^3]
        sage: S.<a,b,c> = PolynomialRing(Qp(7), order='lex')
        sage: toy_matrixF5(list,D,T, matrices = False, use_criteria = True, minimal = True, reduced = True, label = False, tropical = True, tropical_weight = [7,0,-3], tropical_order_ring = S, fast = False)
        [2*x + z, 1/4*x^2 + 1/4*y^2 + x*z, -30*x^2 + x*y, x^3]
        sage: toy_matrixF5(list,D,T, matrices = False, use_criteria = True, minimal = True, reduced = True, label = False, tropical = True, tropical_weight = [0,7,-3], tropical_order_ring = S, fast = False)
        [2*x + z, x^2 - 1/7*y^2, x*y - 30/7*y^2, y^3]
        sage: R.<x,y,z,w> = QQ[]
        sage: I = ideal(x^2 + y^2 - z^2 - w^2, x-y)
        sage: J = I^2
        sage: J.basis
        [x^4 + 2*x^2*y^2 + y^4 - 2*x^2*z^2 - 2*y^2*z^2 + z^4 - 2*x^2*w^2 - 2*y^2*w^2 + 2*z^2*w^2 + w^4, x^3 - x^2*y + x*y^2 - y^3 - x*z^2 + y*z^2 - x*w^2 + y*w^2, x^3 - x^2*y + x*y^2 - y^3 - x*z^2 + y*z^2 - x*w^2 + y*w^2, x^2 - 2*x*y + y^2]
        sage: GB=toy_matrixF5([x^2 - 2*x*y + y^2, x^3 - x^2*y + x*y^2 - y^3 - x*z^2 + y*z^2 - x*w^2 + y*w^2, x^3 - x^2*y + x*y^2 - y^3 - x*z^2 + y*z^2 - x*w^2 + y*w^2, x^4 + 2*x^2*y^2 + y^4 - 2*x^2*z^2 - 2*y^2*z^2 + z^4 - 2*x^2*w^2 - 2*y^2*w^2 + 2*z^2*w^2 +w^4],5,R,reduced=True)
        sage: GB.reverse()
        sage: J.groebner_basis()==GB
        True

    ALGORITHM:

    The following algorithm is adapted from [Vac2014] and [Vac2015]_.

    It consists of the Matrix-F5, its variant for computation over complete discrete valuation fields, and its tropical version. Different options, depending on the way the matrices are built or reduced are available. It produces different behaviour regarding to time-complexity and precision. 

    REFERENCES:

    .. [Vac2014] T. Vaccon. Matrix-F5 algorithms over Finite-Precision Complete Discrete Valuation Fields, proceedings of ISSAC 2014. Extended version in the Journal of Symbolic Computation.
    .. [Vac2015] T. Vaccon. MMatrix-F5 Algorithms and Tropical Gröbner Bases Computation, proceedings of ISSAC 2015.

    """
    pivot_choice = not(fast)
    zero_poly = P(0)
    one_poly = P(1)
    weak_was_enough = True
    tropical_pbm_entry_non_zero = False
    exactness = P.is_exact()

    if prec_theo is not None:
        perte_prec_MF5 = 0
    else:
        perte_prec_MF5 = None

    allmonomials = [ ]
    G = copy(F)
    s = len(F)
    Mac = [ ]
    Mactilda = [ ]
    for d in range(D+1):
        monomials = list_monom(P,d)
        lmon = len(monomials)
        allmonomials.append(monomials)
        l = [ ]
        ltilda = [ ]
        for j in range(s+1):
            if label:
                if fast:
                    mac = MacaulayMatrix(P, monomials, Matrix(P.base_ring(), 0, lmon, []), Matrix(P.base_ring(), 0, lmon, []), [], s, [])
                else:
                    mac = MacaulayMatrix(P, monomials, Matrix(P.base_ring(), 0, lmon, []), Matrix(P.base_ring(), 0, lmon, []), s=s, indexes= [])
            else:
                if fast:
                    mac = MacaulayMatrix(P, monomials, Matrix(P.base_ring(), 0, lmon, []), signatures =[], s=s, indexes= [])
                else:
                    mac = MacaulayMatrix(P, monomials, Matrix(P.base_ring(), 0, lmon, []), indexes = [])
            l.append(mac)
            ltilda.append(mac.copy())
        Mac.append(l)
        Mactilda.append(ltilda)

    B = [ ]
    deglist = [ ]
    for i in range(s):
        if label:
            ei = [ zero_poly ] * s
            ei[i] = one_poly
            if fast:
                B.append([ Matrix(P,1,s,ei), [1,i], F[i] ])
            else:
                B.append([ Matrix(P,1,s,ei), None, F[i] ])
        else:
            if fast:
                B.append([None, [1,i], F[i] ])
            else:
                B.append( [None, None, F[i] ])
        deglist.append(sum((F[i].lm()).degrees()))

    if label_of_output:
        G=[[G[i], B[i][0]] for i in range(s)]

    #Mac is totally set up.
    #now, we begin the loop on the degree.
    #be aware that the i for Mac corresponds to i-1 for F, B, deglist.

    for d in range(deglist[0],D+1):


        for i in range(1,s+1):


            Mac[d][i] = Mactilda[d][i-1].copy()

            if d == deglist[i-1]:
                Mac[d][i].addf(B[i-1][2], B[i-1][0], B[i-1][1])
            elif d > deglist[i-1]:
                if fast:
                    if tropical:
                        Mac[d][i].update_mac(d, i, F[i-1], deglist[i-1], Mactilda[d-deglist[i-1]][i-1], allmonomials,Mactilda[d-1][i], tropical = tropical, tropical_weight = tropical_weight, tropical_order_ring = tropical_order_ring  )    
                    else:
                        Mac[d][i].update_mac(d, i, F[i-1], deglist[i-1], Mactilda[d-deglist[i-1]][i-1], allmonomials,Mactilda[d-1][i], use_criteria)
                else:
                    Mac[d][i].update_mac(d, i, F[i-1], deglist[i-1], Mactilda[d-deglist[i-1]][i-1], allmonomials, None, use_criteria)
            if weak_filling:
                Mactilda[d][i], weak_was_enough_inter, precMF5loc = Mac[d][i].row_echelon(pivot_choice,Mactilda[d-1][i], precMF5loc = prec_theo)
                weak_was_enough = weak_was_enough and weak_was_enough_inter
                if prec_theo is not None:
                    perte_prec_MF5 = max(perte_prec_MF5,precMF5loc)
            elif tropical:
                Mactilda[d][i], tropical_pbm_entry_non_zero0, precMF5loc  = Mac[d][i].row_echelon(pivot_choice, None, tropical = tropical, tropical_weight = tropical_weight, tropical_order_ring = tropical_order_ring, precMF5loc = prec_theo  )    
                tropical_pbm_entry_non_zero = tropical_pbm_entry_non_zero or tropical_pbm_entry_non_zero0
                if prec_theo is not None:
                    perte_prec_MF5 = max(perte_prec_MF5,precMF5loc)
            else:
                Mactilda[d][i], precMF5loc = Mac[d][i].row_echelon(pivot_choice, None, tropical = tropical, tropical_weight = tropical_weight, tropical_order_ring = tropical_order_ring, precMF5loc = prec_theo)
                if prec_theo is not None:
                    perte_prec_MF5 = max(perte_prec_MF5,precMF5loc)



            update_G(Mac[d][i], Mactilda[d][i], G, label_of_output, tropical)



    if minimal or reduced:
        if label_of_output:
            list0 = [[copy(i) for i in y[1]] for y in G]
            list1 = [[copy(G[i][0]),list0[i]] for i in range(len(G))]
            listpol = [list1[i][0] for i in range(len(G))] 
        else:
            listpol = [copy(i) for i in G]
        i = 0
        listlm = []
        if tropical:
            T = tropical_order_ring
            Tbase = T.base_ring()
            nvar = T.ngens()
            for u in range(len(listpol)):
                polo = listpol[u]
                monpolo=polo.monomials()
                lmpolo = monpolo[0]
                max_trop_weight = Tbase(polo.monomial_coefficient(lmpolo)).valuation()+sum([lmpolo.degrees()[ss]*tropical_weight[ss] for ss in range(nvar)])
                for v in range(len(monpolo)-1):
                    mono = monpolo[v+1]
                    term_trop_weight = Tbase(polo.monomial_coefficient(mono)).valuation()+sum([mono.degrees()[ss]*tropical_weight[ss] for ss in range(nvar)])
                    if term_trop_weight < max_trop_weight or (term_trop_weight == max_trop_weight and tropical_order_ring(mono)> tropical_order_ring(lmpolo)):
                        lmpolo = mono
                        max_trop_weight = term_trop_weight
                listlm.append(lmpolo)
        else:
            listlm = [polo.lm() for polo in listpol]    
                

        while i < len(listlm):
            j = 0
            flag = True
            #since we proceed degree by degree, it is enough to only consider j<i, as dividers are necessarily of degree stricly less
            a = listlm[i]
            while j < i:
                b = listlm[j]
                if (i != j) and P.monomial_divides(b,a):
                    del listpol[i]
                    del listlm[i]
                    if label_of_output:
                        del list1[i]
                    flag = False
                    break
                j = j + 1
            if flag:
                i = i + 1
        if label_of_output:
            listpol = list1
        G = listpol
    if reduced:
        if label_of_output:
            G = [ Mactilda[g[0].degree()][s].macaulay_inter_reduction(g[0],g[1],tropical) for g in G]
        else:
            G = [ Mactilda[g.degree()][s].macaulay_inter_reduction(g, tropical = tropical) for g in G]


    #Returning the results

    if prec_theo is not None:
        if weak_filling:
            if matrices:
                return Mac, Mactilda, G, weak_was_enough, perte_prec_MF5
            else:
                return G, weak_was_enough, perte_prec_MF5
        elif tropical:
            if matrices:
                return Mac, Mactilda, G, tropical_pbm_entry_non_zero0,perte_prec_MF5
            else:
                return G, tropical_pbm_entry_non_zero0, perte_prec_MF5
        else:
            if matrices:
                return Mac, Mactilda, G,perte_prec_MF5
            else:
                return G,perte_prec_MF5

    else:
        if weak_filling:
            if matrices:
                return Mac, Mactilda, G, weak_was_enough
            else:
                return G, weak_was_enough
        elif tropical:
            if matrices:
                if exactness:
                    return Mac, Mactilda, G
                else:
                    return Mac, Mactilda, G, tropical_pbm_entry_non_zero0
            else:
                if exactness:
                    return G
                else:
                    return G, tropical_pbm_entry_non_zero0
        else:
            if matrices:
                return Mac, Mactilda, G
            else:
                return G
