r"""
Special methods for matrices over discrete valuation rings/fields.
"""

# ****************************************************************************
#       Copyright (C) 2020 Xavier Caruso <xavier.caruso@normalesup.org>
#                          Raphaël Pagès <raphael.pages@u-bordeaux.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#          https://www.gnu.org/licenses/
# ****************************************************************************


from sage.structure.element cimport RingElement
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.matrix.special import identity_matrix

from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields

from sage.rings.padics.precision_error import PrecisionError
from copy import copy


# Echelon form (Hermite form)
#############################

def echelonize_cdv_nonexact(M, transformation, integral):
    r"""
    Echelonize the matrix `M`.

    This method gets called by
    :meth:`sage.matrix.matrix2.Matrix.echelon_form` and
    :meth:`sage.matrix.matrix2.Matrix.hermite_form`.

    INPUT:

    - ``M`` -- a matrix over this ring

    - ``transformation`` -- a boolean; whether the transformation matrices
      are returned

    - ``integral`` -- a boolean; whether the echelon form should be computed
      over the ring of integers or its fraction field

    OUTPUT:

    A pair ``(pivots, transformation)`` where ``pivots`` is a list of
    pivots and ``transformation`` is the transformation matrix or ``None``
    if not asked.

    In this variant, the transformation matrix is exact (i.e. known at
    the maximal precision) while the echelon form may have inexact entries.

    TESTS::

        sage: R = Zp(3)
        sage: M = matrix(R, 3, 3, [ R(i,5) for i in range(9) ])
        sage: M.echelon_form(exact=False)  # indirect doctest
        [    3 + O(3^5) 1 + 3 + O(3^5) 2 + 3 + O(3^5)]
        [        O(3^5)     1 + O(3^5)     2 + O(3^5)]
        [        O(3^5)         O(3^5)         O(3^5)]

    """
    n = M.nrows()
    m = M.ncols()
    R = M.base_ring()
    if integral:
        if R.is_field():
            Rint = R.integer_ring()
        else:
            Rint = R

    if transformation:
        left = identity_matrix(R,n)
    else:
        left = None

    i = j = 0
    pivots = [ ]
    while i < n and j < m:
        val = Infinity
        pivi = i
        for ii in range(i,n):
            v = M[ii,j].valuation()
            if val is Infinity or v < val:
                pivi = ii
                val = v
        if val is Infinity:
            j += 1
            continue

        if any(M[ii,j].precision_absolute() <= val for ii in range(n)):
            if any(M[ii,j] != 0 for ii in range(i, n)):
                raise PrecisionError("Not enough precision to echelonize")
            else:
                j += 1
                continue

        pivots.append(j)

        M.swap_rows(pivi,i)
        if transformation:
            left.swap_rows(pivi,i)

        if integral:
            inv = (M[i,j] >> val).inverse_of_unit()
            inv = inv.lift_to_precision()
            M.rescale_row(i, inv, j)
            if transformation:
                left.rescale_row(i, inv)
            for ii in range(i+1,n):
                scalar = -(M[ii,j] >> val)
                scalar = scalar.lift_to_precision()
                M.add_multiple_of_row(ii, i, scalar, j)
                if transformation:
                    left.add_multiple_of_row(ii, i, scalar)
        else:
            inv = ~M[i,j]
            inv = inv.lift_to_precision()
            M.rescale_row(i, inv, j)
            if transformation:
                left.rescale_row(i, inv)
            for ii in range(n):
                if ii == i:
                    continue
                scalar = -M[ii,j].lift_to_precision()
                M.add_multiple_of_row(ii, i, scalar, j)
                M[ii,j] = R(0)
                if transformation:
                    left.add_multiple_of_row(ii, i, scalar)

        i += 1; j += 1

    return pivots, left


def echelonize_cdv_exact(M, transformation, integral, secure):
    r"""
    Echelonize the matrix `M`.

    This method gets called by
    :meth:`sage.matrix.matrix2.Matrix.echelon_form` and
    :meth:`sage.matrix.matrix2.Matrix.hermite_form`.

    INPUT:

    - ``M`` -- a matrix over this ring

    - ``transformation`` -- a boolean; whether the transformation matrices
      are returned

    - ``integral`` -- a boolean; whether the echelon form should be computed
      over the ring of integers or its fraction field

    - ``secure`` -- a boolean; if ``True``, raise an error if all the 
      possible pivots on a column are inexact zeroes; if ``False``, ignore
      this and continue with the next column

    OUTPUT:

    a pair ``(pivots, transformation)`` where ``pivots`` is a list of
    pivots and ``transformation`` is the transformation matrix or ``None``
    if not asked.

    In this variant, if `r` is the rank of `M`, the first `r` columns of the
    echelon form are exact, i.e. given at the maximal precision allowed by 
    the parent.

    TESTS::

        sage: R = Zp(3)
        sage: M = matrix(R, 3, 3, [ R(i,5) for i in range(9) ])
        sage: M.echelon_form()  # indirect doctest
        [                 3 + O(3^21)                            0 2*3 + 2*3^2 + 2*3^3 + O(3^4)]
        [                           0                  1 + O(3^20)                   2 + O(3^4)]
        [                           0                            0                       O(3^4)]
    """
    n = M.nrows()
    m = M.ncols()
    R = M.base_ring()
    if integral:
        if R.is_field():
            Rint = R.integer_ring()
        else:
            Rint = R

    if transformation:
        left = identity_matrix(R,n)
    else:
        left = None

    i = j = 0
    pivots = [ ]
    while i < n and j < m:
        val = Infinity
        pivi = i
        for ii in range(i,n):
            v = M[ii,j].valuation()
            if val is Infinity or v < val:
                pivi = ii
                val = v
        if val is Infinity:
            j += 1
            continue

        if any(M[ii,j].precision_absolute() <= val for ii in range(n)):
            if secure or any(M[ii,j] != 0 for ii in range(i, n)):
                raise PrecisionError("Not enough precision to echelonize (try exact=False)")
            else:
                j += 1
                continue

        pivots.append(j)

        M.swap_rows(pivi,i)
        if transformation:
            left.swap_rows(pivi,i)

        if integral:
            inv = (M[i,j] >> val).inverse_of_unit()
            M[i,j] = R(1) << val
            M.rescale_row(i, inv, j+1)
            if transformation:
                left.rescale_row(i, inv)
            for ii in range(i+1,n):
                scalar = -(M[ii,j] >> val)
                M.add_multiple_of_row(ii, i, scalar, j+1)
                M[ii,j] = R(0)
                if transformation:
                    left.add_multiple_of_row(ii, i, scalar)
            for ii in range(i):
                v = min(0, M[ii,j].valuation())
                scalar = -(Rint(M[ii,j] >> v) << (v-val))
                M.add_multiple_of_row(ii, i, scalar, j)
                M[ii,j] = M[ii,j].lift_to_precision()
                if transformation:
                    left.add_multiple_of_row(ii, i, scalar)
        else:
            inv = ~M[i,j]
            M[i,j] = R(1)
            M.rescale_row(i, inv, j+1)
            if transformation:
                left.rescale_row(i, inv)
            for ii in range(n):
                if ii == i:
                    continue
                scalar = -M[ii,j]
                M.add_multiple_of_row(ii, i, scalar, j+1)
                M[ii,j] = R(0)
                if transformation:
                    left.add_multiple_of_row(ii, i, scalar)

        i += 1; j += 1

    return pivots, left


# Smith form
############

def flatten_precision(M):
    r"""
    Rescale rows and columns of ``M`` so that the minimal
    absolute precision of each row and column is equal to `0`.

    This method is useful for increasing the numerical
    stability. It is called by :func:`smith_cdv` and
    :func:`determinant_cdv`

    Only for internal use.

    OUTPUT:

    The lists of valuations by which rows and columns,
    respectively, have been shifted.

    EXAMPLES::

        sage: from sage.matrix.matrix_cdv import flatten_precision
        sage: K = Qp(2, print_mode='digits', prec=10)
        sage: M = matrix(K, 2, 2, [K(1,5),K(2,7),K(3,3),K(5,8)])
        sage: M
        [   ...00001  ...0000010]
        [     ...011 ...00000101]
        sage: flatten_precision(M)
        ([-5, -3], [0, -2])
        sage: M
        [  ...?.00001  ...?.000001]
        [    ...?.011 ...000.00101]
    """
    parent = M.base_ring()
    n = M.nrows()
    m = M.ncols()
    shift_rows = n * [ ZZ(0) ]
    shift_cols = m * [ ZZ(0) ]
    for i in range(n):
        prec = min(M[i,j].precision_absolute() for j in range(m))
        if prec is Infinity:
            continue
        shift_rows[i] = s = -prec
        for j in range(m):
            M[i,j] <<= s
    for j in range(m):
        prec = min(M[i,j].precision_absolute() for i in range(n))
        if prec is Infinity:
            continue
        shift_cols[j] = s = -prec
        for i in range(n):
            M[i,j] <<= s
    return shift_rows, shift_cols


def smith_cdv(M, transformation, integral, exact):
    r"""
    Return the Smith normal form of `M`.

    This method gets called by
    :meth:`sage.matrix.matrix2.Matrix.smith_form` to compute the Smith
    normal form over local rings and fields.

    INPUT:

    - ``M`` -- a matrix over this ring

    - ``transformation`` -- a boolean; whether the transformation matrices
      are returned

    - ``integral`` -- a boolean; whether the Smith form should be computed
      over the ring of integers or its fraction field

    - ``exact`` -- boolean.  If ``True``, the diagonal smith form will
      be exact, or raise a ``PrecisionError`` if this is not possible.
      If ``False``, the diagonal entries will be inexact, but the
      transformation matrices will be exact.

    TESTS::

        sage: A = ZpCR(5, prec=10)
        sage: M = zero_matrix(A, 2)
        sage: M.smith_form(transformation=False)  # indirect doctest
        [0 0]
        [0 0]

        sage: M = matrix(2, 2, [ A(0,10), 0, 0, 0] )
        sage: M.smith_form(transformation=False)  # indirect doctest
        Traceback (most recent call last):
        ...
        PrecisionError: some elementary divisors indistinguishable from zero (try exact=False)
        sage: M.smith_form(transformation=False, exact=False)  # indirect doctest
        [O(5^10) O(5^10)]
        [O(5^10) O(5^10)]
    """
    R = M.base_ring()
    n = M.nrows()
    m = M.ncols()
    if m > n:
        ## It's easier below if we can always deal with precision on left.
        if transformation:
            d, u, v = smith_cdv(M.transpose(), True, integral, exact)
            return d.transpose(), v.transpose(), u.transpose()
        else:
            return smith_cdv(M.transpose(), False, integral, exact).transpose()
    smith = M.parent()(0)
    S = copy(M)
    if R.is_field():
        Z = R.integer_ring()
    else:
        Z = R

    ## the difference between ball_prec and inexact_ring is just for lattice precision.
    if hasattr(R, '_prec_type'):
        ball_prec = R._prec_type() in ['capped-rel', 'capped-abs', 'relaxed']
        inexact_ring = R._prec_type() not in ['fixed-mod', 'floating-point']
    else:
        ball_prec = inexact_ring = True

    if not integral:
        shift_rows, shift_cols = flatten_precision(S)

    precS = min(x.precision_absolute() for x in S.list())
    if transformation:
        left = identity_matrix(R,n)
        right = identity_matrix(R,m)

    if ball_prec and precS is Infinity: # capped-rel and M = 0 exactly
        return (smith, left, right) if transformation else smith

    val = -Infinity
    for piv in range(m): # m <= n
        curval = Infinity
        pivi = pivj = piv
        # allzero tracks whether every possible pivot is zero.
        # if so, we can stop.  allzero is also used in detecting some
        # precision problems: if we can't determine what pivot has
        # the smallest valuation, or if exact=True and some elementary
        # divisor is zero modulo the working precision
        allzero = True
        # allexact is tracked because there is one case where we can correctly
        # deduce the exact smith form even with some elementary divisors zero:
        # if the bottom right block consists entirely of exact zeros.
        allexact = True
        for i in range(piv,n):
            for j in range(piv,m):
                Sij = S[i,j]
                v = Sij.valuation()
                allzero = allzero and Sij.is_zero()
                if exact: # we only care in this case
                    allexact = allexact and Sij.precision_absolute() is Infinity
                if v < curval:
                    pivi = i
                    pivj = j
                    curval = v
                    if v == val:
                        break
            else:
                continue
            break
        val = curval

        if inexact_ring and not allzero and val >= precS:
            if ball_prec:
                raise PrecisionError("not enough precision to compute Smith normal form")
            precS = min([ S[i,j].precision_absolute() for i in range(piv,n) for j in range(piv,m) ])
            if val >= precS:
                raise PrecisionError("not enough precision to compute Smith normal form")

        if allzero:
            if exact:
                if allexact:
                    # We need to finish checking allexact since we broke out of the loop early
                    for i in range(i,n):
                        for j in range(piv,m):
                            allexact = allexact and S[i,j].precision_absolute() is Infinity
                            if not allexact:
                                break
                        else:
                            continue
                        break
                if not allexact:
                    raise PrecisionError("some elementary divisors indistinguishable from zero (try exact=False)")
            break

        # We swap the lowest valuation pivot into position
        S.swap_rows(pivi,piv)
        S.swap_columns(pivj,piv)
        if transformation:
            left.swap_rows(pivi,piv)
            right.swap_columns(pivj,piv)

        # ... and clear out this row and column.  Note that we
        # will deal with precision later, thus the call to lift_to_precision
        smith[piv,piv] = R(1) << val
        inv = (S[piv,piv] >> val).inverse_of_unit()
        if ball_prec:
            inv = inv.lift_to_precision()
        for i in range(piv+1,n):
            scalar = -inv * Z(S[i,piv] >> val)
            if ball_prec:
                scalar = scalar.lift_to_precision()
            S.add_multiple_of_row(i,piv,scalar,piv+1)
            if transformation:
                left.add_multiple_of_row(i,piv,scalar)
        if transformation:
            left.rescale_row(piv,inv)
            for j in range(piv+1,m):
                scalar = -inv * Z(S[piv,j] >> val)
                if ball_prec:
                    scalar = scalar.lift_to_precision()
                right.add_multiple_of_column(j,piv,scalar)
    else:
        # We use piv as an upper bound on a range below, and need to set it correctly
        # in the case that we didn't break out of the loop
        piv = m
    # We update the precision on left
    # The bigoh measures the effect of multiplying by row operations
    # on the left in order to clear out the digits in the smith form
    # with valuation at least precS
    if ball_prec and exact and transformation:
        for j in range(n):
            delta = min(left[i,j].valuation() - smith[i,i].valuation() for i in range(piv))
            if delta is not Infinity:
                for i in range(n):
                    left[i,j] = left[i,j].add_bigoh(precS + delta)
    ## Otherwise, we update the precision on smith
    if ball_prec and not exact:
        smith = smith.apply_map(lambda x: x.add_bigoh(precS))
    ## We now have to adjust the elementary divisors (and precision) in the non-integral case
    if not integral:
        for i in range(piv):
            v = smith[i,i].valuation()
            if transformation:
                for j in range(n):
                    left[i,j] >>= v
            if exact:
                smith[i,i] = R(1)
            else:
                for j in range(n):
                    smith[i,j] = smith[i,j] >> v
        if transformation:
            for i in range(n):
                for j in range(n):
                    left[i,j] <<= shift_rows[j]
            for i in range(m):
                for j in range(m):
                    right[i,j] <<= shift_cols[i]
    if transformation:
        return smith, left, right
    else:
        return smith


# Determinant
#############

def determinant_cdv(M):
    r"""
    Return the determinant of the matrix `M`.

    This method gets called by
    :meth:`sage.matrix.matrix2.Matrix.determinant`.

    INPUT:

    - ``M`` -- a matrix over this ring

    ALGORITHM:

    We flatten the absolute precision in order to increase
    the numerical stability.

    We row-echelonize the matrix by always choosing the
    pivot of smallest valuation and allowing permutations
    of columns.

    Then we compute separately the value of the determinant
    (as the product of the diagonal entries of the row-echelon
    form) and a bound on the precision on it.

    EXAMPLES::

        sage: R = Qp(5,10)
        sage: M = matrix(R, 2, 2, [1, 6, 2, 7])
        sage: M.determinant()  # indirect doctest
        4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)

        sage: (5*M).determinant()  # indirect doctest
        4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12)

    Sometimes, we gain precision on the determinant::

        sage: M = matrix(R, 3, 3,
        ....:             [R(16820,7), R(73642,7), R( 3281,7),
        ....:              R(67830,7), R(63768,7), R(76424,7),
        ....:              R(37790,7), R(38784,7), R(69287,7)])
        sage: M.determinant()  # indirect doctest
        4*5^5 + 4*5^6 + 3*5^7 + 2*5^8 + O(5^9)

    TESTS:

    We check the stability of our algorithm::

        sage: for dim in range(3,10):
        ....:     M = matrix(dim, dim, [ R(1) for _ in range(dim^2) ])
        ....:     print(M.determinant())
        O(5^20)
        O(5^30)
        O(5^40)
        O(5^50)
        O(5^60)
        O(5^70)
        O(5^80)

        sage: A = random_matrix(Qp(5),4)
        sage: B = random_matrix(Qp(5),4)
        sage: (A*B).det() == A.det()*B.det()
        True
        sage: A.change_ring(QQ).det() == A.det()
        True
        sage: matrix(Qp(37),[0]).determinant()
        0
        sage: matrix(Qp(37),[O(37)]).determinant()
        O(37)
    """
    n = M.nrows()

    # For 2x2 matrices, we use the formula
    if n == 2:
        return M[0,0]*M[1,1] - M[0,1]*M[1,0]

    R = M.base_ring()
    if hasattr(R, '_prec_type'):
        track_precision = R._prec_type() in ['capped-rel', 'capped-abs', 'relaxed']
    else:
        track_precision = True

    S = M.change_ring(R.fraction_field())
    shift_rows, shift_cols = flatten_precision(S)
    shift = sum(shift_rows) + sum(shift_cols)
    det = R(1)

    sign = 1
    valdet = 0
    val = -Infinity
    for piv in range(n):
        pivi = pivj = piv
        curval = S[pivi, pivj].valuation()
        for i in range(piv,n):
            for j in range(piv,n):
                v = S[i,j].valuation()
                if v < curval:
                    pivi = i
                    pivj = j
                    curval = v
                    if v == val:
                        break
            else:
                continue
            break
        val = curval
        if S[pivi,pivj] == 0:
            if track_precision:
                return R(0, valdet + (n-piv)*val - shift)
            else:
                return R(0)

        valdet += val
        S.swap_rows(pivi,piv)
        if pivi > piv:
            sign = -sign
        S.swap_columns(pivj,piv)
        if pivj > piv:
            sign = -sign

        det *= S[piv,piv]
        inv = ~(S[piv,piv] >> val)
        for i in range(piv+1,n):
            scalar = -inv * (S[i,piv] >> val)
            if track_precision:
                scalar = scalar.lift_to_precision()
            S.add_multiple_of_row(i,piv,scalar)

    det = sign * det
    if track_precision:
        relprec = +Infinity
        relprec_neg = 0
        for i in range(n):
            prec = Infinity
            for j in range(n):
                prec = min(prec, S[i,j].precision_absolute())
            prec -= S[i,i].valuation()
            if prec < relprec:
                relprec = prec
            if prec < 0:
                relprec_neg += prec
        if relprec_neg < 0:
            relprec = relprec_neg
        if relprec is not Infinity:
            det = det.add_bigoh(valdet + relprec)
    return det >> shift



# Hesserberg form
#################

# We assume that H is square
cpdef hessenbergize_cdvf(Matrix_generic_dense H):
    r"""
    Replace `H` with an Hessenberg form of it.

    .. NOTE::

        This function assumes that H is a matrix over
        a complete discrete valuation field.

        The pivot on each column is always chosen
        with maximal relative precision, which ensures 
        the numerical stability of the algorithm.

    TESTS::

        sage: K = Qp(5, print_mode="digits", prec=5)
        sage: H = matrix(K, 3, 3, range(9))
        sage: H
        [        0  ...00001  ...00002]
        [ ...00003  ...00004 ...000010]
        [ ...00011  ...00012  ...00013]
        sage: H.hessenbergize()  # indirect doctest
        sage: H
        [        0  ...00010  ...00002]
        [ ...00003  ...00024 ...000010]
        [ ...00000  ...44440  ...44443]

    ::

        sage: M = random_matrix(K, 6, 6)
        sage: M.charpoly()[0] == M.determinant()
        True

    We check that :trac:`31753` is resolved::

        sage: R.<t> = GF(5)[[]]
        sage: M = matrix(3, 3, [ 1, t + O(t^3), t^2, 1 + t + O(t^3), 2 + t^2, 3 + 2*t + O(t^3), t - t^2, 2*t, 1 + t ])
        sage: M.charpoly()
        x^3 + (1 + 4*t + 4*t^2 + O(t^3))*x^2 + (t + 2*t^2 + O(t^3))*x + 3 + 2*t^2 + O(t^3)
    """
    cdef Py_ssize_t n, i, j, k
    cdef RingElement entry, pivot, inv, scalar

    n = H.nrows()
    for j in range(n-1):
        k = j + 1
        maxi = H.get_unsafe(k, j).precision_relative()
        i = j + 2
        while maxi is not Infinity and i < n:
            entry = H.get_unsafe(i, j)
            if entry:
                m = entry.precision_relative()
                if m > maxi:
                    maxi = m
                    k = i
            i += 1

        if k != j + 1:
            H.swap_rows_c(j+1, k)
            H.swap_columns_c(j+1, k)
        pivot = H.get_unsafe(j+1, j)
        if pivot:
            inv = ~pivot
            for i in range(j+2, n):
                scalar = inv * H.get_unsafe(i, j)
                H.add_multiple_of_row_c(i, j+1, -scalar, j)
                H.add_multiple_of_column_c(j+1, i, scalar, 0)
