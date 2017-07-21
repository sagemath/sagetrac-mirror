"""
Helper methods for matrices over CDVR/CDVF
"""

from sage.rings.infinity import Infinity
from sage.rings.padics.precision_error import PrecisionError

from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields


def smith_normal_form(M, transformation):
    """
    Helper method for the computation of the Smith normal form

    See also :meth:`_matrix_smith_form`
    """
    n = M.nrows()
    m = M.ncols()
    S = M.parent()(M.list())
    smith = M.parent()(0)
    R = M.base_ring()
    if R.tracks_precision():
        precM = min([ x.precision_absolute() for x in M.list() ])

    if transformation:
        from sage.matrix.special import identity_matrix
        left = identity_matrix(R,n)
        right = identity_matrix(R,m)
    else:
        left = right = None

    val = -Infinity
    for piv in range(min(n,m)):
        curval = Infinity
        pivi = pivj = piv
        for i in range(piv,n):
            for j in range(piv,m):
                v = S[i,j].valuation()
                if v < curval:
                    pivi = i; pivj = j
                    curval = v
                    if v == val: break
            else:
                continue
            break
        val = curval

        if R.tracks_precision() and precM is not Infinity and val >= precM:
            raise PrecisionError("Not enough precision to compute Smith normal form")

        if val is Infinity:
            break

        S.swap_rows(pivi,piv)
        S.swap_columns(pivj,piv)
        if transformation:
            left.swap_rows(pivi,piv)
            right.swap_columns(pivj,piv)

        smith[piv,piv] = R(1) << val
        inv = ~(S[piv,piv] >> val)
        for i in range(piv+1,n):
            scalar = -inv * (S[i,piv] >> val)
            if R.tracks_precision():
                scalar = scalar.lift_to_maximal_precision()
            S.add_multiple_of_row(i,piv,scalar,piv+1)
            if transformation:
                left.add_multiple_of_row(i,piv,scalar)
        if transformation:
            left.rescale_row(piv,inv)
            for j in range(piv+1,m):
                scalar = -inv * (S[piv,j] >> val)
                if R.tracks_precision():
                    scalar = scalar.lift_to_maximal_precision()
                right.add_multiple_of_column(j,piv,scalar)

    if transformation:
        if R.tracks_precision() and precM is not Infinity:
            left = left.apply_map(lambda x: x.add_bigoh(precM-val))
        return smith, left, right
    else:
        return smith


def echelonize(M, transformation, secure):
    """
    Helper method to echelonize matrices over CDVR/CDVF

    See also :meth:`_matrix_echelonize`

    TODO:

    Analyse better precision.
    """
    print "Echelonize"
    print M
    n = M.nrows()
    m = M.ncols()
    R = M.base_ring()
    if R in CompleteDiscreteValuationFields():
        Rint = R.integer_ring()
    else:
        Rint = R

    if transformation:
        from sage.matrix.special import identity_matrix
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
            if v < val:
                pivi = ii
                val = v

        if val is Infinity:
            j += 1
            continue

        if R.tracks_precision():
            not_enough_precision = False
            for ii in range(n):
                if M[ii,j].precision_absolute() <= val:
                    not_enough_precision = True
            if not_enough_precision:
                # In this situation, we do not know for sure what
                # is the pivot. 
                if secure:
                    raise PrecisionError("Not enough precision to echelonize")
                # When secure is False, if all entries on the column
                # are inexact zeroes, we go ahead
                for ii in range(i,n):
                    if M[ii,j] != 0:
                        raise PrecisionError("Not enough precision to echelonize")
                j += 1
                continue

        pivots.append(j)

        M.swap_rows(pivi,i)
        if transformation:
            left.swap_rows(pivi,i)

        inv = ~(M[i,j] >> val); inv = R(inv)
        M.rescale_row(i,inv,j+1)
        M[i,j] = R(1) << val
        if transformation:
            left.rescale_row(i, inv)

        for ii in range(i+1,n):
            scalar = -(M[ii,j] >> val)
            M.add_multiple_of_row(ii,i,scalar,j+1)
            M[ii,j] = R(0)
            if transformation:
                left.add_multiple_of_row(ii,i,scalar)

        for ii in range(i):
            v = min(0, M[ii,j].valuation())
            quo = Rint(M[ii,j] >> v) << (v-val)
            M.add_multiple_of_row(ii,i,-quo,j)
            if R.tracks_precision():
                M[ii,j] = M[ii,j].lift_to_maximal_precision()
            if transformation:
                left.add_multiple_of_row(ii,i,-quo)

        i += 1; j += 1

    return pivots, left
