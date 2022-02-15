"""
Matrices over complete discrete valuation rings/fields

AUTHOR:

- Xavier Caruso
"""

from sage.matrix.matrix_generic_dense import Matrix_generic_dense

from sage.rings.infinity import Infinity
from sage.rings.padics.precision_error import PrecisionError

from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields



def echelonize(M, transformation, secure):
    """
    Helper method to echelonize matrices over CDVR/CDVF

    See also :meth:`_matrix_echelonize`

    TODO:

    Analyse better precision.
    """
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






class Matrix_cdvf_dense(Matrix_generic_dense):
    def integral_smith_form(self, transformation=True):
        """
        Return the integral Smith normal form of this matrix.

        INPUT:

        - ``transformation`` -- a boolean (default: True)
          Indicates whether the transformation matrices are returned

        NOTE:

        The integral Smith decomposition of a matrix `M`
        defined over a complete discrete valuation field
        is a writing of the form `L*M*R = S` where:

        - `L` and `R` are invertible matrices over the ring of
          integers

        - the only non-vanishing entries of `S` are located on
          the diagonal (through `S` might be not a square matrix)

        - if `d_i` denotes the `(i,i)` entry of `S`, then `d_i`
          divides `d_{i+1}` in the ring of integers for all `i`.

        The `d_i`'s are uniquely determined provided that they are
        normalized so that they are all either `0` or a power of the 
        distinguished uniformizer of the base ring.

        EXAMPLES::

            sage: A = Qp(5, prec=10, print_mode="digits")
            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

            sage: S, L, R = M.integral_smith_form()
            sage: S
            [ ...1     0]
            [    0 ...10]
            sage: L
            [...222222223          ...]
            [...444444444         ...2]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        If not needed, it is possible to do not compute the
        transformations matrices L and R as follows::

            sage: M.integral_smith_form(transformation=False)
            [ ...1     0]
            [    0 ...10]

        This method works for rectangular matrices as well::

            sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
            sage: S, L, R = M.integral_smith_form()
            sage: S
            [ ...1     0]
            [    0 ...10]
            [    0     0]
            sage: L
            [...222222223          ...          ...]
            [...444444444         ...2          ...]
            [...444444443         ...1         ...1]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        An error is raised if the precision on the entries is
        not enough to determine the Smith normal form::

            sage: M = matrix(A, 2, 2, [1, 1, 1, 1])
            sage: M.integral_smith_form()
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision to compute Smith normal form

        TESTS::

        We check that Smith decomposition works over various rings::

            sage: from sage.rings.padics.precision_error import PrecisionError
            sage: ring1 = Qp(7,10)
            sage: ring2 = Qq(7^2,names='a')
            sage: ring3 = Qp(7).extension(x^3-7, names='pi')
            sage: ring4 = LaurentSeriesRing(GF(7), name='t')
            sage: for A in [ ring1, ring2, ring4 ]:  # ring3 causes troubles (see ticket #23464)
            ....:     for _ in range(10):
            ....:         M = random_matrix(A,4)
            ....:         try:
            ....:             S, L, R = M.integral_smith_form()
            ....:         except PrecisionError:
            ....:             continue
            ....:         if L*M*R != S: raise RuntimeError
        """
        return smith_normal_form(self, transformation)


    def integral_echelonize(self, transformation=False, secure=False):
        """
        Row-echelonize this matrix using a transformation matrix
        which is invertible over the ring of integers

        INPUT:

        - ``transformation`` -- a boolean (default: False)
          Indicates whether the transformation matrix is returned

        - ``secure`` -- a boolen (default: False)
          When ``secure`` is ``True``, we raise an error whenever
          a pivot cannot be determined for sure.
          Otherwise, if a column contains only non-exact zeroes,
          we go ahead.

        OUTPUT:

        The transformation matrix if asked for.

        EXAMPLES::

            sage: A = Qp(5, prec=10, print_mode="digits")
            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

            sage: M.integral_echelonize()
            sage: M
            [ ...1  ...1]
            [    0 ...10]

            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])
            sage: Mc = M.__copy__()
            sage: L = M.integral_echelonize(transformation=True)
            sage: M
            [ ...1  ...1]
            [    0 ...10]
            sage: L
            [        ...1 ...444444444]
            [...444444444         ...2]
            sage: L*Mc == M
            True

        This method works for rectangular matrices as well::

            sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
            sage: M.integral_echelonize()
            sage: M
            [ ...1  ...1]
            [    0 ...10]
            [    0     0]

        Depending on the value of ``secure``, an error is
        raised if the precision on the entries is not enough 
        to determine the echelon form::

            sage: M = matrix(A, 2, 2, [1, 1, 1, 1])
            sage: M.integral_echelonize()  # by default, secure=False
            sage: M
            [...1 ...1]
            [   0  ...]

            sage: M = matrix(A, 2, 2, [1, 1, 1, 1])
            sage: M.integral_echelonize(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision to echelonize

        TESTS::

        We check that it works over various rings::

            sage: from sage.rings.padics.precision_error import PrecisionError
            sage: ring1 = Qp(7,10)
            sage: ring2 = Qq(7^2,names='a')
            sage: ring3 = Qp(7).extension(x^3-7, names='pi')
            sage: ring4 = LaurentSeriesRing(GF(7), name='t')
            sage: for A in [ ring1, ring2, ring4 ]:  # ring3 causes troubles (see ticket #23464)
            ....:     for _ in range(10):
            ....:         M = random_matrix(A,4)
            ....:         Mc = M.__copy__()
            ....:         try:
            ....:             L = M.integral_echelonize(transformation=True)
            ....:         except PrecisionError:
            ....:             continue
            ....:         if L*Mc != M: raise RuntimeError
        """
        _, left = echelonize(self, transformation, secure=secure)
        if transformation:
            return left

    def integral_echelon_form(self, transformation=False, secure=False):
        """
        Return the row-echelon form of this matrix using a transformation
        matrix which is invertible over the ring of integers

        INPUT:

        - ``transformation`` -- a boolean (default: False)
          Indicates whether the transformation matrix is returned

        - ``secure`` -- a boolen (default: False)
          When ``secure`` is ``True``, we raise an error whenever
          a pivot cannot be determined for sure.
          Otherwise, if a column contains only non-exact zeroes,
          we go ahead.

        OUTPUT:

        The row-echelon from and the transformation matrix if asked for.

        EXAMPLES::

            sage: A = Qp(5, prec=10, print_mode="digits")
            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

            sage: M.integral_echelon_form()
            [ ...1  ...1]
            [    0 ...10]

            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])
            sage: H,L = M.integral_echelon_form(transformation=True)
            sage: H
            [ ...1  ...1]
            [    0 ...10]
            sage: L
            [        ...1 ...444444444]
            [...444444444         ...2]
            sage: L*M == H
            True

        This method works for rectangular matrices as well::

            sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
            sage: M.integral_echelon_form()
            [ ...1  ...1]
            [    0 ...10]
            [    0     0]

        Depending on the value of ``secure``, an error is
        raised if the precision on the entries is not enough 
        to determine the echelon normal form::

            sage: M = matrix(A, 2, 2, [1, 1, 1, 1])
            sage: M.integral_echelon_form()  # by default, secure=False
            [...1 ...1]
            [   0  ...]

            sage: M = matrix(A, 2, 2, [1, 1, 1, 1])
            sage: M.integral_echelon_form(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision to echelonize

        TESTS::

        We check that it works over various rings::

            sage: from sage.rings.padics.precision_error import PrecisionError
            sage: ring1 = Qp(7,10)
            sage: ring2 = Qq(7^2,names='a')
            sage: ring3 = Qp(7).extension(x^3-7, names='pi')
            sage: ring4 = LaurentSeriesRing(GF(7), name='t')
            sage: for A in [ ring1, ring2, ring4 ]:  # ring3 causes troubles (see ticket #23464)
            ....:     for _ in range(10):
            ....:         M = random_matrix(A,4)
            ....:         try:
            ....:             H,L = M.integral_echelon_form(transformation=True)
            ....:         except PrecisionError:
            ....:             continue
            ....:         if L*M != H: raise RuntimeError
        """
        E = self.__copy__()
        _, left = echelonize(E, transformation, secure=secure)
        if transformation:
            return E, left
        else:
            return E
