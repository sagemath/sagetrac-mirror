r"""
Complete Discrete Valuation Rings (CDVR) and Fields (CDVF)
"""
#**************************************************************************
#  Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#**************************************************************************


from sage.misc.abstract_method import abstract_method

from sage.categories.category_singleton import Category_singleton
from .discrete_valuation import DiscreteValuationRings, DiscreteValuationFields
#from sage.misc.cachefunc import cached_method
from copy import copy


class CompleteDiscreteValuationRings(Category_singleton):
    """
    The category of complete discrete valuation rings

    EXAMPLES::

        sage: Zp(7) in CompleteDiscreteValuationRings()
        True
        sage: QQ in CompleteDiscreteValuationRings()
        False
        sage: QQ[['u']] in CompleteDiscreteValuationRings()
        True
        sage: Qp(7) in CompleteDiscreteValuationRings()
        False
        sage: TestSuite(CompleteDiscreteValuationRings()).run()
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationRings().super_categories()
            [Category of discrete valuation rings]
        """
        return [DiscreteValuationRings()]

    class ParentMethods:
        def _matrix_echelonize(self, M, transformation=False, exact=True):
            """
            Row-echelonize this matrix

            INPUT:

            - ``transformation`` -- a boolean; indicates whether the
              transformation matrix is returned

            - ``exact`` -- a boolean (default: ``True``);
              if ``True``, the echelon form will be as exact as possible;
              if ``False``, the transformation matrix will be exact.

            OUTPUT:

            The position of the pivots and the transformation matrix
            if asked for.

            EXAMPLES::

                sage: A = Zp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

                sage: M.echelon_form()  # indirect doctest
                [ ...0000000001  ...0000000001]
                [             0 ...00000000010]

                sage: H, L = M.echelon_form(transformation=True)  # indirect doctest
                sage: H
                [ ...0000000001  ...0000000001]
                [             0 ...00000000010]
                sage: L
                [...000000001 ...444444444]
                [...444444444 ...000000002]
                sage: L*M == H
                True

            This method works for rectangular matrices as well::

                sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
                sage: M.echelon_form()  # indirect doctest
                [ ...0000000001  ...0000000001]
                [             0 ...00000000010]
                [             0              0]

            We illustrate the role of the keyword ``exact``::

                sage: M = matrix(A, 2, 2, [2 + O(5^5), 7 + O(5^5), 1 + O(5^5), 6 + O(5^5)])
                sage: M.echelon_form(transformation=True, exact=True)  # indirect doctest
                (
                [ ...0000000001  ...0000000001]  [...0001 ...4444]
                [             0 ...00000000010], [...4444 ...0002]
                )

                sage: M = matrix(A, 2, 2, [2 + O(5^5), 7 + O(5^5), 1 + O(5^5), 6 + O(5^5)])
                sage: M.echelon_form(transformation=True, exact=False)  # indirect doctest
                (
                [...00001 ...22231]  [...0000022223             0]
                [...00000 ...00010], [...4444444444 ...0000000002]
                )

            An error is raised if the precision on the entries is
            not enough to determine the echelon form::

                sage: M = matrix(A, 2, 2, [A(0,5), 1, 5^8, 1])
                sage: M.echelon_form()  # indirect doctest
                Traceback (most recent call last):
                ...
                PrecisionError: Not enough precision to echelonize (try exact=False)

            TESTS::

            We check that it works over various rings::

                sage: from sage.rings.padics.precision_error import PrecisionError
                sage: ring1 = ZpCA(5,15)
                sage: ring2 = Zq(5^3,names='a')
                sage: ring3 = Zp(5).extension(x^2-5, names='pi')
                sage: ring4 = PowerSeriesRing(GF(5), name='t')
                sage: for A in [ ring1, ring2, ring3, ring4 ]:
                ....:     for _ in range(10):
                ....:         M = random_matrix(A,4)
                ....:         try:
                ....:             H, L = M.echelon_form(transformation=True)
                ....:         except PrecisionError:
                ....:             continue
                ....:         if L*M != H: raise RuntimeError
            """
            if exact:
                from sage.matrix.matrix_cdv import echelonize_cdv_exact
                return echelonize_cdv_exact(M, transformation, True, False)
            else:
                from sage.matrix.matrix_cdv import echelonize_cdv_nonexact
                return echelonize_cdv_nonexact(M, transformation, True)

        def _matrix_charpoly(self, M, var):
            r"""
            Return the characteristic polynomial of `M`.

            EXAMPLES::

                sage: R.<t> = PowerSeriesRing(GF(5))
                sage: M = matrix(4, 4, [ (t^(i+j)).add_bigoh(10)
                ....:                    for i in range(4) for j in range(4) ])
                sage: M
                [  1 + O(t^10)   t + O(t^10) t^2 + O(t^10) t^3 + O(t^10)]
                [  t + O(t^10) t^2 + O(t^10) t^3 + O(t^10) t^4 + O(t^10)]
                [t^2 + O(t^10) t^3 + O(t^10) t^4 + O(t^10) t^5 + O(t^10)]
                [t^3 + O(t^10) t^4 + O(t^10) t^5 + O(t^10) t^6 + O(t^10)]
                sage: M.charpoly()   # indirect doctest
                x^4 + (4 + 4*t^2 + 4*t^4 + 4*t^6 + O(t^10))*x^3

            Note that this function uses a Hessenberg-like algorithm
            that performs divisions. Hence, truncations may show up
            even if the input matrix is exact::

                sage: M = matrix(3, 3, [ 1, t, t^2, 1+t, t^2, t^3, t^2, t^3, t^4 ])
                sage: M
                [    1     t   t^2]
                [1 + t   t^2   t^3]
                [  t^2   t^3   t^4]
                sage: M.charpoly()
                x^3 + (4 + 4*t^2 + 4*t^4 + O(t^25))*x^2 + (4*t + O(t^24))*x

            Another example over the p-adics::

                sage: R = Zp(5, print_mode="digits", prec=5)
                sage: M = matrix(R, 3, 3, range(9))
                sage: M
                [        0  ...00001  ...00002]
                [ ...00003  ...00004 ...000010]
                [ ...00011  ...00012  ...00013]
                sage: M.charpoly()
                ...00001*x^3 + ...44423*x^2 + ...44412*x + ...00000
            """
            return M._charpoly_hessenberg(var)

        def _matrix_smith_form(self, M, transformation, integral, exact=True):
            r"""
            Return the Smith normal form of `M`.

            This method gets called by
            :meth:`sage.matrix.matrix2.Matrix.smith_form` to compute the Smith
            normal form over local rings and fields.

            The entries of the Smith normal form are normalized such that non-zero
            entries of the diagonal are powers of the distinguished uniformizer.

            INPUT:

            - ``M`` -- a matrix over this ring

            - ``transformation`` -- a boolean; whether the transformation matrices
              are returned

            - ``integral`` -- a subring of the base ring or ``True`` or ``None``;
              the entries of the transformation matrices are in this ring.
              If ``True``, the entries are in the ring of integers of the base
              ring; if ``None``, they are in the base ring.

            - ``exact`` -- a boolean (default: ``True``);
              if ``True``, the diagonal smith form will be exact, or raise a
              ``PrecisionError`` if this is not possible;
              if ``False``, the diagonal entries will be inexact, but the
              transformation matrices will be exact.

            EXAMPLES::

                sage: A = Zp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

                sage: S, L, R = M.smith_form()  # indirect doctest
                sage: S
                [ ...0000000001              0]
                [             0 ...00000000010]
                sage: L
                [...222222223 ...000000000]
                [...444444444 ...000000002]
                sage: R
                [...0000000001 ...2222222214]
                [            0 ...0000000001]

            If not needed, it is possible to avoid the computation of
            the transformations matrices `L` and `R`::

                sage: M.smith_form(transformation=False)  # indirect doctest
                [ ...0000000001              0]
                [             0 ...00000000010]

            This method works for rectangular matrices as well::

                sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
                sage: S, L, R = M.smith_form()  # indirect doctest
                sage: S
                [ ...0000000001              0]
                [             0 ...00000000010]
                [             0              0]
                sage: L
                [ ...222222223  ...000000000             0]
                [ ...444444444  ...000000002             0]
                [ ...444444443  ...000000001 ...0000000001]
                sage: R
                [...0000000001 ...2222222214]
                [            0 ...0000000001]

            If some of the elementary divisors have valuation larger than the
            minimum precision of any entry in the matrix, then they are
            reported as an inexact zero::

                sage: A = ZpCA(5, prec=10)
                sage: M = matrix(A, 2, 2, [5, 5, 5, 5])
                sage: M.smith_form(transformation=False, exact=False)  # indirect doctest
                [5 + O(5^10)     O(5^10)]
                [    O(5^10)     O(5^10)]

            However, an error is raised if the precision on the entries is
            not enough to determine which column to use as a pivot at some point::

                sage: M = matrix(A, 2, 2, [A(0,5), A(5^6,10), A(0,8), A(5^7,10)]); M
                [       O(5^5) 5^6 + O(5^10)]
                [       O(5^8) 5^7 + O(5^10)]
                sage: M.smith_form(transformation=False, exact=False)  # indirect doctest
                Traceback (most recent call last):
                ...
                PrecisionError: not enough precision to compute Smith normal form
            """
            from sage.matrix.matrix_cdv import smith_cdv
            integral = (integral is True or integral is None or integral is self)
            return smith_cdv(M, transformation, integral, exact)

        def _test_matrix_smith(self, **options):
            r"""
            Test that :meth:`_matrix_smith_form` works correctly.

            EXAMPLES::

                sage: ZpCA(5, 15)._test_matrix_smith()
            """
            tester = self._tester(**options)
            tester.assertEqual(self.residue_field().characteristic(), self.residue_characteristic())

            from itertools import chain
            from sage.all import MatrixSpace
            from sage.rings.padics.precision_error import PrecisionError
            matrices = chain(*[MatrixSpace(self, n, m).some_elements() for n in (1,3,7) for m in (1,4,7)])
            for M in tester.some_elements(matrices):
                bases = [self]
                if self is not self.integer_ring():
                    bases.append(self.integer_ring())
                for base in bases:
                    try:
                       S,U,V = M.smith_form(integral=base)
                    except PrecisionError:
                        continue

                    if self.is_exact() or self._prec_type() not in ['fixed-mod','floating-point']:
                        tester.assertEqual(U*M*V, S)

                    tester.assertEqual(U.nrows(), U.ncols())
                    tester.assertEqual(U.base_ring(), base)

                    tester.assertEqual(V.nrows(), V.ncols())
                    tester.assertEqual(V.base_ring(), base)

                    for d in S.diagonal():
                        if not d.is_zero():
                            tester.assertTrue(d.unit_part().is_one())

                    for (d,dd) in zip(S.diagonal(), S.diagonal()[1:]):
                        tester.assertTrue(d.divides(dd))

        def _matrix_hermite_form(self, M, include_zero_rows, transformation, integral=True, exact=True):
            """
            Return the Hermite normal form of this matrix.

            INPUT:

            - ``transformation`` -- a boolean; indicates whether the
              transformation matrix is returned

            - ``include_zero_rows`` -- a boolean

            - ``integral`` -- a subring of the base ring or ``True`` or ``None``;
              the entries of the transformation matrices are in this ring.
              If ``True``, the entries are in the ring of integers of the base
              ring; if ``None``, they are in the base ring.

            - ``exact`` -- a boolean (default: ``True``);
              if ``True``, the echelon form will be as exact as possible;
              if ``False``, the transformation matrix will be exact.

            EXAMPLES::

                sage: A = Zp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

                sage: M.hermite_form()  # indirect doctest
                [ ...0000000001  ...0000000001]
                [             0 ...00000000010]

                sage: H, L = M.hermite_form(transformation=True)  # indirect doctest
                sage: H
                [ ...0000000001  ...0000000001]
                [             0 ...00000000010]
                sage: L
                [...000000001 ...444444444]
                [...444444444 ...000000002]
                sage: L*M == H
                True

            With ``integral=False``, we get an echelon form over the fraction field::

                sage: M.hermite_form(integral=False)  # indirect doctest
                [...0000000001            0]
                [            0 ...000000001]

            We illustrate the role of the keyword ``exact``::

                sage: M = matrix(A, 2, 2, [2 + O(5^5), 7 + O(5^5), 1 + O(5^5), 6 + O(5^5)])
                sage: M.hermite_form(transformation=True, exact=True)  # indirect doctest
                (
                [ ...0000000001  ...0000000001]  [...0001 ...4444]
                [             0 ...00000000010], [...4444 ...0002]
                )
                sage: M.hermite_form(transformation=True, exact=False)  # indirect doctest
                (
                [...00001 ...22231]  [...0000022223             0]
                [...00000 ...00010], [...4444444444 ...0000000002]
                )

            .. SEEALSO::

                meth:`_matrix_echelonize`

            TESTS::

                sage: M = matrix(A, 3, 3, range(9))
                sage: M.hermite_form(exact=False)   # indirect doctest
                [ ...0000000001  ...3131313133 ...31313131320]
                [             0  ...0000000001  ...0000000002]
                [ ...0000000000  ...0000000000  ...0000000000]
                sage: M.hermite_form(include_zero_rows=False, exact=False)   # indirect doctest
                [ ...0000000001  ...3131313133 ...31313131320]
                [             0  ...0000000001  ...0000000002]

            """
            integral = (integral is True or integral is None or integral is self)
            if integral:
                M = copy(M)
            else:
                M = M.change_ring(self.fraction_field())
            if exact:
                from sage.matrix.matrix_cdv import echelonize_cdv_exact
                pivots, L = echelonize_cdv_exact(M, transformation, integral, True)
            else:
                from sage.matrix.matrix_cdv import echelonize_cdv_nonexact
                pivots, L = echelonize_cdv_nonexact(M, transformation, integral)
            rk = len(pivots)
            if not include_zero_rows and rk < M.nrows():
                M = M.submatrix(0, 0, rk)
                if transformation:
                    L = L.submatrix(0, 0, rk)
            if transformation:
                return M, L
            else:
                return M


    class ElementMethods:
        @abstract_method
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: R = Zp(7)
                sage: x = R(7); x
                7 + O(7^21)
                sage: x.valuation()
                1
            """

        def denominator(self):
            """
            Return the denominator of this element normalized
            as a power of the uniformizer

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(1/21)
                sage: x.denominator()
                7 + O(7^21)

                sage: x = K(7)
                sage: x.denominator()
                1 + O(7^20)

            Note that the denominator lives in the ring of integers::

                sage: x.denominator().parent()
                7-adic Ring with capped relative precision 20

            When the denominator is indistinguishable from 0 and the
            precision on the input is `O(p^n)`, the return value is `1`
            if `n` is nonnegative and `p^(-n)` otherwise::

                sage: x = K(0,5); x
                O(7^5)
                sage: x.denominator()
                1 + O(7^20)

                sage: x = K(0,-5); x
                O(7^-5)
                sage: x.denominator()
                7^5 + O(7^25)
            """
            return self.parent()(1)

        def numerator(self):
            """
            Return the numerator of this element, normalized in such a
            way that `x = x.numerator() / x.denominator()` always holds
            true.

            EXAMPLES::

                sage: K = Qp(7, 5)
                sage: x = K(1/21)
                sage: x.numerator()
                5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + O(7^5)

                sage: x == x.numerator() / x.denominator()
                True

            Note that the numerator lives in the ring of integers::

                sage: x.numerator().parent()
                7-adic Ring with capped relative precision 5

            TESTS::

                sage: x = K(0,-5); x
                O(7^-5)
                sage: x.numerator()
                O(7^0)
                sage: x.denominator()
                7^5 + O(7^10)
            """
            return self

        @abstract_method
        def lift_to_precision(self, absprec=None):
            """
            Return another element of the same parent with absolute precision
            at least ``absprec``, congruent to this element modulo the
            precision of this element.

            INPUT:

            - ``absprec`` -- an integer or ``None`` (default: ``None``), the
              absolute precision of the result. If ``None``, lifts to the maximum
              precision allowed.

            .. NOTE::

                If setting ``absprec`` that high would violate the precision cap,
                raises a precision error.  Note that the new digits will not
                necessarily be zero.

            EXAMPLES::

                sage: R = ZpCA(17)
                sage: R(-1,2).lift_to_precision(10)
                16 + 16*17 + O(17^10)
                sage: R(1,15).lift_to_precision(10)
                1 + O(17^15)
                sage: R(1,15).lift_to_precision(30)
                Traceback (most recent call last):
                ...
                PrecisionError: precision higher than allowed by the precision cap
                sage: R(-1,2).lift_to_precision().precision_absolute() == R.precision_cap()
                True

                sage: R = Zp(5); c = R(17,3); c.lift_to_precision(8)
                2 + 3*5 + O(5^8)
                sage: c.lift_to_precision().precision_relative() == R.precision_cap()
                True

            """


class CompleteDiscreteValuationFields(Category_singleton):
    """
    The category of complete discrete valuation fields

    EXAMPLES::

        sage: Zp(7) in CompleteDiscreteValuationFields()
        False
        sage: QQ in CompleteDiscreteValuationFields()
        False
        sage: LaurentSeriesRing(QQ,'u') in CompleteDiscreteValuationFields()
        True
        sage: Qp(7) in CompleteDiscreteValuationFields()
        True
        sage: TestSuite(CompleteDiscreteValuationFields()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationFields().super_categories()
            [Category of discrete valuation fields]
        """
        return [DiscreteValuationFields()]

    class ParentMethods:
        @abstract_method
        def integer_ring(self):
            """
            Return the integer ring of this CDVF

            EXAMPLES::

                sage: K = Qp(5)
                sage: K.integer_ring()
                5-adic Ring with capped relative precision 20

                sage: L.<t> = LaurentSeriesRing(GF(5))
                sage: L.integer_ring()
                Power Series Ring in t over Finite Field of size 5
            """

        def _matrix_hessenbergize(self, H):
            r"""
            Replace `H` with an Hessenberg form of it.

            EXAMPLES::

                sage: R.<t> = PowerSeriesRing(GF(5))
                sage: K = R.fraction_field()
                sage: H = matrix(K, 4, 4, [ (t^(i+j)).add_bigoh(10)
                ....:                       for i in range(4) for j in range(4) ])
                sage: H
                [  1 + O(t^10)   t + O(t^10) t^2 + O(t^10) t^3 + O(t^10)]
                [  t + O(t^10) t^2 + O(t^10) t^3 + O(t^10) t^4 + O(t^10)]
                [t^2 + O(t^10) t^3 + O(t^10) t^4 + O(t^10) t^5 + O(t^10)]
                [t^3 + O(t^10) t^4 + O(t^10) t^5 + O(t^10) t^6 + O(t^10)]
                sage: H.hessenbergize()
                sage: H
                [              1 + O(t^10)   t + t^3 + t^5 + O(t^10)             t^2 + O(t^10)             t^3 + O(t^10)]
                [              t + O(t^10) t^2 + t^4 + t^6 + O(t^10)             t^3 + O(t^10)             t^4 + O(t^10)]
                [                  O(t^10)                   O(t^10)                   O(t^10)                   O(t^10)]
                [                  O(t^10)                   O(t^10)                   O(t^10)                   O(t^10)]

            Another example over the p-adics::

                sage: K = Qp(5, print_mode="digits", prec=5)
                sage: H = matrix(K, 3, 3, range(9))
                sage: H
                [        0  ...00001  ...00002]
                [ ...00003  ...00004 ...000010]
                [ ...00011  ...00012  ...00013]
                sage: H.hessenbergize()
                sage: H
                [        0  ...00010  ...00002]
                [ ...00003  ...00024 ...000010]
                [ ...00000  ...44440  ...44443]
            """
            from sage.matrix.matrix_cdv import hessenbergize_cdvf
            hessenbergize_cdvf(H)

        def _matrix_echelonize(self, M, transformation=False, exact=True):
            """
            Row-echelonize this matrix

            INPUT:

            - ``transformation`` -- a boolean; indicates whether the
              transformation matrix is returned

            - ``exact`` -- a boolean (default: ``True``);
              if ``True``, the echelon form will be as exact as possible;
              if ``False``, the transformation matrix will be exact.

            OUTPUT:

            The position of the pivots and the transformation matrix
            if asked for.

            EXAMPLES::

                sage: A = Qp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 3, [2, 7, 1, 6, 3, 8])
                sage: M.echelon_form()  # indirect doctest
                [...0000000001             0 ...0441230443]
                [            0 ...0000000001 ...4330114330]

            We illustrate the role of the keyword ``exact``::

                sage: M = M.apply_map(lambda x: x.add_bigoh(3))
                sage: M.echelon_form(transformation=True, exact=True)  # indirect doctest
                (
                [...0000000001             0        ...443]  [...202 ...322]
                [            0 ...0000000001        ...330], [...041 ...433]
                )

                sage: M = M.apply_map(lambda x: x.add_bigoh(3))
                sage: M.echelon_form(transformation=True, exact=False)  # indirect doctest
                (
                [...001      0 ...443]  [...1223201202 ...4444222322]
                [     0 ...001 ...330], [...4442033041 ...0000000433]
                )

            An error is raised if the precision on the entries is
            not enough to determine the echelon form::

                sage: M = matrix(A, 2, 2, [A(0,5), 1, 5^8, 1])
                sage: M.echelon_form()  # indirect doctest
                Traceback (most recent call last):
                ...
                PrecisionError: Not enough precision to echelonize (try exact=False)
            """
            if exact:
                from sage.matrix.matrix_cdv import echelonize_cdv_exact
                return echelonize_cdv_exact(M, transformation, False, False)
            else:
                from sage.matrix.matrix_cdv import echelonize_cdv_nonexact
                return echelonize_cdv_nonexact(M, transformation, False)

        def _matrix_smith_form(self, M, transformation, integral, exact=True):
            r"""
            Return the Smith normal form of `M`.

            This method gets called by
            :meth:`sage.matrix.matrix2.Matrix.smith_form` to compute the Smith
            normal form over local rings and fields.

            The entries of the Smith normal form are normalized such that non-zero
            entries of the diagonal are powers of the distinguished uniformizer.

            INPUT:

            - ``M`` -- a matrix over this ring

            - ``transformation`` -- a boolean; whether the transformation matrices
              are returned

            - ``integral`` -- a subring of the base ring or ``True`` or ``None``;
              the entries of the transformation matrices are in this ring.
              If ``True``, the entries are in the ring of integers of the base
              ring; if ``None``, they are in the base ring.

            - ``exact`` -- a boolean (default: ``True``);
              if ``True``, the diagonal smith form will be exact, or raise a
              ``PrecisionError`` if this is not possible;
              if ``False``, the diagonal entries will be inexact, but the
              transformation matrices will be exact.

            EXAMPLES::

                sage: A = Qp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

                sage: S, L, R = M.smith_form()  # indirect doctest
                sage: S
                [...0000000001             0]
                [            0 ...0000000001]
                sage: L
                [ ...222222223  ...000000000]
                [...44444444.4 ...00000000.2]
                sage: R
                [...0000000001 ...2222222214]
                [            0 ...0000000001]

            The above is the Smith form of ``M`` over `\QQ_p`. If we desire
            a Smith form over `\ZZ_p`, we need to pass in explicitely
            ``integral=True``::

                sage: S, L, R = M.smith_form(integral=True)  # indirect doctest
                sage: S
                [...0000000001             0]
                [            0 ...0000000010]
                sage: L
                [...222222223 ...000000000]
                [...444444444 ...000000002]
                sage: R
                [...0000000001 ...2222222214]
                [            0 ...0000000001]

            If not needed, it is possible to avoid the computation of
            the transformations matrices `L` and `R`::

                sage: M.smith_form(transformation=False)  # indirect doctest
                [...0000000001             0]
                [            0 ...0000000001]

            This method works for rectangular matrices as well::

                sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
                sage: S, L, R = M.smith_form(integral=True)  # indirect doctest
                sage: S
                [ ...0000000001              0]
                [             0 ...00000000010]
                [             0              0]
                sage: L
                [ ...222222223  ...000000000             0]
                [ ...444444444  ...000000002             0]
                [ ...444444443  ...000000001 ...0000000001]
                sage: R
                [...0000000001 ...2222222214]
                [            0 ...0000000001]

            If some of the elementary divisors have valuation larger than the
            minimum precision of any entry in the matrix, then they are
            reported as an inexact zero::

                sage: A = ZpCA(5, prec=10)
                sage: M = matrix(A, 2, 2, [5, 5, 5, 5])
                sage: M.smith_form(transformation=False, exact=False)  # indirect doctest
                [5 + O(5^10)     O(5^10)]
                [    O(5^10)     O(5^10)]

            However, an error is raised if the precision on the entries is
            not enough to determine which column to use as a pivot at some point::

                sage: M = matrix(A, 2, 2, [A(0,5), A(5^6,10), A(0,8), A(5^7,10)]); M
                [       O(5^5) 5^6 + O(5^10)]
                [       O(5^8) 5^7 + O(5^10)]
                sage: M.smith_form(transformation=False, exact=False)  # indirect doctest
                Traceback (most recent call last):
                ...
                PrecisionError: not enough precision to compute Smith normal form
            """
            from sage.matrix.matrix_cdv import smith_cdv
            integral = (integral is True or integral is self.integer_ring())
            return smith_cdv(M, transformation, integral, exact)

        def _matrix_hermite_form(self, M, include_zero_rows, transformation, integral=None, exact=True):
            """
            Return the Hermite normal form of this matrix.

            INPUT:

            - ``transformation`` -- a boolean; indicates whether the
              transformation matrix is returned

            - ``include_zero_rows`` -- a boolean

            - ``integral`` -- a subring of the base ring or ``True`` or ``None``;
              the entries of the transformation matrices are in this ring.
              If ``True``, the entries are in the ring of integers of the base
              ring; if ``None``, they are in the base ring.

            - ``exact`` -- a boolean (default: ``True``);
              if ``True``, the echelon form will be as exact as possible;
              if ``False``, the transformation matrix will be exact.

            EXAMPLES::

                sage: A = Qp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 3, [1, 2, 3, 6, 7, 8])

                sage: M.hermite_form()  # indirect doctest
                [...0000000001             0  ...444444444]
                [            0 ...0000000001  ...000000002]

                The user should be very careful that the above is the
                Hermite form over `\QQ_p`; in order to get a Hermite normal
                form over `\ZZ_p`, we must pass in the argument ``integral=True``::

                sage: M.hermite_form(integral=True)  # indirect doctest
                [ ...0000000001  ...0000000002  ...0000000003]
                [             0 ...00000000010  ...0000000020]

            We can compute in addition the transformation matrix::

                sage: H, L = M.hermite_form(integral=True, transformation=True)  # indirect doctest
                sage: H
                [ ...0000000001  ...0000000002  ...0000000003]
                [             0 ...00000000010  ...0000000020]
                sage: L
                [...000000001 ...000000000]
                [...000000011 ...444444444]
                sage: L*M == H
                True

            .. SEEALSO::

                meth:`_matrix_echelonize`

            TESTS::

                sage: M = matrix(A, 3, 3, range(9))
                sage: M.hermite_form(exact=False)   # indirect doctest
                [...0000000001             0 ...4444444444]
                [            0 ...0000000001 ...0000000002]
                [            0             0 ...0000000000]
                sage: M.hermite_form(include_zero_rows=False, exact=False)   # indirect doctest
                [...0000000001             0 ...4444444444]
                [            0 ...0000000001 ...0000000002]

            """
            integral = (integral is True or integral is self.integer_ring())
            M = copy(M)
            if exact:
                from sage.matrix.matrix_cdv import echelonize_cdv_exact
                pivots, L = echelonize_cdv_exact(M, transformation, integral, True)
            else:
                from sage.matrix.matrix_cdv import echelonize_cdv_nonexact
                pivots, L = echelonize_cdv_nonexact(M, transformation, integral)
            rk = len(pivots)
            if not include_zero_rows and rk < M.nrows():
                M = M.submatrix(0, 0, rk)
                if transformation:
                    L = L.submatrix(0, 0, rk)
            if transformation:
                return M, L
            else:
                return M

    class ElementMethods:
        @abstract_method
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(7); x
                7 + O(7^21)
                sage: x.valuation()
                1
            """

        def denominator(self):
            """
            Return the denominator of this element normalized
            as a power of the uniformizer

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(1/21)
                sage: x.denominator()
                7 + O(7^21)

                sage: x = K(7)
                sage: x.denominator()
                1 + O(7^20)

            Note that the denominator lives in the ring of integers::

                sage: x.denominator().parent()
                7-adic Ring with capped relative precision 20

            When the denominator is indistinguishable from 0 and the
            precision on the input is `O(p^n)`, the return value is `1`
            if `n` is nonnegative and `p^(-n)` otherwise::

                sage: x = K(0,5); x
                O(7^5)
                sage: x.denominator()
                1 + O(7^20)

                sage: x = K(0,-5); x
                O(7^-5)
                sage: x.denominator()
                7^5 + O(7^25)
            """
            val = self.valuation()
            R = self.parent().integer_ring()
            if val >= 0:
                return R(1)
            else:
                return R(1) << (-val)

        def numerator(self):
            """
            Return the numerator of this element, normalized in such a
            way that `x = x.numerator() / x.denominator()` always holds
            true.

            EXAMPLES::

                sage: K = Qp(7, 5)
                sage: x = K(1/21)
                sage: x.numerator()
                5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + O(7^5)

                sage: x == x.numerator() / x.denominator()
                True

            Note that the numerator lives in the ring of integers::

                sage: x.numerator().parent()
                7-adic Ring with capped relative precision 5

            TESTS::

                sage: x = K(0,-5); x
                O(7^-5)
                sage: x.numerator()
                O(7^0)
                sage: x.denominator()
                7^5 + O(7^10)
            """
            R = self.parent().integer_ring()
            return R(self * self.denominator())
