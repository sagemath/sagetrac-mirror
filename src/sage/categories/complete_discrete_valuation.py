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
        def _matrix_echelonize(self, M, transformation=False, exact=False):
            """
            Row-echelonize this matrix

            INPUT:

            - ``transformation`` -- a boolean (default: True)
              Indicates whether the transformation matrix is returned

            OUTPUT:

            The position of the pivots and the transformation matrix
            if asked for.

            EXAMPLES::

                sage: A = Zp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

                sage: M.echelon_form()  # indirect doctest
                [ ...1  ...1]
                [    0 ...10]

                sage: H,L = M.echelon_form(transformation=True)  # indirect doctest
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
                sage: M.echelon_form()  # indirect doctest
                [ ...1  ...1]
                [    0 ...10]
                [    0     0]

            An error is raised if the precision on the entries is
            not enough to determine the echelon form::

                sage: M = matrix(A, 2, 2, [A(0,5), 1, 5^8, 1])
                sage: M.echelon_form()  # indirect doctest
                Traceback (most recent call last):
                ...
                PrecisionError: Not enough precision to echelonize

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
            from sage.matrix.matrix_cdv import echelonize_cdv
            return echelonize_cdv(M, transformation, True, exact)

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

        def _matrix_smith_form(self, M, transformation, integral=True, exact=True):
            from sage.matrix.matrix_cdv import smith_cdv
            return smith_cdv(M, transformation, integral, exact)

        def _matrix_hermite_form(self, M, transformation, integral=True, exact=True):
            from sage.matrix.matrix_cdv import echelonize_cdv
            return echelonize_cdv(copy(M), transformation, integral, exact)

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
            from .precision_error import PrecisionError
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

        def _matrix_echelonize(self, M, transformation=False, integral=False, exact=False):
            from sage.matrix.matrix_cdv import echelonize_cdv
            return echelonize_cdv(M, transformation, integral, exact)

        def _matrix_smith_form(self, M, transformation, integral=True, exact=True):
            from sage.matrix.matrix_cdv import smith_cdv
            return smith_cdv(M, transformation, integral, exact)

        def _matrix_hermite_form(self, M, transformation, integral=True, exact=True):
            from sage.matrix.matrix_cdv import echelonize_cdv
            return echelonize_cdv(copy(M), transformation, integral, exact)

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
