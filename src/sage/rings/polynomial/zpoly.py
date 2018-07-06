r"""
Polynomial rings with variables indexed by integers

This class implements endomorphisms that are
induced from functions from the integers to themselves.
"""
from functools import partial
from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial_dense
from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing_dense
from sage.rings.integer import Integer
from sage.misc.cachefunc import cached_method
from sage.rings.fraction_field import FractionField_generic
from sage.rings.fraction_field_element import FractionFieldElement
from sage.categories.quotient_fields import QuotientFields

class PolynomialRingIntegerGenerators(InfinitePolynomialRing_dense):
    r"""
    Returns the polynomial algebra with generators indexed by the integers.

    INPUT:

    - ``base_ring`` -- (default: ZZ) a commutative base ring
    - ``prefix`` -- (default: None) a pair of strings for variables indexed by positive and
    nonpositive integers respectively
    - ``order`` -- (default: None) ordering for polynomial ring

    EXAMPLES::

        sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
        sage: R = PolynomialRingIntegerGenerators(QQ,('a','b'))
        sage: R.an_element()
        a_1
        sage: [R.igen(i) for i in range(-2,3)]
        [b_2, b_1, b_0, a_1, a_2]
    """

    def __init__(self, base_ring, prefix=None, order=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: PolynomialRingIntegerGenerators(QQ,('a','b'))
            Polynomial Ring over Rational Field with variables indexed by integers

        Issues: TestSuite fails since pure functions cannot be pickled.
        """
        if not prefix:
            prefix = ['a','b']
        if not order:
            order = 'degrevlex'
        super(PolynomialRingIntegerGenerators, self).__init__(base_ring, prefix, order)
        #InfinitePolynomialRing_sparse.__init__(self, base_ring, prefix, order)
        self._positive_prefix = prefix[0]
        self._nonpositive_prefix = prefix[1]

    def __repr__(self):
        r"""
        A string representing ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: PolynomialRingIntegerGenerators(QQ,['x','y']) # indirect doctest
            Polynomial Ring over Rational Field with variables indexed by integers
        """
        return "Polynomial Ring over {} with variables indexed by integers".format(self.base_ring())

    def igen(self, i):
        r"""
        The generator of ``self`` indexed by the integer `i`.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['x','y'])
            sage: [R.igen(i) for i in range(-2,3)]
            [y_2, y_1, y_0, x_1, x_2]
        """
        if i > 0:
            return self.gen(0)[i]
        return self.gen(1)[-i]

    @cached_method
    def variable_index(self, v):
        r"""
        Given the variable `v` of ``self``, return its integer index.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['x','y'])
            sage: R.variable_index(R.igen(-3))
            -3
            sage: R.variable_index(R.igen(0))
            0
            sage: R.variable_index(R.igen(2))
            2
        """
        s = str(v)
        i = s.find('_')
        if i == -1:
            raise ValueError("Invalid variable: no underscore")
        name = s[:i]
        value = Integer(s[(i+1):])
        if name == self._positive_prefix:
            return value
        if name == self._nonpositive_prefix:
            return -value
        raise ValueError("Invalid prefix")

    def apply_map_induced_by_func(self, func, f):
        r"""
        Apply to `f`, the endomorphism of ``self`` induced by the function `func` from integers to integers.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['a','b'])
            sage: from sage.combinat.affine_permutation import AffinePermutationTypeA
            sage: A=AffinePermutationGroup(['A',3,1]) 
            sage: p=A([-1, 0, 5, 6])
            sage: R.apply_map_induced_by_func(p, R.igen(2)**2 + R.igen(-1)*R.igen(-3))
            a_1*b_5 + b_0^2
        """
        f = self(f)
        par = f.polynomial().parent()
        if par == self.base_ring():
            return f
        return par.hom([self.igen(func(self.variable_index(v))) for v in par.gens()])(f)

    def map_induced_by_func(self, func):
        r"""
        The endomorphism of ``self`` induced by the function `func` from integers to integers.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['a','b'])
            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
            sage: p = PermutationOfIntegers(-2, 4, [0, 1, -2, 4, -1, 3, 2])
            sage: f = R.map_induced_by_func(p.func())
            sage: f(R.igen(2)**2 + R.igen(-1)*R.igen(-3))
            a_1*b_3 + b_1^2
            sage: f(R.igen(-2))
            b_0
            sage: f(R.igen(1)**2)
            a_4^2
        """
        return partial(self.apply_map_induced_by_func, func)

    @cached_method
    def shift_func(self, r):
        r"""
        The function that adds `r` to an integer.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['x','y'])
            sage: sh = R.shift_func(3)
            sage: sh(5)
            8
            sage: sh(-5)
            -2
        """
        return lambda x: x+r

    @cached_method
    def shift_auto(self, r):
        r"""
        The automorphism of ``self`` which shifts the indices of variables by `r`.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['a','b'])
            sage: sh = R.shift_auto(3)
            sage: sh(R.igen(-6))
            b_3
            sage: sh(R.igen(2)**2 + R.igen(-1)*R.igen(-3))
            a_5^2 + a_2*b_0
        """
        return self.map_induced_by_func(self.shift_func(r))

    def shift_auto_on_element(self, r, elt):
        r"""
        Compute the automorphism of shifting `elt` by `r` times.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: P = PolynomialRingIntegerGenerators(QQ, ('a','b'))
            sage: P.shift_auto_on_element(-3, P.igen(2)**2*P.igen(-3)+P.igen(0))
            b_6*b_1^2 + b_3
        """
        return self.shift_auto(r)(elt)

    @cached_method
    def swap_func(self, i):
        r"""
        The function that swaps i and i+1.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['x','y'])
            sage: sw = R.swap_func(-1)
            sage: [(i,sw(i)) for i in range(-2,2)]
            [(-2, -2), (-1, 0), (0, -1), (1, 1)]
            sage: sw = R.swap_func(0)
            sage: [(i,sw(i)) for i in range(-2,2)]
            [(-2, -2), (-1, -1), (0, 1), (1, 0)]
        """
        return lambda x: i+1 if x==i else i if x==i+1 else x

    @cached_method
    def swap_auto(self, i):
        r"""
        The automorphism of ``self`` that swaps the i-th and i+1-th variables.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['a','b'])
            sage: sw = R.swap_auto(-1)
            sage: sw(R.igen(-1)**2 + R.igen(0)*R.igen(3))
            a_3*b_1 + b_0^2
            sage: sw = R.swap_auto(3)
            sage: sw(R.igen(-1)**2 + R.igen(0)*R.igen(3)+R.igen(4))
            b_1^2 + a_4*b_0 + a_3
        """
        return self.map_induced_by_func(self.swap_func(i))

    def swap_auto_on_element(self, i, elt):
        r"""
        Compute the automorphism of swapping i and i+1.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: P = PolynomialRingIntegerGenerators(QQ, ('a','b'))
            sage: P.swap_auto_on_element(-3, P.igen(2)**2*P.igen(-3)+P.igen(0))
            a_2^2*b_2 + b_0
        """
        return self.swap_auto(i)(elt)

    def ddiff_on_element(self, i, elt):
        r"""
        Compute the divided difference on the element.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: P = PolynomialRingIntegerGenerators(QQ, ('a','b'))
            sage: P.ddiff_on_element(1, P.igen(1)**2)
            a_2 + a_1
            sage: P.ddiff_on_element(2, P.igen(1)**3*P.igen(2)**2*P.igen(3))
            a_3*a_2*a_1^3
            sage: P.ddiff_on_element(1,_)
            a_3*a_2^2*a_1 + a_3*a_2*a_1^2
        """
        num = elt - self.swap_auto(i)(elt)
        den = self.igen(i) - self.igen(i+1)
        par = (num*den).polynomial().parent()
        nump = par(num.polynomial())
        denp = par(den.polynomial())
        return self(par(nump/denp))

    @cached_method
    def dynkin_reversal_func(self):
        r"""
        The function that sends i to 1-i.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['x','y'])
            sage: df = R.dynkin_reversal_func()
            sage: [(i,df(i)) for i in range(-2,2)]
            [(-2, 3), (-1, 2), (0, 1), (1, 0)]
        """
        return lambda x: 1-x

    @cached_method
    def dynkin_reversal_auto(self):
        r"""
        The automorphism of ``self`` induced by the dynkin_reversal_func.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: R = PolynomialRingIntegerGenerators(QQ,['a','b'])
            sage: df = R.dynkin_reversal_auto()
            sage: df(R.igen(-1)**2 + R.igen(0)*R.igen(3))
            a_2^2 + a_1*b_2
        """
        return self.map_induced_by_func(self.dynkin_reversal_func())

    def dynkin_reversal_auto_on_element(self, elt):
        r"""
        This does the above automorphism but also sends each variable to its negative.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: P = PolynomialRingIntegerGenerators(QQ, ('a','b'))
            sage: P.dynkin_reversal_auto_on_element(P.igen(2)**2*P.igen(-3)+P.igen(0))
            -a_4*b_1^2 - a_1
            sage: P.dynkin_reversal_auto_on_element(P.igen(2)**2+P.igen(0))
            b_1^2 - a_1
        """
        return self.negate_variables(self.dynkin_reversal_auto()(elt))

    def negate_variables(self, f):
        r"""
        Negate all variables in f.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
            sage: P = PolynomialRingIntegerGenerators(QQ, ('a','b'))
            sage: P.negate_variables(P.igen(2)**2*P.igen(-3)+P.igen(0)+P.igen(-2)**2)
            -a_2^2*b_3 + b_2^2 - b_0
        """
        fp = f.polynomial()
        if fp.parent() == self.base_ring():
            return f
        vars = fp.variables()
        return self(fp.subs(dict([[v,-v] for v in vars])))

    def forget_to_loop_rotation(self, d, f, positive=True):
        r"""
        Apply to `f` the forgetful map to the polynomial ring with one generator.

        INPUT:

        - `f` -- a polynomial in ``self``
        - `d` -- a parameter
        - ``positive`` -- optional (default: True) If True, send positive variables to `-d` and others to 0.
          If False, send positive variables to 0 and others to `d`.

        sage: from sage.rings.polynomial.zpoly import PolynomialRingIntegerGenerators
        sage: R = PolynomialRingIntegerGenerators(QQ,['a','b'])
        sage: S = PolynomialRing(QQ, ['d'])
        sage: f = (R.igen(2)-R.igen(-1))**2; f
        a_2^2 - 2*a_2*b_1 + b_1^2
        sage: R.forget_to_loop_rotation(S.gen(0), f)
        d^2
        sage: R.forget_to_loop_rotation(S.gen(0), f, positive=False)
        d^2
        """
        vars = f.variables()
        f = f.polynomial()
        if len(vars) == 0:
            return f
        zed = d.parent().zero()
        if positive:
            pos_image = -d
            nonpos_image = zed
        else:
            pos_image = zed
            nonpos_image = d
        return f(tuple([nonpos_image if self.variable_index(self(v)) <= 0 else pos_image for v in f.parent().gens()]))

class FunctionFieldIntegerGeneratorsElement(FractionFieldElement):
    def act_by_func(self, func):
        r"""
        Act on ``self`` by the automorphism induced by the function ``func`` on integers.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: Q = FunctionFieldIntegerGenerators(QQ, ('a','b'))
            sage: c = Q.igen(2)/(Q.igen(-2)**2+1); c
            a_2/(b_2^2 + 1)
            sage: sh = Q.polynomial_ring().shift_func(3)
            sage: c.act_by_func(sh)
            a_5/(a_1^2 + 1)
        """
        return self.parent().apply_map_induced_by_func(func, self)

    def fix(self):
        denom = self.denominator()
        if len(denom.variables()) == 0:
            return (1/denom)*self.numerator()
        return self

    def forget_to_loop_rotation(self, d, positive=True):
        r"""
        Apply to `f` the forgetful map to the polynomial ring with one generator.

        INPUT:

        - `d` -- a parameter
        - ``positive`` -- optional (default: True) If True, send positive variables to `-d` and others to 0.
          If False, send positive variables to 0 and others to `d`.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: Q = FunctionFieldIntegerGenerators(QQ,('a','b'))
            sage: S = PolynomialRing(QQ, ['d'])
            sage: ((Q.igen(4)**2-Q.igen(-3))/(2+Q.igen(2))).forget_to_loop_rotation(S.gen(0))
            d^2/(-d + 2)
            sage: ((Q.igen(4)-Q.igen(-3))/(2+Q.igen(2))).forget_to_loop_rotation(S.gen(0), positive=False)
            -1/2*d
        """
        P = self.parent().polynomial_ring()
        return P.forget_to_loop_rotation(d, P(self.numerator()), positive=positive)/P.forget_to_loop_rotation(d, P(self.denominator()), positive=positive)

class FunctionFieldIntegerGenerators(FractionField_generic):
    r"""
    Fraction field of Polynomial ring with variables indexed by integers.

    This class supports endomorphisms induced by maps from the integers to themselves.
    """
    Element = FunctionFieldIntegerGeneratorsElement

    def __init__(self, base_ring, prefix):
        r"""
        Initialize ``self``.
        """
        R = PolynomialRingIntegerGenerators(base_ring, prefix)
        FractionField_generic.__init__(self, R, element_class=FunctionFieldIntegerGeneratorsElement,category=QuotientFields())
        self._R = R
        
    def polynomial_ring(self):
        r"""
        The polynomial ring of which ``self`` is the fraction field.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: FunctionFieldIntegerGenerators(QQ, ('a','b')).polynomial_ring()
            Polynomial Ring over Rational Field with variables indexed by integers
        """
        return self._R

    def __repr__(self):
        r"""
        A string representing ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: FunctionFieldIntegerGenerators(QQ, ('a','b')) # indirect doctest
            Function Field over Rational Field with generators indexed by the integers
        """
        return "Function Field over {} with generators indexed by the integers".format(self.polynomial_ring().base_ring())

    @cached_method
    def an_element(self):
        r"""
        An element of ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: FunctionFieldIntegerGenerators(QQ, ('a','b')).an_element()
            1/a_1
        """
        return self.one()/self(self.polynomial_ring().an_element())

    @cached_method
    def igen(self, i):
        r"""
        Return the `i`-th generator of ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: Q = FunctionFieldIntegerGenerators(QQ, ('a','b'))
            sage: [Q.igen(-3), Q.igen(0), Q.igen(4)]
            [b_3, b_0, a_4]
        """
        return self(self.polynomial_ring().igen(i))

    @cached_method
    def variable_index(self, v):
        r"""
        Given a variable `v`, return its integer index.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: Q = FunctionFieldIntegerGenerators(QQ, ('a','b'))
            sage: Q.variable_index(Q.igen(4))
            4
            sage: Q.variable_index(Q.igen(-3))
            -3
        """
        return self.polynomial_ring().variable_index(v)

    def apply_map_induced_by_func(self, func, f):
        r"""
        Act by the endomorphism of ``self`` induced by the function ``func`` from integers to integers.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: Q = FunctionFieldIntegerGenerators(QQ, ('a','b'))
            sage: c = Q.igen(2)/(Q.igen(-2)**2+1)
            sage: Q.apply_map_induced_by_func(Q.polynomial_ring().shift_func(3), c)
            a_5/(a_1^2 + 1)
        """
        R = self.polynomial_ring()
        Rmap = R.map_induced_by_func(func)
        return self(Rmap(f.numerator()))/self(Rmap(f.denominator()))

    def map_induced_by_func(self, func):
        r"""
        The endomorphism of ``self`` induced by the function ``func`` from integers to integers.
        """
        return partial(self.apply_map_induced_by_func, func)

    def shift_auto_on_element(self, r, elt):
        r"""
        Compute the automorphism of shifting `elt` by `r` times.

        EXAMPLES::

            sage: from sage.rings.polynomial.zpoly import FunctionFieldIntegerGenerators
            sage: Q = FunctionFieldIntegerGenerators(QQ, ('a','b'))
            sage: Q.shift_auto_on_element(-3, Q.igen(2)/Q.igen(-5))
            b_1/b_8
        """
        rshift = lambda i: i+r
        return self.apply_map_induced_by_func(rshift, elt)
