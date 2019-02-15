r"""
Colored Permutations

.. TODO::

    Much of the colored permutations (and element) class can be
    generalized to `G \wr S_n`
"""
import itertools
from six.moves import range

from sage.structure.element import parent, Element
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod

from sage.arith.all import factorial
from sage.combinat.permutation import Permutations
from sage.matrix.constructor import diagonal_matrix
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.all import ZZ


class DecoratedPermutation(Element):
    """
    A decorated permutation.
    """
    def __init__(self, parent, colors, perm):
        """
        Initialize ``self``.
        
        TESTS::
            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: TestSuite(s1*s2*t).run()
        """
        self._colors = tuple(colors)
        self._perm = perm
        self._numfp = perm.number_of_fixed_points()
        self._fp = perm.fixed_points()
        Element.__init__(self, parent=parent)

    def __hash__(self):
        r"""
        TESTS::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: hash(s1), hash(s2), hash(t)
            (2666658751600856334, 3639282354432100950, 3639281107336048003) # 64-bit
            (-1973744370, 88459862, -1467077245)                            # 32-bit
        """
        return hash(self._perm) ^ hash(self._colors)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: s1*s2*t
            [[1, 0, 0], [3, 1, 2]]
        """
        return repr([list(self._colors), self._perm])

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: latex(s1*s2*t)
            [3_{1}, 1_{0}, 2_{0}]
        """
        ret = "["
        ret += ", ".join("{}_{{{}}}".format(x, self._colors[i])
                         for i, x in enumerate(self._perm))
        return ret + "]"

#    def _mul_(self, other):
#        """
#        Multiply ``self`` and ``other``.
#
#        EXAMPLES::
#
#            sage: C = ColoredPermutations(4, 3)
#            sage: s1,s2,t = C.gens()
#            sage: s1*s2*s1 == s2*s1*s2
#            True
#        """
#        colors = tuple(self._colors[i] + other._colors[val - 1]  # -1 for indexing
#                       for i, val in enumerate(self._perm))
#        p = self._perm._left_to_right_multiply_on_right(other._perm)
#        return self.__class__(self.parent(), colors, p)

#    def inverse(self):
#        """
#        Return the inverse of ``self``.
#
#        EXAMPLES::
#
#            sage: C = ColoredPermutations(4, 3)
#            sage: s1,s2,t = C.gens()
#            sage: ~t
#            [[0, 0, 3], [1, 2, 3]]
#            sage: all(x * ~x == C.one() for x in C.gens())
#            True
#        """
#        ip = ~self._perm
#        return self.__class__(self.parent(),
#                              tuple([-self._colors[i - 1] for i in ip]),  # -1 for indexing
#                              ip)
#
#    __invert__ = inverse
#
    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: s1*s2*s1 == s2*s1*s2
            True
            sage: t^4 == C.one()
            True
            sage: s1*s2 == s2*s1
            False
        """
        if not isinstance(other, ColoredPermutation):
            return False
        return (self.parent() is other.parent()
                and self._colors == other._colors
                and self._perm == other._perm)

    def __ne__(self, other):
        """
        Check inequality.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: s1*s2*s1 != s2*s1*s2
            False
            sage: s1*s2 != s2*s1
            True
        """
        return not self == other

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: x = s1*s2*t
            sage: list(x)
            [(1, 3), (0, 1), (0, 2)]
        """
        for i, p in enumerate(self._perm):
            yield (self._colors[i], p)

    def one_line_form(self):
        """
        Return the one line form of ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: x = s1*s2*t
            sage: x
            [[1, 0, 0], [3, 1, 2]]
            sage: x.one_line_form()
            [(1, 3), (0, 1), (0, 2)]
        """
        return list(self)

    def colors(self):
        """
        Return the colors of ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: x = s1*s2*t
            sage: x.colors()
            [1, 0, 0]
        """
        return list(self._colors)

    def permutation(self):
        """
        Return the permutation of ``self``.

        This is obtained by forgetting the colors.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: x = s1*s2*t
            sage: x.permutation()
            [3, 1, 2]
        """
        return self._perm

#    def to_matrix(self):
#        """
#        Return a matrix of ``self``.
#
#        The colors are mapped to roots of unity.
#
#        EXAMPLES::
#
#            sage: C = ColoredPermutations(4, 3)
#            sage: s1,s2,t = C.gens()
#            sage: x = s1*s2*t*s2; x.one_line_form()
#            [(1, 2), (0, 1), (0, 3)]
#            sage: M = x.to_matrix(); M
#            [    0     1     0]
#            [zeta4     0     0]
#            [    0     0     1]
#
#        The matrix multiplication is in the *opposite* order::
#
#            sage: M == s2.to_matrix()*t.to_matrix()*s2.to_matrix()*s1.to_matrix()
#            True
#        """
#        Cp = CyclotomicField(self.parent()._m)
#        g = Cp.gen()
#        D = diagonal_matrix(Cp, [g ** i for i in self._colors])
#        return self._perm.to_matrix() * D


class DecoratedPermutations(Parent, UniqueRepresentation):
    r"""
    The group of `m`-colored permutations on `\{1, 2, \ldots, n\}`.

    Let `S_n` be the symmetric group on `n` letters and `C_m` be the cyclic
    group of order `m`. The `m`-colored permutation group on `n` letters
    is given by `P_n^m = C_m \wr S_n`. This is also the complex reflection
    group `G(m, 1, n)`.

    We define our multiplication by

    .. MATH::

        ((s_1, \ldots s_n), \sigma) \cdot ((t_1, \ldots, t_n), \tau)
        = ((s_1 t_{\sigma(1)}, \ldots, s_n t_{\sigma(n)}), \tau \sigma).

    EXAMPLES::

        sage: C = ColoredPermutations(4, 3); C
        4-colored permutations of size 3
        sage: s1,s2,t = C.gens()
        sage: (s1, s2, t)
        ([[0, 0, 0], [2, 1, 3]], [[0, 0, 0], [1, 3, 2]], [[0, 0, 1], [1, 2, 3]])
        sage: s1*s2
        [[0, 0, 0], [3, 1, 2]]
        sage: s1*s2*s1 == s2*s1*s2
        True
        sage: t^4 == C.one()
        True
        sage: s2*t*s2
        [[0, 1, 0], [1, 2, 3]]

    We can also create a colored permutation by passing
    either a list of tuples consisting of ``(color, element)``::

        sage: x = C([(2,1), (3,3), (3,2)]); x
        [[2, 3, 3], [1, 3, 2]]

    or a list of colors and a permutation::

        sage: C([[3,3,1], [1,3,2]])
        [[3, 3, 1], [1, 3, 2]]

    There is also the natural lift from permutations::

        sage: P = Permutations(3)
        sage: C(P.an_element())
        [[0, 0, 0], [3, 1, 2]]

    REFERENCES:

    - :wikipedia:`Generalized_symmetric_group`
    - :wikipedia:`Complex_reflection_group`
    """
    def __init__(self, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: TestSuite(C).run()
            sage: C = ColoredPermutations(2, 3)
            sage: TestSuite(C).run()
            sage: C = ColoredPermutations(1, 3)
            sage: TestSuite(C).run()
        """
        self._n = ZZ(n)
#        self._C = IntegerModRing(self._m)
        self._P = Permutations(self._n)

#        if self._m == 1 or self._m == 2:
#            from sage.categories.finite_coxeter_groups import FiniteCoxeterGroups
#            category = FiniteCoxeterGroups().Irreducible()
#        else:
#            from sage.categories.complex_reflection_groups import ComplexReflectionGroups
#            category = ComplexReflectionGroups().Finite().Irreducible().WellGenerated()
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ColoredPermutations(4, 3)
            4-colored permutations of size 3
        """
        return "Decorated permutations of size {}".format(self._n)



    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        INPUT:

        Either a list of pairs (color, element)
        or a pair of lists (colors, elements).

        TESTS::

            sage: C = ColoredPermutations(4, 3)
            sage: x = C([(2,1), (3,3), (3,2)]); x
            [[2, 3, 3], [1, 3, 2]]
            sage: x == C([[2,3,3], [1,3,2]])
            True
        """
        if isinstance(x, list):
            if isinstance(x[0], tuple):
                c = []
                p = []
                for k in x:
                    if len(k) != 2:
                        raise ValueError("input must be pairs (color, element)")
                    if k[0]!=0 and k[0]!=1 and k[0]!=-1:
                        raise ValueError("colors can only be -1,0,1")                    
                    c.append(k[0])
                    p.append(k[1])
                    for i in p:
                        if i not in Permutation(p).fixed(point) and c[i]!=0:
                            raise ValueError("non-fixed points must be colored with 0")
                return self.element_class(self, c, self._P(p))

            if len(x) != 2:
                raise ValueError("input must be a pair of a list of colors and a permutation")
            return self.element_class(self, x[0], self._P(x[1]))

#Consider revisiting if there is need for a coercion from Permutations to decorated permutations (by choosing all fixed points to be colored up)

#    def _coerce_map_from_(self, C):
#        """
#        Return a coerce map from ``C`` if it exists and ``None`` otherwise.
#
#        EXAMPLES::
#
#            sage: C = ColoredPermutations(2, 3)
#            sage: S = SignedPermutations(3)
#            sage: C.has_coerce_map_from(S)
#            True
#
#            sage: C = ColoredPermutations(4, 3)
#            sage: C.has_coerce_map_from(S)
#            False
#            sage: S = SignedPermutations(4)
#            sage: C.has_coerce_map_from(S)
#            False
#
#            sage: P = Permutations(3)
#            sage: C.has_coerce_map_from(P)
#            True
#            sage: P = Permutations(4)
#            sage: C.has_coerce_map_from(P)
#            False
#        """
#        if isinstance(C, Permutations) and C.n == self._n:
#            return lambda P, x: P.element_class(P, [P._C.zero()]*P._n, x)
#        if self._m == 2 and isinstance(C, SignedPermutations) and C._n == self._n:
#            return lambda P, x: P.element_class(P,
#                                                [P._C.zero() if v == 1 else P._C.one()
#                                                 for v in x._colors],
#                                                x._perm)
#        return super(ColoredPermutations, self)._coerce_map_from_(C)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(2, 2)
            sage: [x for x in C]
            [[[0, 0], [1, 2]],
             [[0, 1], [1, 2]],
             [[1, 0], [1, 2]],
             [[1, 1], [1, 2]],
             [[0, 0], [2, 1]],
             [[0, 1], [2, 1]],
             [[1, 0], [2, 1]],
             [[1, 1], [2, 1]]]
        """
        for p in self._P:
            for c in itertools.product({-1,1}, repeat=p.number_of_fixed_points()):
                d=[0]*len(p)
                for i in p.fixed_points():
                    d[i-1]=c[p.fixed_points().index(i)]
                yield self.element_class(self, d, p)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.cardinality()
            384
            sage: C.cardinality() == 4**3 * factorial(3)
            True
        """
        return sum (factorial(self._n)/factorial(k) for k in range(self._n + 1))

    order = cardinality

    def rank(self):
        """
        Return the rank of ``self``.

        The rank of a complex reflection group is equal to the dimension
        of the complex vector space the group acts on.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 12)
            sage: C.rank()
            12
            sage: C = ColoredPermutations(7, 4)
            sage: C.rank()
            4
            sage: C = ColoredPermutations(1, 4)
            sage: C.rank()
            3
        """
        if self._m == 1:
            return self._n - 1
        return self._n

    def degrees(self):
        r"""
        Return the degrees of ``self``.

        The degrees of a complex reflection group are the degrees of
        the fundamental invariants of the ring of polynomial invariants.

        If `m = 1`, then we are in the special case of the symmetric group
        and the degrees are `(2, 3, \ldots, n, n+1)`. Otherwise the degrees
        are `(m, 2m, \ldots, nm)`.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.degrees()
            (4, 8, 12)
            sage: S = ColoredPermutations(1, 3)
            sage: S.degrees()
            (2, 3)

        We now check that the product of the degrees is equal to the
        cardinality of ``self``::

            sage: prod(C.degrees()) == C.cardinality()
            True
            sage: prod(S.degrees()) == S.cardinality()
            True
        """
        # For the usual symmetric group (self._m=1) we need to start at 2
        start = 2 if self._m == 1 else 1
        return tuple(self._m * i for i in range(start, self._n + 1))

    def codegrees(self):
        r"""
        Return the codegrees of ``self``.

        Let `G` be a complex reflection group. The codegrees
        `d_1^* \leq d_2^* \leq \cdots \leq d_{\ell}^*` of `G` can be
        defined by:

        .. MATH::

            \prod_{i=1}^{\ell} (q - d_i^* - 1)
            = \sum_{g \in G} \det(g) q^{\dim(V^g)},

        where `V` is the natural complex vector space that `G` acts on
        and `\ell` is the :meth:`rank`.

        If `m = 1`, then we are in the special case of the symmetric group
        and the codegrees are `(n-2, n-3, \ldots 1, 0)`. Otherwise the degrees
        are `((n-1)m, (n-2)m, \ldots, m, 0)`.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.codegrees()
            (8, 4, 0)
            sage: S = ColoredPermutations(1, 3)
            sage: S.codegrees()
            (1, 0)

        TESTS:

        We check the polynomial identity::

            sage: R.<q> = ZZ[]
            sage: C = ColoredPermutations(3, 2)
            sage: f = prod(q - ds - 1 for ds in C.codegrees())
            sage: d = lambda x: sum(1 for e in x.to_matrix().eigenvalues() if e == 1)
            sage: g = sum(det(x.to_matrix()) * q**d(x) for x in C)
            sage: f == g
            True
        """
        # Special case for the usual symmetric group
        last = self._n-1 if self._m == 1 else self._n
        return tuple(self._m * i for i in reversed(range(last)))

    def number_of_reflection_hyperplanes(self):
        """
        Return the number of reflection hyperplanes of ``self``.

        The number of reflection hyperplanes of a complex reflection
        group is equal to the sum of the codegrees plus the rank.

        EXAMPLES::

            sage: C = ColoredPermutations(1, 2)
            sage: C.number_of_reflection_hyperplanes()
            1
            sage: C = ColoredPermutations(1, 3)
            sage: C.number_of_reflection_hyperplanes()
            3
            sage: C = ColoredPermutations(4, 12)
            sage: C.number_of_reflection_hyperplanes()
            276
        """
        return sum(self.codegrees()) + self.rank()

    def fixed_point_polynomial(self, q=None):
        r"""
        The fixed point polynomial of ``self``.

        The fixed point polynomial `f_G` of a complex reflection group `G`
        is counting the dimensions of fixed points subspaces:

        .. MATH::

            f_G(q) = \sum_{w \in W} q^{\dim V^w}.

        Furthermore, let `d_1, d_2, \ldots, d_{\ell}` be the degrees of `G`,
        where `\ell` is the :meth:`rank`. Then the fixed point polynomial
        is given by

        .. MATH::

            f_G(q) = \prod_{i=1}^{\ell} (q + d_i - 1).

        INPUT:

        - ``q`` -- (default: the generator of ``ZZ['q']``) the parameter `q`

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.fixed_point_polynomial()
            q^3 + 21*q^2 + 131*q + 231

            sage: S = ColoredPermutations(1, 3)
            sage: S.fixed_point_polynomial()
            q^2 + 3*q + 2

        TESTS:

        We check the against the degrees and codegrees::

            sage: R.<q> = ZZ[]
            sage: C = ColoredPermutations(4, 3)
            sage: C.fixed_point_polynomial(q) == prod(q + d - 1 for d in C.degrees())
            True
        """
        if q is None:
            q = PolynomialRing(ZZ, 'q').gen(0)
        return prod(q + d - 1 for d in self.degrees())

    def is_well_generated(self):
        r"""
        Return if ``self`` is a well-generated complex reflection group.

        A complex reflection group `G` is well-generated if it is
        generated by `\ell` reflections. Equivalently, `G` is well-generated
        if `d_i + d_i^* = d_{\ell}` for all `1 \leq i \leq \ell`.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.is_well_generated()
            True
            sage: C = ColoredPermutations(2, 8)
            sage: C.is_well_generated()
            True
            sage: C = ColoredPermutations(1, 4)
            sage: C.is_well_generated()
            True
        """
        deg = self.degrees()
        dstar = self.codegrees()
        return all(deg[-1] == d + dstar[i] for i, d in enumerate(deg))

    Element = DecoratedPermutation
