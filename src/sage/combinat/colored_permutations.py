r"""
Colored Permutations

.. TODO::

    Much of the colored permutations (and element) class can be
    generalized to `G \wr S_n`
"""
from sage.categories.groups import Groups
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.finite_coxeter_groups import FiniteCoxeterGroups
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod

from sage.combinat.permutation import Permutations
from sage.combinat.cartesian_product import CartesianProduct
from sage.matrix.constructor import diagonal_matrix
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.all import ZZ


class ColoredPermutation(MultiplicativeGroupElement):
    """
    A colored permutation.
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
        MultiplicativeGroupElement.__init__(self, parent=parent)

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

    def _mul_(self, other):
        """
        Multiply ``self`` and ``other``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: s1*s2*s1 == s2*s1*s2
            True
        """
        colors = tuple(self._colors[i] + other._colors[val - 1]  # -1 for indexing
                       for i, val in enumerate(self._perm))
        p = self._perm._left_to_right_multiply_on_right(other._perm)
        return self.__class__(self.parent(), colors, p)

    def __invert__(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: ~t
            [[0, 0, 3], [1, 2, 3]]
            sage: all(x * ~x == C.one() for x in C.gens())
            True
        """
        ip = ~self._perm
        return self.__class__(self.parent(),
                              tuple([-self._colors[i - 1] for i in ip]),  # -1 for indexing
                              ip)

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
        return not self.__eq__(other)

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

    def to_matrix(self):
        """
        Return a matrix of ``self``.

        The colors are mapped to roots of unity.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: s1,s2,t = C.gens()
            sage: x = s1*s2*t*s2; x.one_line_form()
            [(1, 2), (0, 1), (0, 3)]
            sage: M = x.to_matrix(); M
            [    0     1     0]
            [zeta4     0     0]
            [    0     0     1]

        The matrix multiplication is in the *opposite* order::

            sage: M == s2.to_matrix()*t.to_matrix()*s2.to_matrix()*s1.to_matrix()
            True
        """
        Cp = CyclotomicField(self.parent()._m)
        g = Cp.gen()
        D = diagonal_matrix(Cp, [g ** i for i in self._colors])
        return self._perm.to_matrix() * D


# TODO: Parts of this should be put in the category of complex
# reflection groups
class ColoredPermutations(Parent, UniqueRepresentation):
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

    REFERENCES:

    - :wikipedia:`Generalized_symmetric_group`
    - :wikipedia:`Complex_reflection_group`
    """
    def __init__(self, m, n, category=None):
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
        if m <= 0:
            raise ValueError("m must be a positive integer")
        self._m = m
        self._n = n
        self._C = IntegerModRing(self._m)
        self._P = Permutations(self._n)
        if category is None:
            category = (Groups(), FiniteEnumeratedSets())
        Parent.__init__(self, category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ColoredPermutations(4, 3)
            4-colored permutations of size 3
        """
        return "{}-colored permutations of size {}".format(self._m, self._n)

    @cached_method
    def one(self):
        """
        Return the identity element of ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.one()
            [[0, 0, 0], [1, 2, 3]]
        """
        return self.element_class(self, [self._C.zero()] * self._n,
                                  self._P.identity())

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.gens()
            ([[0, 0, 0], [2, 1, 3]],
             [[0, 0, 0], [1, 3, 2]],
             [[0, 0, 1], [1, 2, 3]])
        """
        zero = [self._C.zero()] * self._n
        g = []
        for i in range(self._n - 1):
            p = range(1, self._n + 1)
            p[i] = i + 2
            p[i + 1] = i + 1
            g.append(self.element_class(self, zero, self._P(p)))
        zero[-1] = self._C.one()
        g.append(self.element_class(self, zero, self._P.identity()))
        return tuple(g)

    def matrix_group(self):
        """
        Return the matrix group corresponding to ``self``.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.matrix_group()
            Matrix group over Cyclotomic Field of order 4 and degree 2 with 3 generators (
            [0 1 0]  [1 0 0]  [    1     0     0]
            [1 0 0]  [0 0 1]  [    0     1     0]
            [0 0 1], [0 1 0], [    0     0 zeta4]
            )
        """
        from sage.groups.matrix_gps.finitely_generated import MatrixGroup
        return MatrixGroup([g.to_matrix() for g in self.gens()])

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        INPUT:

        either a list of pairs (color, element)
        or a pair of lists (colors, elements)

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
                    c.append(self._C(k[0]))
                    p.append(k[1])
                return self.element_class(self, c, self._P(p))

            if len(x) != 2:
                raise ValueError("input must be a pair of a list of colors and a permutation")
            return self.element_class(self, map(self._C, x[0]), self._P(x[1]))

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
        C = CartesianProduct(*[self._C] * self._n)
        for p in self._P:
            for c in C:
                yield self.element_class(self, c, p)

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
        return self._m ** self._n * self._P.cardinality()

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
        """
        Return the degrees of ``self``.

        The degrees of a complex reflection group are the degrees of
        the fundamental invariants of the ring of polynomial invariants.

        If `m = 1`, then we are in the special case of the symmetric group
        and the degrees are `(2, 3, \ldots, n, n+1)`. Otherwise the degrees
        are `(m, 2m, \ldots, nm)`.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.degrees()
            [4, 8, 12]
            sage: S = ColoredPermutations(1, 3)
            sage: S.degrees()
            [2, 3]

        We now check that the product of the degrees is equal to the
        cardinality of ``self``::

            sage: prod(C.degrees()) == C.cardinality()
            True
            sage: prod(S.degrees()) == S.cardinality()
            True
        """
        if self._m == 1:  # Special case for the usual symmetric group
            return range(2, self._n + 1)
        return [self._m * i for i in range(1, self._n + 1)]

    def codegrees(self):
        """
        Return the codegrees of ``self``.

        Let `G` be a complex reflection group. The codegrees
        `d_1^* \leq d_2^* \leq \cdots \leq d_{\ell}^*` of `G` can be
        defined by:

            \prod_{i=1}^{\ell} (q - d_i^* - 1)
            = \sum_{g \in G} \det(g) q^{\dim(V^g)},

        where `V` is the natural complex vector space that `G` acts on.

        If `m = 1`, then we are in the special case of the symmetric group
        and the codegrees are `(n-2, n-3, \ldots 1, 0)`. Otherwise the degrees
        are `((n-1)m, (n-2)m, \ldots, m, 0)`.

        EXAMPLES::

            sage: C = ColoredPermutations(4, 3)
            sage: C.codegrees()
            [8, 4, 0]
            sage: S = ColoredPermutations(1, 3)
            sage: S.codegrees()
            [1, 0]

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
        if self._m == 1:  # Special case for the usual symmetric group
            return list(reversed(range(self._n - 1)))
        return [self._m * i for i in reversed(range(self._n))]

    def number_reflection_hyperplanes(self):
        """
        Return the number of reflection hyperplanes of ``self``.

        The number of reflection hyperplanes of a complex reflection
        group is equal to the sum of the codegrees plus the rank.

        EXAMPLES::

            sage: C = ColoredPermutations(1, 2)
            sage: C.number_reflection_hyperplanes()
            1
            sage: C = ColoredPermutations(1, 3)
            sage: C.number_reflection_hyperplanes()
            3
            sage: C = ColoredPermutations(4, 12)
            sage: C.number_reflection_hyperplanes()
            276
        """
        return sum(self.codegrees()) + self.rank()

    def fixed_point_polynomial(self, q=None):
        r"""
        The fixed point polynomial of ``self``.

        The fixed point polynomial `f_G` of a complex reflection group `G` is
        counting the dimensions of fixed points subspaces:

        .. MATH::

            f_G(q) = \sum_{w \in W} q^{\dim V^w}.

        Furthermore, let `d_1, d_2, \ldots, d_{\ell}` be the degrees of `G`,
        then the fixed point polynomial is given by

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
        """
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

    Element = ColoredPermutation

#####################################################################
## Signed permutations


class SignedPermutation(ColoredPermutation):
    """
    A signed permutation.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: s1,s2,s3,s4 = S.gens()
            sage: s4*s1*s2*s3*s4
            [-4, 1, 2, -3]
        """
        return repr(list(self))

    _latex_ = _repr_

    def _mul_(self, other):
        """
        Multiply ``self`` and ``other``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: s1,s2,s3,s4 = S.gens()
            sage: x = s4*s1*s2*s3*s4; x
            [-4, 1, 2, -3]
            sage: x * x
            [3, -4, 1, -2]

        ::

            sage: s1*s2*s1 == s1*s2*s1
            True
            sage: s3*s4*s3*s4 == s4*s3*s4*s3
            True
        """
        colors = tuple(self._colors[i] * other._colors[val - 1]  # -1 for indexing
                       for i, val in enumerate(self._perm))
        p = self._perm._left_to_right_multiply_on_right(other._perm)
        return self.__class__(self.parent(), colors, p)

    def inverse(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: s1,s2,s3,s4 = S.gens()
            sage: x = s4*s1*s2*s3*s4
            sage: ~x
            [2, 3, -4, -1]
            sage: x * ~x == S.one()
            True
        """
        ip = ~self._perm
        return self.__class__(self.parent(),
                              tuple([self._colors[i - 1] for i in ip]),  # -1 for indexing
                              ip)

    __invert__ = inverse

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: s1,s2,s3,s4 = S.gens()
            sage: x = s4*s1*s2*s3*s4
            sage: [a for a in x]
            [-4, 1, 2, -3]
        """
        for i, p in enumerate(self._perm):
            yield self._colors[i] * p

    def to_matrix(self):
        """
        Return a matrix of ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: s1,s2,s3,s4 = S.gens()
            sage: x = s4*s1*s2*s3*s4
            sage: M = x.to_matrix(); M
            [ 0  1  0  0]
            [ 0  0  1  0]
            [ 0  0  0 -1]
            [-1  0  0  0]

        The matrix multiplication is in the *opposite* order::

            sage: m1,m2,m3,m4 = [g.to_matrix() for g in S.gens()]
            sage: M == m4 * m3 * m2 * m1 * m4
            True
        """
        return self._perm.to_matrix() * diagonal_matrix(self._colors)

    def has_left_descent(self, i):
        """
        Return ``True`` if ``i`` is a left descent of ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: s1,s2,s3,s4 = S.gens()
            sage: x = s4*s1*s2*s3*s4
            sage: [x.has_left_descent(i) for i in S.index_set()]
            [True, False, False, True]
        """
        n = self.parent()._n
        if i == n:
            return self._colors[-1] == -1
        if self._colors[i - 1] == -1:
            return self._colors[i] == 1 or self._perm[i - 1] < self._perm[i]
        return self._colors[i] == 1 and self._perm[i - 1] > self._perm[i]


class SignedPermutations(ColoredPermutations):
    r"""
    Group of signed permutations.

    The group of signed permutations is also known as the hyperoctahedral
    group, the Coxeter group of type `B_n`, and the 2-colored permutation
    group. Thus it can be constructed as the wreath product `S_2 \wr S_n`.

    EXAMPLES::

        sage: S = SignedPermutations(4)
        sage: s1,s2,s3,s4 = S.group_generators()
        sage: x = s4*s1*s2*s3*s4; x
        [-4, 1, 2, -3]
        sage: x^4 == S.one()
        True

    This is a finite Coxeter group of type `B_n`::

        sage: S.canonical_representation()
        Finite Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 4]
        [2 2 4 1]
        sage: S.long_element()
        [-4, -3, -2, -1]
        sage: S.long_element().reduced_word()
        [4, 3, 4, 2, 3, 4, 1, 2, 3, 4]

    REFERENCES:

    - :wikipedia:`Hyperoctahedral_group`
    """
    def __init__(self, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: TestSuite(S).run()
        """
        ColoredPermutations.__init__(self, 2, n, FiniteCoxeterGroups())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SignedPermutations(4)
            Signed permutations of 4
        """
        return "Signed permutations of {}".format(self._n)

    @cached_method
    def one(self):
        """
        Return the identity element of ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: S.one()
            [1, 2, 3, 4]
        """
        return self.element_class(self, [ZZ.one()] * self._n,
                                  self._P.identity())

    def simple_reflection(self, i):
        r"""
        Return the ``i``-th simple refection of ``self``.

        Return the generators of ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: S.simple_reflection(1)
            [2, 1, 3, 4]
            sage: S.simple_reflection(4)
            [1, 2, 3, -4]
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        if i < self._n:
            p = range(1, self._n + 1)
            p[i - 1] = i + 1
            p[i] = i
            return self.element_class(self, [ZZ.one()] * self._n, self._P(p))
        temp = [ZZ.one()] * self._n
        temp[-1] = -ZZ.one()
        return self.element_class(self, temp, self._P.identity())

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: S.gens()
            ([2, 1, 3, 4], [1, 3, 2, 4], [1, 2, 4, 3], [1, 2, 3, -4])
        """
        return tuple(self.simple_reflection(i) for i in self.index_set())

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        TESTS::

            sage: S = SignedPermutations(3)
            sage: x = S([(+1,1), (-1,3), (-1,2)]); x
            [1, -3, -2]
            sage: x == S([[+1,-1,-1], [1,3,2]])
            True
            sage: x == S([1, -3, -2])
            True
        """
        if isinstance(x, list):
            if isinstance(x[0], tuple):
                c = []
                p = []
                for k in x:
                    if len(k) != 2:
                        raise ValueError("input must be pairs (sign, element)")
                    if k[0] != 1 and k[0] != -1:
                        raise ValueError("the sign must be +1 or -1")
                    c.append(ZZ(k[0]))
                    p.append(k[1])
                return self.element_class(self, c, self._P(p))

            if len(x) == self._n:
                c = []
                p = []
                one = ZZ.one()
                for v in x:
                    if v > 0:
                        c.append(one)
                        p.append(v)
                    else:
                        c.append(-one)
                        p.append(-v)
                return self.element_class(self, c, self._P(p))

            if len(x) != 2:
                raise ValueError("input must be a pair of a list of signs and a permutation")
            if any(s != 1 and s != -1 for s in x[0]):
                raise ValueError("the sign must be +1 or -1")
            return self.element_class(self, map(ZZ, x[0]), self._P(x[1]))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(2)
            sage: [x for x in S]
            [[1, 2], [1, -2], [-1, 2], [-1, -2],
             [2, 1], [2, -1], [-2, 1], [-2, -1]]
        """
        one = ZZ.one()
        C = CartesianProduct(*[[one, -one]] * self._n)
        for p in self._P:
            for c in C:
                yield self.element_class(self, c, p)

    @cached_method
    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: S.index_set()
            (1, 2, 3, 4)
        """
        return tuple(range(1, self._n + 1))

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: S.coxeter_matrix()
            [1 3 2 2]
            [3 1 3 2]
            [2 3 1 4]
            [2 2 4 1]
        """
        from sage.combinat.root_system.cartan_type import CartanType
        return CartanType(['B', self._n]).coxeter_matrix()

    def long_element(self, index_set=None):
        """
        Return the longest element of ``self``, or of the
        parabolic subgroup corresponding to the given ``index_set``.

        INPUT:

        - ``index_set`` -- (optional) a subset (as a list or iterable)
          of the nodes of the indexing set

        EXAMPLES::

            sage: S = SignedPermutations(4)
            sage: S.long_element()
            [-4, -3, -2, -1]
        """
        if index_set is not None:
            return super(SignedPermutations, self).long_element()
        p = range(self._n, 0, -1)
        return self.element_class(self, [-ZZ.one()] * self._n, self._P(p))

    Element = SignedPermutation

# TODO: Make this a subgroup
#class EvenSignedPermutations(SignedPermutations):
#    """
#    Group of even signed permutations.
#    """
#    def _repr_(self):
#        """
#        Return a string representation of ``self``.
#        """
#        return "Even signed permtuations of {}".format(self._n)
#
#    def __iter__(self):
#        """
#        Iterate over ``self``.
#        """
#        for s in SignedPermutations.__iter__(self):
#            total = 0
#            for pm in s._colors:
#                if pm == -1:
#                    total += 1
#
#            if total % 2 == 0:
#                yield s
