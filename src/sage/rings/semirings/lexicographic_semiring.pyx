r"""
Lexiographic semirings

AUTHORS:

- Darij Grinberg (2014-02-14) - Initial version, based on work by
  Travis Scrimshaw
"""
#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>,
#                     2014 Darij Grinberg <darij@mit.edu>
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

include "../../ext/stdsage.pxi"

from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element cimport RingElement, Element, ModuleElement
from sage.categories.semirings import Semirings
from sage.categories.map cimport Map
from sage.sets.family import Family
from sage.rings.all import ZZ

import operator

cdef class LexicographicSemiringElement(RingElement):
    r"""
    An element in the lexicographic semiring of degree `n` over an ordered
    cancellative additive semigroup `R`. Either a tuple of `n` elements of
    `R`, or `\infty` (which is the zero of the semiring!). The operators
    `+` and `\cdot` are defined as the lexicographic operators `\oplus`
    and `\odot`, respectively (see :class:`LexicographicSemiring` for
    their definitions).
    """
    cdef tuple _val

    cdef LexicographicSemiringElement _new(self):
        """
        Return a new lexicographic semiring element with parent ``self`.
        """
        cdef LexicographicSemiringElement x
        x = PY_NEW(LexicographicSemiringElement)
        x._parent = self._parent
        x._val = self._val
        return x

    def __init__(self, parent, tuple val=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 5)
            sage: elt = T((3, 1, -2, 4, 5))
            sage: TestSuite(elt).run()
        """
        RingElement.__init__(self, parent)
        self._val = val

    def __reduce__(self):
        """
        Used in pickling lexicographic semiring elements.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 3)
            sage: elt = T((2, 0, 16))
            sage: elt.__reduce__()
            (<type 'sage.rings.semirings.lexicographic_semiring.LexicographicSemiringElement'>,
             (Lexicographic semiring of degree 3 over Rational Field, (2, 0, 16)))
        """
        return (LexicographicSemiringElement, (self.parent(), self._val))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 4)
            sage: T((1, 0, 3, -1))
            (1, 0, 3, -1)
            sage: T.infinity()
            +infinity
            sage: T = LexicographicSemiring(QQ, 5, False)
            sage: T.infinity()
            -infinity
        """
        if self._val is None:
            if self.parent()._use_min:
                return "+infinity"
            return "-infinity"
        return repr(self._val)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 2)
            sage: T((1, 6))
            (1, 6)
            sage: latex(T((1, -2)))
            \left(1, -2\right)
            sage: latex(T.infinity())
            \infty
            sage: T = LexicographicSemiring(QQ, 2, False)
            sage: latex(T.infinity())
            -\infty
            sage: T = LexicographicSemiring(QQ, 1)
            sage: latex(T((1,)))
            \left(1\right)
            sage: T = LexicographicSemiring(QQ, 0)
            sage: latex(T(()))
            \left(\right)
        """
        if self._val is None:
            if self.parent()._use_min:
                return "\\infty"
            return "-\\infty"
        res = "\\left("
        for i in self._val:
            res += i._latex_() + ", "
        if len(self._val) > 0:
            res = res[:-2]
        res += "\\right)"
        return res

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 4)
            sage: hash(T((2, 1, -2, 4))) == hash(T((2, 1, -2, 4)))
            True
            sage: hash(T.infinity()) == hash(T.infinity())
            True
        """
        return hash(self._val)

    # Comparisons

    def __richcmp__(left, right, int op):
        """
        Rich comparisons.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 3)
            sage: T((2, 1, -2)) == T((2, 1, -2))
            True
            sage: T.infinity() == T.infinity()
            True
            sage: T((2, 1, -2)) != T((2, 1, 3))
            True
            sage: T.infinity() != T.infinity()
            False
            sage: T((2, 1, -2)) < T((2, 1, 3))
            False
            sage: T((2, 2, -2)) < T((2, 1, 3))
            True
            sage: T((2, 1, -2)) < T((3, 1, 3))
            False
            sage: T((3, 2, 1)) < T((1, 2, 3))
            True
            sage: T((3, 2, 1)) < T.infinity()
            False
            sage: T.infinity() < T.infinity()
            False
            sage: T((1, 5, 1)) <= T((1, 6, 1))
            False
            sage: T.infinity() <= T.infinity()
            True
            sage: T((0, 1, 6)) > T((1, 0, 0))
            True
            sage: T((2, 2, 2)) > T.infinity()
            True
            sage: T.infinity() > T.infinity()
            False
            sage: T((1, 2, 5)) >= T((1, 0, 3))
            False
            sage: T((1, 2, 5)) >= T((1, 2, 1))
            False
            sage: T.infinity() >= T.infinity()
            True
        """
        return (<RingElement>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Return ``-1`` if ``left`` is less than ``right``, ``0`` if
        ``left`` and ``right`` are equal, and ``1`` if ``left`` is
        greater than ``right``.

        This uses the ordering which makes the lexicographic
        semiring into an ordered semiring. This is not the lexicographic
        ordering but its reverse!

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 3)
            sage: T((2, 1, 2)) == T((2, 1, 2))
            True
            sage: T((2, 0, 6)) != T((1, 2, 4))
            True
            sage: T((1, 4, 2)) != T((1, 3, 2))
            True
            sage: T((-1, 0, 1)) != T((-1, 0, -1))
            True
            sage: T((-1, 0, 1)) < T((-1/2, 1, 3))
            False
            sage: T((-1, 0, 1)) < T((-1, 1/2, 3))
            False
            sage: T((-1, 0, 1)) < T((-1, 0, 1/2))
            True
            sage: T((3/2, 1, 1/2)) > T((3/1, 1, 1))
            True
            sage: T.infinity() == T.infinity()
            True
            sage: T((1/2, 1, 1)) <= T.infinity()
            False
            sage: T((1, 2, 1)) > T.infinity()
            True

        Using the `\max` definition::

            sage: T = LexicographicSemiring(QQ, 4, False)
            sage: T((1, 3, 2)) == T((1, 3, 2))
            True
            sage: T((1, 2, 1)) != T((1, 3/2, 1))
            True
            sage: T((1, 3, 2)) < T((2, -10, 1/2))
            True
            sage: T((1, 2, 1)) > T((2, 0, 1))
            False
            sage: T.infinity() == T.infinity()
            True
            sage: T((2, 1, 1)) <= T.infinity()
            True
        """
        cdef LexicographicSemiringElement self, x
        self = left
        x = right

        if self._val is None:
            if x._val is None:
                return 0
            if self.parent()._use_min:
                return -1
            return 1

        if x._val is None:
            if self.parent()._use_min:
                return 1
            return -1

        if self.parent()._use_min:
            if self._val > x._val:
                return -1
            if self._val < x._val:
                return 1
        if self._val > x._val:
            return 1
        if self._val < x._val:
            return -1
        return 0

    cpdef ModuleElement _add_(left, ModuleElement right):
        """
        Add ``left`` to ``right``.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 2)
            sage: T((2, 1)) + T((3, 5))
            (2, 1)
            sage: T((-2, 3)) + T.infinity()
            (-2, 3)
            sage: T((2, -4)) + T((2, 5))
            (2, -4)
            sage: T.infinity() + T((1, 2))
            (1, 2)
            sage: T = LexicographicSemiring(QQ, 3, False)
            sage: T((2, 1, -2)) + T((2, 0, 3))
            (2, 1, -2)
            sage: T((1, 4, 9)) + T((2, 0, 0))
            (2, 0, 0)
            sage: T((1, 2, 3)) + T.infinity()
            (1, 2, 3)
            sage: T.infinity() + T((1, 0, 1))
            (1, 0, 1)
        """
        cdef LexicographicSemiringElement self, rhs
        self = left
        rhs = right
        if self._val is None:
            return rhs
        if rhs._val is None:
            return self
        cdef LexicographicSemiringElement x
        x = self._new()
        if self.parent()._use_min:
            x._val = min(self._val, rhs._val)
        else:
            x._val = max(self._val, rhs._val)
        return x

    def __neg__(self):
        """
        Return the additive inverse, which only exists for `\infty`.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 4)
            sage: -T.infinity()
            +infinity
            sage: -T((1, 2, 1, 2))
            Traceback (most recent call last):
            ...
            ArithmeticError: cannot negate any non-infinite element
        """
        if self._val is None:
            return self
        raise ArithmeticError("cannot negate any non-infinite element")

    cpdef RingElement _mul_(left, RingElement right):
        """
        Multiply ``left`` and ``right``.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 3)
            sage: T((2, 1, 0)) * T((-1, -1, -3))
            (1, 0, -3)
            sage: T((2, 1, 0)) * T.infinity()
            +infinity
            sage: T((-4, 1, 2)) * T((3, 2, 3))
            (-1, 3, 5)

        The min/max choice determining the addition has no effect on the
        multiplication::

            sage: T = LexicographicSemiring(QQ, 3, False)
            sage: T((2, 1, 0)) * T((-1, -1, -3))
            (1, 0, -3)
            sage: T((2, 1, 0)) * T.infinity()
            -infinity
            sage: T((-4, 1, 2)) * T((3, 2, 3))
            (-1, 3, 5)
        """
        cdef LexicographicSemiringElement self, rhs
        self = left
        rhs = right
        if self._val is None:
            return self
        if rhs._val is None:
            return rhs
        cdef LexicographicSemiringElement x
        x = self._new()
        # bugger cython
        # x._val = tuple(a + b for (a, b) in zip(self._val, rhs._val))
        x._val = tuple([a + b for (a, b) in zip(self._val, rhs._val)])
        return x

    cpdef RingElement _div_(left, RingElement right):
        """
        Divide ``left`` by ``right``.

        This is well-defined if the semigroup `R` is a group,
        and occasionally in other cases too.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 3)
            sage: T((2, 1/2, 5)) / T((1/2, 5, 1))
            (3/2, -9/2, 4)
            sage: T.infinity() / T((3, 0, -1))
            +infinity
        """
        cdef LexicographicSemiringElement self, rhs
        self = left
        rhs = right

        if rhs._val is None:
            raise ZeroDivisionError("Tropical division by infinity")
        if self._val is None:
            return self
        cdef LexicographicSemiringElement x
        x = self._new()
        # bugger cython
        # x._val = tuple(a - b for (a, b) in zip(self._val, rhs._val))
        x._val = tuple([a - b for (a, b) in zip(self._val, rhs._val)])
        return x

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        This is well-defined if the semigroup `R` is a group,
        and occasionally in other cases too.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 4)
            sage: ~T((2, 1, -1/2, 0))
            (-2, -1, 1/2, 0)
        """
        if self.is_one():
            return self
        if self._val is None:
            raise ZeroDivisionError("Tropical division by infinity")
        cdef LexicographicSemiringElement x
        x = self._new()
        x._val = tuple(-a for a in self._val)
        return x

    def __pow__(base, exp, dummy):
        """
        Return ``self`` to ``exp``.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 3)
            sage: elt = T((2, 1, -1/3))
            sage: elt**3
            (6, 3, -1)
            sage: elt**-2
            (-4, -2, 2/3)
            sage: elt**(3/7)
            (6/7, 3/7, -1/7)

            sage: elt = T.infinity()
            sage: elt**0
            (0, 0, 0)
            sage: elt**(1/2)
            +infinity
            sage: elt*33
            +infinity
        """
        cdef LexicographicSemiringElement self, x
        self = base
        if self._val is None:
            if exp > 0:
                return self
            elif exp == 0:
                return self.parent().one()
            raise ZeroDivisionError("Tropical division by infinity")
        x = self._new()
        x._val = tuple(exp*a for a in self._val)
        return x

    def multiplicative_order(self):
        """
        Return the multiplicative order of ``self``.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 3)
            sage: T.one().multiplicative_order()
            1
            sage: T.zero().multiplicative_order()
            +Infinity
        """
        if self.is_one():
            return ZZ.one()
        from sage.rings.infinity import infinity
        return infinity

class LexicographicSemiring(Parent, UniqueRepresentation):
    r"""
    The lexicographic semiring of a given degree over a given
    semigroup.

    Given a totally ordered cancellative additive (in particular, abelian)
    semigroup `R` and a nonnegative integer `n`, we define the
    *lexicographic semiring of degree `n` over `R`* as the set
    `R^n \cup \{+\infty\}` equipped with the following semiring structure:
    We order `R^n` lexicographically, and extend this to an ordering of
    `R^n \cup \{+\infty\}` by making `+\infty` bigger than everything
    else. We define the addition `\oplus` and the multiplication `\odot`
    on our lexicographic semiring by

    .. MATH::

        a \oplus b = \min\{a, b\}, \quad \quad a \odot b = a + b

    (where the minimum is with respect to the extended lexicographic
    order, and the sum `a + b` is componentwise). The element `+\infty`
    becomes the zero of this semiring, and the `n`-tuple
    `(0, 0, \ldots, 0)` formed out of the neutral element of the semigroup
    `R` becomes the unity of this semiring.

    We denote the lexicographic semiring of degree `n` over `R` by
    `\mbox{Lex}_n(R)`. We can turn it into a totally ordered commutative
    semiring by defining its order to be the *reverse* of the (extended)
    lexicographic order; all elements of `\mbox{Lex}_n(R)` are then
    `\geq \infty`. Note that no element of `\mbox{Lex}_n(R)` other than
    `\infty` has an additive inverse.

    When the semigroup `R` is a group, `\mbox{Lex}_n(R)` is a semifield.

    The ordered semiring `\mbox{Lex}_n(R)` has a distinguished `n`-tuple
    of elements, which we call `e_1, e_2, \ldots, e_n`, and which are
    defined by `e_i` being the `n`-tuple
    `(0, 0, \ldots, 0, 1, 0, 0, \ldots, 0)` with the `1` at the `i`-th
    position. When `R = \ZZ` (as additive group), these elements have
    the following universal property: If `B` is any totally ordered
    commutative semiring satisfying `b + c = \max\{b, c\}` for all
    `b, c \in B`, and if `b_1, b_2, \ldots, b_n` are `n` multiplicatively
    invertible elements of `B` satisfying `b_i < b_{i+1}^N` for all
    `i \in \{1, 2, \ldots, n-1\}` and `N \in \NN`, then there exists a
    unique semiring homomorphism `A \to B` sending each `e_i` to `b_i`.
    This semiring homomorphism is (automatically) order-preserving.
    (The intuition behind this is that computing in `\mbox{Lex}_n(R)`
    models working with piecewise linear functions in `n` variables,
    each of which is much larger than the next.)

    The semiring `\mbox{Lex}_n(R)` generalizes both the tropical
    semiring
    (:class:`~sage.rings.semirings.tropical_semiring.TropicalSemiring`,
    obtained for `n = 1`) and the Boolean semiring (obtained for
    `n = 0`).

    There is an alternative definition where we define
    `\mbox{Lex}'_n(R) = R \cup \{-\infty\}`
    and alter the addition to be defined by

    .. MATH::

        a \oplus b = \max\{a, b\}.

    To use the `\max` definition, set the argument ``use_min = False``.

    .. TODO::

        Check that everything is both implemented and documented
        correctly for this alternative definition!

    .. WARNING::

        :meth:`zero` and :meth:`one` refer to the additive and the
        multiplicative identities of `\mbox{Lex}_n(R)`. These are **not**
        usually `n`-tuples formed of zeroes resp. ones (even when such
        things are well-defined).

        Do not use ``sum(...)`` to sum an array of elements of a
        lexicographic semiring, as this tries to start summing from the
        integer `0` rather than from the zero of the semiring. Instead
        use the ``sum`` method of the lexicographic semiring::

            sage: T = LexicographicSemiring(QQ, 2)
            sage: T.sum([T((1, 2)), T((1, 1))])
            (1, 1)

    .. WARNING::

        Comparing elements of a lexicographic semiring is 

    .. TODO::

        Explain this and why.

    .. TODO::

        Ensure that comparison is sane with the opposite convention.
        (I think it is not.)

    INPUT:

    - ``base`` -- The base totally ordered cancellative additive semigroup
      `R`.
    - ``n`` -- the nonnegative integer `n`.
    - ``use_min`` -- (Default: ``True``) If ``True``, then the semiring uses
      `a \oplus b = \min(a, b)`; otherwise uses `a \oplus b = \max(a, b)`

    EXAMPLES::

        sage: T = LexicographicSemiring(QQ, 3)
        sage: elt = T((2, 3, 1)); elt
        (2, 3, 1)

    The addition in the lexicographic semiring is the lexicographic minimum
    of two elements::

        sage: T((-1, 2/3, 1)) + T((-1, 1/3, 1))
        (-1, 1/3, 1)

    The multiplication in the lexicographic semiring is coordinatewise
    addition::

        sage: T((1, 1, 2)) * T((-1, 2/3, 1))
        (0, 5/3, 3)
        sage: T((0, 1, 0)) * T((-2, 1, 3))
        (-2, 2, 3)

    We can also do division and exponentiation::

        sage: T((10, 11, 12)) / T((1, 3, 5))
        (9, 8, 7)
        sage: T((1, 3, 7))^(-3/7)
        (-3/7, -9/7, -3)

    Note that "zero" and "one" are the additive and multiplicative
    identities of the lexicographic semiring. They are not the `n`-tuples
    formed out of `0`s resp. `1`s (even when there is such a thing as `1`
    in `R`); instead they are the additive identity `+\infty` and the
    multiplicative identity `(0, 0, \ldots, 0)`, respectively, of the
    lexicographic semiring::

        sage: T.zero() + T((2, -1, 2)) == T((2, -1, 2))
        True
        sage: T.one() * T((2, -1, 2)) == T((2, -1, 2))
        True

    For `n = 1`, the lexicographic semiring is simply the tropical
    semiring (with a bulkier syntax)::

        sage: T = LexicographicSemiring(ZZ, 1)
        sage: T((2,))
        (2,)
        sage: T((2,)) + T((1,))
        (1,)
        sage: T((2,)) * T((-1,))
        (1,)

    For `n = 0`, we get the Boolean semiring::

        sage: T = LexicographicSemiring(ZZ, 0)
        sage: T(())
        ()
    """
    def __init__(self, base, n, use_min=True):
        r"""
        Initialize ``self``.

        TESTS::

            sage: T = LexicographicSemiring(QQ, 5); T
            Lexicographic semiring of degree 5 over Rational Field
            sage: TestSuite(T).run()
        """
        self._n = n
        self._use_min = use_min
        self._names = ('x', 'infty')
        Parent.__init__(self, base=base, n=n, category=Semirings())

    def n(self):
        r"""
        Return the nonnegative integer `n` such that ``self`` is a
        lexicographic semiring of degree `n` over a semigroup (which
        semigroup is returned by the :meth:`base` method).

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 2); T.n()
            2
        """
        return self._n

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LexicographicSemiring(QQ, 7)
            Lexicographic semiring of degree 7 over Rational Field
        """
        return "Lexicographic semiring of degree %s over %s"%(self.n(), self.base())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(LexicographicSemiring(QQ, 5))
            \mathrm{Lex}_{5}\left( \Bold{Q} \right)
            sage: latex(LexicographicSemiring(QQ, 5, False))
            \mathrm{Lex}^{\prime}_{5}\left( \Bold{Q} \right)
        """
        res = "\\mathrm{Lex}"
        if self._use_min is False:
            res += "^{\\prime}"
        res += "_{" + self._n._latex_() + "}\left( " + self.base()._latex_() + " \\right)"
        return res

    def _coerce_map_from_(self, S):
        """
        Canonical coercion of into ``self`` from ``S``.

        The only objects that canonically coerce to a lexicographic
        semiring are lexicographic semirings of the same ring whose base
        semigroups have a coercion.

        EXAMPLES::

            sage: TR = LexicographicSemiring(RR, 2)
            sage: T60 = LexicographicSemiring(RealField(60), 2)
            sage: TR.has_coerce_map_from(T60)
            True
            sage: TQ = LexicographicSemiring(QQ, 2)
            sage: TQ.has_coerce_map_from(LexicographicSemiring(ZZ, 2))
            True
            sage: TR.has_coerce_map_from(TR)
            True
            sage: TQ.has_coerce_map_from(TQ)
            True
            sage: TR.has_coerce_map_from(TQ)
            True
            sage: TR.has_coerce_map_from(float)
            False
            sage: TR.has_coerce_map_from(RR)
            False
            sage: TR.has_coerce_map_from(QQ)
            False
            sage: TR.coerce_map_from(T60)(T60((2, 3)))
            (2.000000000..., 3.0000000000...)
            sage: TR.coerce(T60((3.4, 1.5)))
            (3.400000000..., 1.5000000000...)
            sage: TR.coerce(T60.infinity())
            +infinity
            sage: TQ.coerce(TR((3.4, 1.22)))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Lexicographic semiring of degree 2 over
             Real Field with 53 bits of precision to Lexicographic semiring of degree 2
             over Rational Field
            sage: TQ3 = LexicographicSemiring(QQ, 3)
            sage: TQ3.has_coerce_map_from(TR)
            False
            sage: TQ3.has_coerce_map_from(TQ)
            False
            sage: TR.has_coerce_map_from(TQ3)
            False
            sage: TQ.has_coerce_map_from(TQ3)
            False
        """
        if isinstance(S, LexicographicSemiring) and self._use_min == S._use_min \
                and self._n == S._n and self.base().has_coerce_map_from(S.base()):
            return LexicographicToLexicographic(S, self)

    def _element_constructor_(self, val):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 3)
            sage: T((1, -1/5, 2))
            (1, -1/5, 2)
        """
        if val is not None:
            val = tuple(map(self.base(), val))
        return self.element_class(self, val)

    Element = LexicographicSemiringElement

    @cached_method
    def zero_element(self):
        """
        Return the additive identity element `+\infty` (or `-\infty`, if
        the alternative definition has been used) of the lexicographic
        semiring.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 5)
            sage: T.zero_element()
            +infinity
            sage: T = LexicographicSemiring(QQ, 5, False)
            sage: T.zero_element()
            -infinity
        """
        return self.element_class(self, None)

    zero = zero_element
    infinity = zero_element
    additive_identity = zero_element

    @cached_method
    def one_element(self):
        """
        Return the multiplicative identity element `(0, 0, \ldots, 0)` of
        the lexicographic semiring.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 4)
            sage: T.one_element()
            (0, 0, 0, 0)
        """
        base_zero = self.base().zero_element()
        return self.element_class(self, tuple(base_zero for _ in range(self._n)))

    one = one_element
    multiplicative_identity = one_element

    def _an_element_(self):
        """
        Return some element of ``self``.

        EXAMPLES::

            sage: T = LexicographicSemiring(ZZ, 3)
            sage: T.an_element()
            (0, 0, 0)
        """
        return self.one()

    def gens(self, with_infty=True):
        r"""
        Return the generators of ``self``.

        This makes sense only if `R` is a ring. In this case, the
        generators returned are the `n`-tuples
        `e_i = (0, 0, \ldots, 0, 1, 0, 0, \ldots, 0)` (with `1` in `i`-th
        position) and the infinity element `+\infty` (or `-\infty`
        depending on the construction of ``self``), in this order.
        If the optional keyword variable ``with_infty`` is set to
        ``False``, the infinity element is missing from the list.

        EXAMPLES::

            sage: T = LexicographicSemiring(QQ, 4)
            sage: T.gens()
            ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), +infinity)
            sage: T.gens(with_infty=False)
            ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))
        """
        x = []
        base = self.base()
        one = base.one()
        zero = base.zero()
        n = self._n
        for i in range(n):
            x.append(self.element_class(self, tuple([zero for _ in range(i)] + [one]
                                                    + [zero for _ in range(n - 1 - i)])))
        if with_infty:
            x.append(self.infinity())
        return tuple(x)

    def ith_monomial(self, i):
        r"""
        Return the ``i``-th monomial in ``self``.

        This makes sense only if `R` is a ring. In this case, the
        `i`-th monomial in the lexicographic semiring of degree `n` over
        `R` is defined as the element
        `e_i = (0, 0, \ldots, 0, 1, 0, 0, \ldots, 0)`, which is an
        `n`-tuple with `1` in `i`-th position and `0`'s everywhere else.

        EXAMPLES::

            sage: T = LexicographicSemiring(ZZ, 4)
            sage: T.ith_monomial(3)
            (0, 0, 1, 0)
            sage: T.ith_monomial(1)
            (1, 0, 0, 0)
        """
        n = self._n
        base = self.base()
        one = base.one()
        zero = base.zero()
        return self.element_class(self, tuple([zero for _ in range(i - 1)] + [one]
                                                    + [zero for _ in range(n - i)]))

cdef class LexicographicToLexicographic(Map):
    """
    Map from the Lexicographic semiring to itself (possibly with
    different bases, but with same `n`).
    Used in coercion.
    """
    cpdef LexicographicSemiringElement _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.semirings.lexicographic_semiring import LexicographicToLexicographic
            sage: TZ = LexicographicSemiring(ZZ, 2)
            sage: TQ = LexicographicSemiring(QQ, 2)
            sage: f = LexicographicToLexicographic(TZ, TQ)
            sage: a = TZ((1, 5))
            sage: f(a)
            (1, 5)
            sage: parent(f(a))
            Lexicographic semiring of degree 2 over Rational Field
            sage: f(TZ.infinity())
            +infinity
            sage: parent(f(TZ.infinity()))
            Lexicographic semiring of degree 2 over Rational Field
        """
        dom = self._domain
        cod = self._codomain
        cod_base = cod.base()
        if (<LexicographicSemiringElement>x)._val is None:
            return cod.infinity()
        return cod(tuple([cod_base(i) for i in (<LexicographicSemiringElement>x)._val]))

