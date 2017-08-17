r"""
Infinite Multivariate Polynomial Ring with Sparse Exponents

Examples
========

Creating an infinite polynomial ring is done by
::

    sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
    sage: P
    Infinite polynomial ring in x, y over Rational Field

Note that a ``'sparse'`` and a ``'dense'`` implementation are
available as well; see :mod:`sage.rings.polynomial.infinite_polynomial_ring`.

Using variables is easy::

    sage: x[0]
    x_0
    sage: y[42]
    y_42

Arithmetic operations are as usual::

    sage: x[1]^2 + y[2]^3*x[0]
    x_0*y_2^3 + x_1^2

Non-integer Indices
-------------------

Indices can be anything (hashable, totally orderable)::

    sage: from functools import total_ordering
    sage: @total_ordering
    ....: class a(object):
    ....:     def __init__(self, i, j):
    ....:         self.i = i
    ....:         self.j = j
    ....:     def __repr__(self):
    ....:         return '{i}_{j}'.format(i=self.i, j=self.j)
    ....:     def __key__(self):
    ....:         return (self.i + self.j, self.j)
    ....:     def __hash__(self):
    ....:         return hash(self.__key__())
    ....:     def __eq__(self, other):
    ....:         return self.__key__() == other.__key__()
    ....:     def __lt__(self, other):
    ....:         return self.__key__() < other.__key__()
    sage: P.<A> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
    sage: A[a(1,1)]
    A_1_1
    sage: A[a(1,1)] + A[a(2,0)] + A[a(0,2)]
    A_0_2 + A_1_1 + A_2_0
    sage: A[a(1,1)]^2 * A[a(2,1)] * A[a(1,2)]^3 + 3*A[a(4,4)]^7
    3*A_4_4^7 + A_1_1^2*A_2_1*A_1_2^3

::

    sage: class b(a):
    ....:     def __neg__(self):
    ....:         return self.__class__(-self.i, -self.j)
    sage: Q.<B> = InfinitePolynomialRing(QQ, order='degrevlex', implementation='sparse_exponents')
    sage: B[b(1,1)]
    B_1_1
    sage: B[b(1,1)] + B[b(2,0)] + B[b(0,2)]
    B_2_0 + B_1_1 + B_0_2
    sage: B[b(1,1)]^2 * B[b(2,1)] * B[b(1,2)]^3 + 3*B[b(4,4)]^7
    3*B_4_4^7 + B_1_1^2*B_2_1*B_1_2^3

Various
=======

AUTHORS:

- Daniel Krenn (2017)

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the
  Austrian Science Fund (FWF): P 24644-N26.

Classes and Methods
===================
"""
# *****************************************************************************
# Copyright (C) 2017 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

from six import iteritems
from six import itervalues
from six import string_types

from sage.misc.cachefunc import cached_method
from sage.rings.ring import Algebra
from sage.structure.element import CommutativeAlgebraElement
from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialGen as InfinitePolynomialGen_generic


def updated_by_adding_values(D, E):
    r"""
    Return a copy of dictionary ``D`` which is updated from
    dictionary/iterable ``E`` by adding values.

    INPUT:

    - ``D`` -- dictionary

    - ``E`` -- dictionary

    OUTPUT:

    dictionary

    EXAMPLES::

        sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import updated_by_adding_values
        sage: updated_by_adding_values({1: 2, 3: 4}, {2: 5, 3: 6})
        {1: 2, 2: 5, 3: 10}
    """
    from copy import copy
    DD = copy(D)
    for key, value in iteritems(E):
        try:
            DD[key] += value
            if not DD[key]:
                del DD[key]
        except KeyError:
            DD[key] = value
    return DD


def monomial_factory(data):
    r"""
    Convert ``data`` to :class:`Monomial`.

    INPUT:

    - ``data`` -- object

    OUTPUT:

    :class:`Monomial`

    EXAMPLES::

        sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import monomial_factory
        sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
        sage: mx = next(iter(x[0]._summands_))
        sage: monomial_factory(mx)
        x0_0
        sage: monomial_factory(x[13])
        x0_13
        sage: monomial_factory(({3: 4}, {5: 6, 7: 8}))
        x0_3^4*x1_5^6*x1_7^8

        sage: S.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse')
        sage: monomial_factory(x[13])
        x0_13
        sage: S.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='dense')
        sage: monomial_factory(y[31])
        x1_31
    """
    from .infinite_polynomial_element import InfinitePolynomial_sparse
    from .infinite_polynomial_element import InfinitePolynomial_dense

    if isinstance(data, Monomial):
        return data

    elif isinstance(data, InfinitePolynomial_sparse_exponents):
        summands = data._summands_
        if len(summands) != 1:
            raise ValueError('{} is not monomial'.format(data))
        monomial, coefficient = next(iteritems(summands))
        if coefficient != 1:
            raise ValueError('{} is not normalized monomial')
        return monomial

    elif isinstance(data, (InfinitePolynomial_sparse,
                           InfinitePolynomial_dense)):
        summands = data.dict()
        if len(summands) != 1:
            raise ValueError('{} is not monomial'.format(data))
        monomial, coefficient = next(iteritems(summands))
        if coefficient != 1:
            raise ValueError('{} is not normalized monomial')
        name_to_component = {name: i
                             for i, name in
                             enumerate(data.parent()._names)}
        indices = tuple((name_to_component[name], int(index))
            for name, index in
            (repr(gen).rsplit('_', 1)
             for gen in data.polynomial().parent().gens()))

        exponents = tuple({} for _ in data.parent()._names)
        for (component, index), exponent in zip(indices, monomial):
            exponents[component][index] = exponent
        return Monomial(exponents)

    else:
        return Monomial(data)


class Monomial(object):
    r"""
    Datastructure of a monomial in
    :class:`InfinitePolynomialRing_sparse_exponents`.
    """

    def __init__(self, exponents):
        r"""
        INPUT:

        - ``exponents`` -- tuple (or other iterable) of
          dictionaries mapping indices to exponents

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial
            sage: Monomial(({0: 2, 2: 3}, {0: 3, 3: 1}))
            x0_0^2*x0_2^3*x1_0^3*x1_3
        """
        self._exponents_ = tuple({index: int(exponent)
                                  for index, exponent in iteritems(component)
                                  if exponent != 0}
                                 for component in exponents)

    def __repr__(self, names=None):
        r"""
        Return a representation string of this monomial.

        INPUT:

        - ``names`` -- tuple (or other iterable) of strings

          If ``names`` is ``None``, then ``'x0'``, ``'x1'``, etc. are used.

        OUTPUT:

        string

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial
            sage: a = Monomial(({0: 2, 2: 3}, {0: 3, 3: 1}))
            sage: a  # indirect doctest
            x0_0^2*x0_2^3*x1_0^3*x1_3
            sage: a.__repr__(('x', 'y'))
            'x_0^2*x_2^3*y_0^3*y_3'
        """
        if names is None:
            names = tuple('x{}'.format(c) for c in range(len(self._exponents_)))
        return '*'.join('{}_{}'.format(name, index)
                        + ('^{}'.format(exponent) if exponent > 1 else '')
                        for component, name in zip(self._exponents_, names)
                        for index, exponent in sorted(iteritems(component)))

    __str__ = __repr__

    @cached_method
    def degree(self):
        r"""
        Return the degree of this monomial.

        OUTPUT:

        nonnegative integer

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial
            sage: a = Monomial(({0: 2, 2: 3}, {0: 3, 3: 1})); a
            x0_0^2*x0_2^3*x1_0^3*x1_3
            sage: a.degree()
            9
        """
        return sum(exponent
                   for component in self._exponents_
                   for exponent in itervalues(component))

    @cached_method
    def _sorting_key_lex_(self):
        r"""
        Return a key for sorting in lexicographic ordering.

        OUTPUT:

        tuple

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial

            sage: x_0 = Monomial(({0: 1},))
            sage: x_1 = Monomial(({1: 1},))
            sage: u_0 = Monomial(({0: 2},))
            sage: u_1 = Monomial(({1: 2},))
            sage: x_0._sorting_key_lex_() < x_1._sorting_key_lex_()
            True
            sage: u_0._sorting_key_lex_() < u_1._sorting_key_lex_()
            True
            sage: x_0._sorting_key_lex_() < u_1._sorting_key_lex_()
            True
            sage: u_0._sorting_key_lex_() < x_1._sorting_key_lex_()
            True
            sage: x_0._sorting_key_lex_() < u_0._sorting_key_lex_()
            True
            sage: x_1._sorting_key_lex_() < u_1._sorting_key_lex_()
            True

        ::

            sage: names = ('x', 'y')
            sage: x_0 = Monomial(({0: 1}, {}))
            sage: x_1 = Monomial(({1: 1}, {}))
            sage: y_0 = Monomial(({}, {0: 1}))
            sage: y_1 = Monomial(({}, {1: 1}))
            sage: x_0._sorting_key_lex_() < x_1._sorting_key_lex_()
            True
            sage: y_0._sorting_key_lex_() < y_1._sorting_key_lex_()
            True
            sage: y_0._sorting_key_lex_() < x_0._sorting_key_lex_()
            True
            sage: y_1._sorting_key_lex_() < x_0._sorting_key_lex_()
            True
            sage: y_0._sorting_key_lex_() < x_1._sorting_key_lex_()
            True
            sage: y_1._sorting_key_lex_() < x_1._sorting_key_lex_()
            True

            sage: sorted([x_0, x_1, y_0, y_1], key=lambda k: k._sorting_key_lex_())
            [x1_0, x1_1, x0_0, x0_1]

            sage: x_2 = Monomial(({2: 1}, {}))
            sage: x_3 = Monomial(({3: 1}, {}))
            sage: y_2 = Monomial(({}, {2: 1}))
            sage: y_3 = Monomial(({}, {3: 1}))
            sage: sorted([x_2, x_3, y_0, y_1], key=lambda k: k._sorting_key_lex_())
            [x1_0, x1_1, x0_2, x0_3]
            sage: sorted([x_0, x_1, y_2, y_3], key=lambda k: k._sorting_key_lex_())
            [x1_2, x1_3, x0_0, x0_1]

            sage: sorted([x_0, x_1, x_2, x_3, y_0, y_1, y_2, y_3],
            ....:        key=lambda k: k._sorting_key_lex_())
            [x1_0, x1_1, x1_2, x1_3, x0_0, x0_1, x0_2, x0_3]

        ::

            sage: a = Monomial(({0: 2, 2: 3}, {0: 3, 3: 1})); a.__repr__(names)
            'x_0^2*x_2^3*y_0^3*y_3'
            sage: b = Monomial(({0: 2, 2: 3}, {0: 2, 3: 2})); b.__repr__(names)
            'x_0^2*x_2^3*y_0^2*y_3^2'
            sage: c = Monomial(({0: 3, 2: 2}, {0: 2, 3: 2})); c.__repr__(names)
            'x_0^3*x_2^2*y_0^2*y_3^2'
            sage: d = Monomial(({0: 1, 2: 4}, {0: 2, 3: 2})); d.__repr__(names)
            'x_0*x_2^4*y_0^2*y_3^2'
            sage: a._sorting_key_lex_() < b._sorting_key_lex_()
            True
            sage: c._sorting_key_lex_() < a._sorting_key_lex_()
            True
            sage: a._sorting_key_lex_() < d._sorting_key_lex_()
            True
            sage: c._sorting_key_lex_() < b._sorting_key_lex_()
            True
            sage: b._sorting_key_lex_() < d._sorting_key_lex_()
            True
            sage: c._sorting_key_lex_() < d._sorting_key_lex_()
            True
            sage: e = Monomial(({0: 4, 2: 2}, {0: 2, 3: 2})); e.__repr__(names)
            'x_0^4*x_2^2*y_0^2*y_3^2'
            sage: c._sorting_key_lex_() < e._sorting_key_lex_()
            True
        """
        return tuple(sum(sorted(iteritems(component), reverse=True), ())
                     for component in self._exponents_)

    def _sorting_key_revlex_(self):
        r"""
        Return a key for sorting in reverse lexicographic ordering.

        OUTPUT:

        tuple

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial

            sage: x_0 = Monomial(({0: 1},))
            sage: x_0._sorting_key_lex_()
            ((0, 1),)
            sage: x_0._sorting_key_revlex_()
            ((0, -1),)
        """
        return tuple(tuple(-t for t in T) for T in self._sorting_key_lex_())

    def _sorting_key_deglex_(self):
        r"""
        Return a key for sorting in degree lexicographic ordering.

        OUTPUT:

        tuple

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial

            sage: x_0 = Monomial(({0: 1},))
            sage: x_0._sorting_key_deglex_()
            (1, ((0, 1),))
        """
        return (self.degree(), self._sorting_key_lex_())

    def _sorting_key_degrevlex_(self):
        r"""
        Return a key for sorting in degree reverse lexicographic ordering.

        OUTPUT:

        tuple

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='degrevlex', implementation='sparse_exponents')
            sage: x[0] + y[1]
            x_0 + y_1
        """
        return (self.degree(), self._sorting_key_revlex_())

    def __hash__(self):
        r"""
        Return a hash value of this monomial.

        OUTPUT:

        integer

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: mx = next(iter(x[0]._summands_))
            sage: hash(mx)  # random
            -42
            sage: my = next(iter(y[0]._summands_))
            sage: hash(mx) == hash(my)
            False
            sage: mx2 = Monomial(({0: 1}, {}))
            sage: mx is mx2, mx == mx2, hash(mx) == hash(mx2)
            (False, True, True)
        """
        return hash(self._sorting_key_lex_())

    def __eq__(self, other):
        r"""
        Return whether this monomial equals the monomial ``other``.

        INPUT:

        - ``other`` -- :class:`Monomial`

        OUTPUT:

        A boolean.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: mx = next(iter(x[0]._summands_))
            sage: mx == mx
            True
            sage: my = next(iter(y[0]._summands_))
            sage: mx == my
            False
        """
        return self._sorting_key_lex_() == other._sorting_key_lex_()

    def __ne__(self, other):
        r"""
        Return whether this monomial does not equal the monomial ``other``.

        INPUT:

        - ``other`` -- :class:`Monomial`

        OUTPUT:

        A boolean.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: mx = next(iter(x[0]._summands_))
            sage: mx != mx
            False
            sage: my = next(iter(y[0]._summands_))
            sage: mx != my
            True
        """
        return not self.__eq__(other)

    def __mul__(self, other):
        r"""
        Return the product of this monomial with the monomial ``other``.

        INPUT:

        - ``other`` -- :class:`Monomial`

        OUTPUT:

        :class:`Monomial`

        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial
            sage: x_0 = Monomial(({0: 1}, {}))
            sage: x_1 = Monomial(({1: 1}, {}))
            sage: y_0 = Monomial(({}, {0: 1}))
            sage: y_1 = Monomial(({}, {1: 1}))
            sage: x_0*x_0
            x0_0^2
            sage: x_0*x_1
            x0_0*x0_1
            sage: x_0*y_0
            x0_0*x1_0
            sage: x_0*y_1
            x0_0*x1_1
        """
        return Monomial(
            tuple(updated_by_adding_values(scomponent, ocomponent)
                  for scomponent, ocomponent
                  in zip(self._exponents_, other._exponents_)))

    def subs(self, rules):
        r"""
        Substitute according to the given rules.

        This is called by :meth:`InfinitePolynomial_sparse_exponents.subs`.

        INPUT:

        - ``rules`` -- a tuple (or other iterable) of the same length
          as the number of generators with indexable
          entries

        OUTPUT:

        object

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial
            sage: a = Monomial(({0: 5, 1: 6}, {2: 7}))
            sage: print(a.__repr__(('x', 'y')))
            x_0^5*x_1^6*y_2^7
            sage: a.subs([x, y])
            x_0^5*x_1^6*y_2^7
        """
        from sage.misc.misc_c import prod
        return prod(rule[index]**exponent
                    for rule, component in zip(rules, self._exponents_)
                    for index, exponent in iteritems(component))

    def indices_of_variables(self):
        r"""
        Return the indices of the variables occurring in this monomial.

        OUTPUT:

        tuple of sets of indices

        EXAMPLES::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial
            sage: a = Monomial(({0: 5, 1: 6}, {2: 7}))
            sage: a.indices_of_variables()
            ({0, 1}, {2})
        """
        return tuple(set(component)
                     for component in self._exponents_)

    def indices_nonempty_components(self):
        r"""
        Return the indices of the components occurring in this monomial.

        OUTPUT:

        tuple

        EXAMPLES::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparse_exponents import Monomial
            sage: a = Monomial(({0: 5, 1: 6}, {2: 7}))
            sage: a.indices_nonempty_components()
            (0, 1)
        """
        return tuple(c
                     for c, component in enumerate(self._exponents_)
                     if component)


class InfinitePolynomial_sparse_exponents(CommutativeAlgebraElement):
    r"""
    A polynomial with sparse exponents of an infinite polynomial ring.
    """

    def __init__(self, parent, data):
        r"""
        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: x[0]  # indirect doctest
            x_0
        """
        super(InfinitePolynomial_sparse_exponents, self).__init__(parent=parent)

        coefficient_ring = parent.coefficient_ring()
        if isinstance(data, Monomial):
            self._summands_ = {data: coefficient_ring(1)}
        elif isinstance(data, dict):
            self._summands_ = {monomial_factory(monomial):
                               coefficient_ring(coefficient)
                               for monomial, coefficient in iteritems(data)
                               if coefficient != 0}
        else:
            raise TypeError('cannot create polynomial out of {}'.format(data))

    def __iter__(self):
        r"""
        Return an iterator over all pairs ``(coefficient, monomial)``
        of this polynomial.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: a = x[0] + 2*x[1] + y[0]*y[1]
            sage: sorted(a, key=lambda cm: repr(cm[1]))
            [(1, x_0), (2, x_1), (1, y_0*y_1)]
        """
        parent = self.parent()
        cls = parent.element_class
        return iter((coefficient, cls(parent, monomial))
                    for monomial, coefficient in iteritems(self._summands_))

    @cached_method
    def _sorted_monomials_and_coefficients_(self, reverse=False):
        r"""
        Return a sorted list over all pairs ``(monomial, coefficient)``
        of this polynomial.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: a = x[0] + 2*x[1] + y[0]*y[1]
            sage: a._sorted_monomials_and_coefficients_()
            [(x0_0, 1), (x0_1, 2), (x1_0*x1_1, 1)]
        """
        def key(monomial_and_coefficient):
            monomial = monomial_and_coefficient[0]
            return self.parent()._sorting_key_monomial_(monomial)
        return sorted(iteritems(self._summands_), key=key, reverse=reverse)

    def _repr_(self):
        r"""
        Return the representation string of this polynomial.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: 2*x[0]
            2*x_0
            sage: x[0] - x[0]
            0
            sage: P(1)
            1
            sage: P(42)
            42
            sage: 1 + y[0]
            y_0 + 1
            sage: y[42]-2*x[13]
            -2*x_13 + y_42
            sage: x[1] - y[2]
            x_1 - y_2
            sage: y[1] - x[2]
            -x_2 + y_1

        ::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='degrevlex', implementation='sparse_exponents')
            sage: x[0] + y[1]
            x_0 + y_1
        """
        def summand(monomial, coefficient):
            if coefficient == 1:
                sc = ''
                sign = ''
            elif coefficient == -1:
                sc = ''
                sign = '-'
            else:
                sc = str(coefficient)
                sign = ''
            factors = (sc, monomial.__repr__(self.parent()._names_))
            s = sign + '*'.join(f for f in factors if f)
            return s or '1'

        r = ' + '.join(summand(monomial, coefficient)
                       for monomial, coefficient
                       in self._sorted_monomials_and_coefficients_(reverse=True))
        r = r.replace(' + -', ' - ')
        return r or '0'

    def __hash__(self):
        r"""
        Return a hash value of this polynomial.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='degrevlex', implementation='sparse_exponents')
            sage: hash(x[0] + y[1])  # random output
            42
        """
        return hash((self.parent(),
                     tuple(self._sorted_monomials_and_coefficients_(reverse=True))))

    def __bool__(self):
        r"""
        Return whether this polynomial is nonzero.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='degrevlex', implementation='sparse_exponents')
            sage: bool(P(0))
            False
            sage: bool(x[1])
            True
        """
        return bool(self._summands_)

    __nonzero__ = __bool__

    def __eq__(self, other):
        r"""
        Return whether this polynomial is equal to ``other``.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: x[0] == y[0], x[0] != y[0]
            (False, True)
            sage: P(1) == 1, P(1) != 1
            (True, False)

        ::

            sage: P(1) == P(0)
            False
            sage: P(0) == P(0)
            True
        """
        if other is None:
            return False
        if (isinstance(other, InfinitePolynomial_sparse_exponents)
            and not other._summands_):
                return not bool(self)
        try:
            return not bool(self - other)
        except (TypeError, ValueError):
            return False

    def __ne__(self, other):
        r"""
        Return whether this polynomial is not equal to ``other``.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: x[0] == y[0], x[0] != y[0]
            (False, True)
            sage: P(1) == 1, P(1) != 1
            (True, False)
        """
        return not self.__eq__(other)

    def _add_(self, other):
        r"""
        Return the sum of this polynomial and ``other``.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: x[0] + y[0]
            x_0 + y_0
            sage: x[0] + y[0] + x[0]
            2*x_0 + y_0

        TESTS::

            sage: x[0] + 1
            x_0 + 1
        """
        summands = updated_by_adding_values(
            self._summands_, other._summands_)
        return self.parent().element_class(self.parent(), summands)

    def _lmul_(self, other):
        r"""
        Return the product of this polynomial and the scalar ``other``.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: (42/37)*x[0]
            42/37*x_0
        """
        if not other:
            return self.parent().zero()
        summands = {monomial: other*coefficient
                    for monomial, coefficient in iteritems(self._summands_)}
        return self.parent().element_class(self.parent(), summands)

    def _sub_(self, other):
        r"""
        Return the difference of this polynomial and ``other``.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: x[0] - x[0]
            0
            sage: x[1] - y[2]
            x_1 - y_2
            sage: y[1] - x[2]
            -x_2 + y_1
        """
        return self + self.parent().coefficient_ring()(-1) * other

    def _mul_by_summand_(self, other_monomial, other_coefficient):
        r"""
        Return the product of this polynomial and the given data.

        This method is called in :meth:`_mul_`.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: x[0]*y[0]  # indirect doctest
            x_0*y_0
        """
        summands = {monomial*other_monomial: coefficient*other_coefficient
                    for monomial, coefficient in iteritems(self._summands_)}
        return self.parent().element_class(self.parent(), summands)

    def _mul_(self, other):
        r"""
        Return the product of this polynomial and ``other``.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: x[0]*y[0]
            x_0*y_0
        """
        return sum((self._mul_by_summand_(monomial, coefficient)
                    for monomial, coefficient in iteritems(other._summands_)),
                   self.parent().zero())

    def degree(self):
        r"""
        Return the total degree of this polynomial.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: x[0].degree()
            1
            sage: (x[0]*x[2]^3).degree()
            4
            sage: (x[1]^4*y[2]^5 + x[3]^5*y[3]^6).degree()
            11

        TESTS::

            sage: P(1/2).degree()
            0
            sage: P(0).degree()
            -1
        """
        if not self._summands_:
            return -1
        return max(monomial.degree() for monomial in self._summands_)

    def subs(self, **rules):
        r"""
        Substitute the specified generators.

        INPUT:

        - ``rules`` -- keyword arguments ``generator=replacement``

          A generator can be indexed (e.g. ``x_3``) or not (e.g. ``x``).

          .. NOTE::

              A more specific rules overrides a more general rules.
              For example, if ``x=a`` as well as ``x_3=b`` are both specified,
              then ``x_3`` is replaced according to its rule by ``b``
              and ``x_i`` for `i\neq3` are replaced according to the
              rule of ``x``, namely by ``a[i]``.

        OUTPUT:

        object

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: Q.<a, b> = InfinitePolynomialRing(QQ)
            sage: (2*x[4] + 3*y[5]).subs(x=b, y=a)
            3*a_5 + 2*b_4

            sage: (2*x[4] + 3*x[5]).subs(x_5=15)
            2*x_4 + 45
            sage: (2*x[4] + 3*x[5]).subs(x=y, x_5=x[5])
            3*x_5 + 2*y_4
        """
        class Rule(object):
            def __init__(self, default):
                self.default = default
                self.rules = {}
            def __getitem__(self, index):
                try:
                    return self.rules[index]
                except KeyError:
                    return self.default[index]

        rules_monomial = tuple(Rule(gen) for gen in self.parent().gens())
        for name, replacement in iteritems(rules):
            try:
                rule = rules_monomial[self.parent()._index_by_name_(name)]
                rule.default = replacement
                continue
            except ValueError:
                pass

            gen, index = name.rsplit('_', 1)
            try:
                rule = rules_monomial[self.parent()._index_by_name_(gen)]
                rule.rules[int(index)] = replacement
                continue
            except ValueError:
                pass

            raise ValueError("'{}' does not specify a generator of {}".format(
                name, self.parent()))

        return sum(coefficient * monomial.subs(rules_monomial)
                   for monomial, coefficient in iteritems(self._summands_))

    def variables(self, skip_indices=False):
        r"""
        Return the variables occurring in this polynomial.

        INPUT:

        - ``skip_indices`` -- boolean

          If set, then the indices are skipped, i.e. only the generators
          occurring in this polynomial are returned.

        OUTPUT:

        A tuple of generators/variables.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')

            sage: (x[2]*y[3] + x[4]).variables()
            (x_2, x_4, y_3)
            sage: (x[2]*y[3] + x[4]).variables(skip_indices=True)
            (x_*, y_*)

            sage: (x[2] + x[4]).variables()
            (x_2, x_4)
            sage: (x[2] + x[4]).variables(skip_indices=True)
            (x_*,)

        TESTS::

            sage: P(1).variables()
            ()
            sage: P(1).variables(skip_indices=True)
            ()
            sage: P(0).variables()
            ()
            sage: P(0).variables(skip_indices=True)
            ()
        """
        if skip_indices:
            indices = set()
            for monomial, _ in iteritems(self._summands_):
                indices.update(monomial.indices_nonempty_components())
            return tuple(g
                         for i, g in enumerate(self.parent().gens())
                         if i in indices)
        else:
            indices = tuple(set() for _ in self.parent()._names_)
            for monomial, _ in iteritems(self._summands_):
                for I, J in zip(indices, monomial.indices_of_variables()):
                    I.update(J)
            return tuple(gen[i]
                         for gen, I in zip(self.parent().gens(), indices)
                         for i in sorted(I))


class InfinitePolynomialGen_sparse_exponents(InfinitePolynomialGen_generic):
    r"""
    A infinite generator for infinite polynomial rings with sparse exponents.
    """
    def __init__(self, parent, name, index):
        r"""
        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P.gen()  # indirect doctest
            x_*
        """
        self._index = index
        super(InfinitePolynomialGen_sparse_exponents, self).__init__(parent, name)

    def __hash__(self):
        r"""
        Return a hash value of this infinite polynomial generator.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: hash(x)  # random output
            42
        """
        return hash((self._parent, self._name, self._index))

    @cached_method
    def __getitem__(self, i):
        r"""
        Return the ``i``th variable associated to this infinite generator.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: x[0]  # indirect doctest
            x_0
        """
        monomial = Monomial(
            tuple(({i: 1} if index == self._index else {})
                  for index in range(len(self._parent._names_))))
        return InfinitePolynomial_sparse_exponents(self._parent, monomial)


class InfinitePolynomialRing_sparse_exponents(Algebra, UniqueRepresentation):

    Element = InfinitePolynomial_sparse_exponents

    @staticmethod
    def __classcall__(cls, coefficient_ring, names,
                      order='lex',
                      category=None):
        r"""
        Create a unique key.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P
            Infinite polynomial ring in x, y over Rational Field

        ::

            sage: InfinitePolynomialRing(QQ, names=('x', 'y'), order='deglex', implementation='sparse_exponents') is P
            True
            sage: InfinitePolynomialRing(QQ, ('x', 'y'), order='deglex', implementation='sparse_exponents') is P
            True
        """

        names = tuple(names)

        if category is None:
            from sage.categories.commutative_algebras import CommutativeAlgebras
            from sage.categories.rings import Rings
            category = CommutativeAlgebras(Rings())

        return super(InfinitePolynomialRing_sparse_exponents,
                     cls).__classcall__(cls, coefficient_ring, names,
                                        order=order,
                                        category=category)

    def __init__(self, coefficient_ring, names, order, category):
        r"""
        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P
            Infinite polynomial ring in x, y over Rational Field

        ::

            sage: P.category()
            Category of commutative algebras over rings
        """
        from sage.categories.rings import Rings
        from sage.rings.polynomial.term_order import TermOrder
        from sage.symbolic.ring import isidentifier

        if coefficient_ring not in Rings():
            raise ValueError('%s is not a ring. Cannot continue.' % (coefficient_ring,))
        self._coefficient_ring_ = coefficient_ring

        if order not in ('lex', 'deglex', 'degrevlex'):
            raise ValueError("wrong order '{}'".format(order))
        self._order_ = TermOrder(order)

        if not isinstance(names, tuple):
            raise TypeError('wrong names {}'.format(names))
        if not all(isinstance(name, string_types) and isidentifier(name)
                   for name in names):
            raise ValueError('wrong names {}'.format(names))
        self._names_ = names

        super(InfinitePolynomialRing_sparse_exponents, self).__init__(
            base_ring=coefficient_ring,
            category=category)

    def _repr_(self):
        r"""
        Return the representation string of this infinite polynomial ring.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P  # indirect doctest
            Infinite polynomial ring in x, y over Rational Field
        """
        return 'Infinite polynomial ring in {} over {}'.format(
            ', '.join(self._names_),
            self.coefficient_ring())

    @cached_method
    def __hash__(self):
        r"""
        Return a hash value of this infinite polynomial ring.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: hash(P)  # random output
            42
        """
        return hash((self.coefficient_ring(),
                     self._names_,
                     self.term_order()))

    def coefficient_ring(self):
        r"""
        Return the coefficient ring of this infinite polynomial ring.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P.coefficient_ring()
            Rational Field
        """
        return self._coefficient_ring_

    def term_order(self):
        r"""
        Return the term order of this infinite polynomial ring.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P.term_order()
            Degree lexicographic term order
        """
        return self._order_

    @cached_method
    def gens(self):
        r"""
        Return the infinite (indexable) generators of this
        infinite polynomial ring.

        OUTPUT:

        tuple of :class:`InfinitePolynomialRing_sparse_exponents`

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')  # indirect doctest
            sage: P.gens()
            (x_*, y_*)
            sage: x
            x_*
            sage: x[3]
            x_3
            sage: y
            y_*
            sage: y[42]
            y_42
        """
        return tuple(InfinitePolynomialGen_sparse_exponents(self, name, index)
                     for index, name in enumerate(self._names_))

    def gen(self, n=0):
        r"""
        Return the `n`th infinite (indexable) generators of this
        infinite polynomial ring.

        OUTPUT:

        :class:`InfinitePolynomialRing_sparse_exponents`

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P.gen()
            x_*
            sage: P.gen(1)
            y_*
        """
        return self.gens()[n]

    def _index_by_name_(self, name):
        r"""
        Return the index of the generator ``name`` in all generators.

        INPUT:

        - ``name`` -- string

        OUTPUT:

        integer

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P._index_by_name_('y')
            1
        """
        try:
            return self._names_.index(name)
        except ValueError:
            raise ValueError("'{}' does not specify a generator of {}".format(
                name, self))

    def gen_by_name(self, name):
        r"""
        Return the generator with given ``name``.

        INPUT:

        - ``name`` -- string

        OUTPUT:

        :class:`InfinitePolynomialRing_sparse_exponents`

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P.gen_by_name('y')
            y_*
        """
        return self.gen(self._index_by_name_(name))

    def ngens(self):
        r"""
        Return the number of infinite (indexable) generators of this
        infinite polynomial ring.

        OUTPUT:

        integer

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P.ngens()
            2
        """
        return len(self.gens())

    def _sorting_key_monomial_(self, monomial):
        r"""
        Return a key for sorting monomials according to the rings order.

        INPUT:

        - ``monomial`` -- :class:`Monomial`

        OUTPUT:

        tuple

        TESTS::

            sage: P.<x> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: monomial = next(iter(x[0]._summands_))
            sage: P._sorting_key_monomial_(monomial)
            (1, ((0, 1),))

        ::

            sage: P.<x> = InfinitePolynomialRing(QQ, order='degrevlex', implementation='sparse_exponents')
            sage: monomial = next(iter(x[0]._summands_))
            sage: P._sorting_key_monomial_(monomial)
            (1, ((0, -1),))
        """
        return getattr(monomial,
                       '_sorting_key_{}_'.format(self.term_order().name()))()

    def _monomial_one_(self):
        r"""
        Return the monomial `1` associated to this infinite polynomial ring.

        OUTPUT:

        :class:`Monomial`

        TESTS::

            sage: P.<x> = InfinitePolynomialRing(QQ, order='degrevlex', implementation='sparse_exponents')
            sage: repr(P._monomial_one_()) == ''
            True
        """
        return Monomial(tuple({} for _ in self._names_))

    def _element_constructor_(self, data):
        r"""
        Convert ``data`` to a polynomial.

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')

            sage: P(x[0])
            x_0

            sage: P({x[13]: 3, y[14]: 4})
            3*x_13 + 4*y_14
            sage: mx = next(iter(x[0]._summands_))
            sage: P({mx: 2})
            2*x_0

            sage: P(int(0))
            0
            sage: P(QQ(0))
            0
            sage: P(3/2)
            3/2

            sage: Q = InfinitePolynomialRing(QQ, names=('y', 'x'), order='deglex', implementation='sparse_exponents')
            sage: Q(x[0]), Q(y[0])
            (x_0, y_0)
            sage: Q(x[0] * y[0])
            y_0*x_0
            sage: Q(x[0] + y[0])
            y_0 + x_0
            sage: Q(x[42]^3 * y[24]^5)
            y_24^5*x_42^3

            sage: R = InfinitePolynomialRing(QQ, names=('x',), order='deglex', implementation='sparse_exponents')
            sage: R(x[1])
            x_1
            sage: R(y[3])
            Traceback (most recent call last):
            ...
            ValueError: cannot convert y_3
            to Infinite polynomial ring in x over Rational Field
            > *previous* ValueError: 'y' does not specify a generator
            of Infinite polynomial ring in x over Rational Field
            sage: R(x[3]^2 + x[2]^5 * x[1])
            x_1*x_2^5 + x_3^2

            sage: Z.<z> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: R(z[3])
            x_3

            sage: S = InfinitePolynomialRing(QQ, names=('x', 'y'), order='deglex', implementation='sparse')
            sage: sx, sy = S.gens()
            sage: P(sx[4] + 2*sy[3]^5)
            2*y_3^5 + x_4

            sage: D = InfinitePolynomialRing(QQ, names=('x', 'y'), order='deglex', implementation='dense')
            sage: dx, dy = D.gens()
            sage: P(dx[4] + 2*dy[3]^5)
            2*y_3^5 + x_4

            sage: F = PolynomialRing(QQ, names=('x_4', 'y_3'), order='deglex')
            sage: fx, fy = F.gens()
            sage: P(fx + 2*fy^5)
            2*y_3^5 + x_4
        """
        from sage.rings.asymptotic.misc import combine_exceptions

        def rewire_summands(data, summands, names):
            if self.ngens() == 1 and len(names) == 1:
                return summands

            def map_index(index):
                try:
                    return self._index_by_name_(names[index])
                except ValueError as e:
                    raise combine_exceptions(
                        ValueError('cannot convert {} to {}'.format(
                            data, self)), e)

            def rewire(monomial):
                rewired_monomial = [{} for _ in self._names_]
                for index, component in enumerate(monomial._exponents_):
                    if component:
                        rewired_monomial[map_index(index)] = component
                return Monomial(rewired_monomial)

            return {rewire(monomial): coefficient
                    for monomial, coefficient
                    in iteritems(summands)}

        if isinstance(data, dict):
            return self.element_class(self, data)

        elif isinstance(data, InfinitePolynomial_sparse_exponents):
            summands = rewire_summands(data,
                                       data._summands_,
                                       data.parent()._names_)
            return self.element_class(self, summands)

        elif data == 0:
            return self.element_class(self, {})

        try:
            return self.element_class(self, {self._monomial_one_(): data})
        except (TypeError, ValueError):
            pass

        F, B = self.construction(_implementation='sparse')
        P = F(B)
        try:
            data = P(data)
        except (TypeError, ValueError) as e:
            raise combine_exceptions(
                ValueError('cannot convert {} to {}'.format(
                    data, self)), e)

        summands = rewire_summands(data,
                                   {monomial_factory(monomial): coefficient
                                    for coefficient, monomial in data},
                                   data.parent()._names)
        return self.element_class(self, summands)

    def _coerce_map_from_(self, R):
        r"""
        Does ``R`` coerce to this infinite polynomial ring?

        INPUT:

        - ``R`` -- parent

        TESTS::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P.has_coerce_map_from(P)
            True
            sage: P.has_coerce_map_from(QQ)
            True
            sage: P.has_coerce_map_from(ZZ)
            True

            sage: Q.<y, z, x> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P.has_coerce_map_from(Q)
            False
            sage: Q.has_coerce_map_from(P)
            True

            sage: S = InfinitePolynomialRing(QQ, names=('x', 'y'), order='deglex', implementation='sparse')
            sage: P.has_coerce_map_from(S)
            True
            sage: D = InfinitePolynomialRing(QQ, names=('x', 'y'), order='deglex', implementation='dense')
            sage: P.has_coerce_map_from(D)
            True
            sage: F = PolynomialRing(QQ, names=('x_4', 'y_3'), order='deglex')
            sage: P.has_coerce_map_from(F)  # not tested; see :trac:`23632`.
            True
        """
        if isinstance(R, InfinitePolynomialRing_sparse_exponents):
            return all(name in self._names_ for name in R._names_)

        F, B = self.construction(_implementation='sparse')
        if F(B).has_coerce_map_from(R):
            return True

    def construction(self, _implementation='sparse_exponents'):
        r"""
        Return the construction of this infinite polynomial ring.

        OUTPUT:

        A pair ``(F, R)`` with
        - a construction functor ``F`` and
        - a ring ``R``
        such that ``F(R)`` is this infinite polynomial ring.

        EXAMPLES::

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex', implementation='sparse_exponents')
            sage: P.construction()
            (InfPoly{[x,y], "deglex", "sparse_exponents"},
             Rational Field)
        """
        from sage.categories.pushout import InfinitePolynomialFunctor
        return (InfinitePolynomialFunctor(self._names_,
                                          self._order_.name(),
                                          _implementation),
                self.coefficient_ring())
