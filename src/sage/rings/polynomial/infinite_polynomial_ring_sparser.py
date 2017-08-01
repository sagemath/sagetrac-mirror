r"""
Infinite Multivariate Polynomial Ring with Sparse Exponents

::

    sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing

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
    sage: P.<A> = InfinitePolynomialRing(QQ, order='deglex')
    sage: A[a(1,1)]
    A_1_1
    sage: A[a(1,1)]^2 * A[a(2,1)] * A[a(1,2)]^3 + 3*A[a(4,4)]^7
    3*A_4_4^7 + A_1_1^2*A_2_1*A_1_2^3

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


def updated_by_adding_coefficients(D, E):
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
    TESTS::

        sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
        sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import monomial_factory
        sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')

        sage: mx = next(iter(x[0]._summands_))
        sage: monomial_factory(mx)
        x0_0
        sage: monomial_factory(x[13])
        x0_13
        sage: monomial_factory(({3: 4}, {5: 6, 7: 8}))
        x0_3^4*x1_5^6*x1_7^8
    """
    if isinstance(data, Monomial):
        return data
    elif isinstance(data, InfinitePolynomial_sparser):
        summands = data._summands_
        if len(summands) != 1:
            raise ValueError('{} is not monomial'.format(data))
        monomial, coefficient = next(iteritems(summands))
        if coefficient != 1:
            raise ValueError('{} is not normalized monomial')
        return monomial
    else:
        return Monomial(data)


class Monomial(object):

    def __init__(self, exponents):
        self._exponents_ = tuple({index: int(exponent)
                                  for index, exponent in iteritems(component)
                                  if exponent != 0}
                                 for component in exponents)

    def __repr__(self, names=None):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import Monomial
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
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import Monomial
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
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import Monomial

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
            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import Monomial

            sage: x_0 = Monomial(({0: 1},))
            sage: x_0._sorting_key_lex_()
            ((0, 1),)
            sage: x_0._sorting_key_revlex_()
            ((0, -1),)
        """
        return tuple(tuple(-t for t in T) for T in self._sorting_key_lex_())

    def _sorting_key_deglex_(self):
        return (self.degree(), self._sorting_key_lex_())

    def _sorting_key_degrevlex_(self):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import Monomial
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='degrevlex')
            sage: x[0] + y[1]
            x_0 + y_1
        """
        return (self.degree(), self._sorting_key_revlex_())

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import Monomial
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
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
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
            sage: mx = next(iter(x[0]._summands_))
            sage: mx == mx
            True
            sage: my = next(iter(y[0]._summands_))
            sage: mx == my
            False
        """
        return self._sorting_key_lex_() == other._sorting_key_lex_()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __mul__(self, other):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import Monomial
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
            tuple(updated_by_adding_coefficients(scomponent, ocomponent)
                  for scomponent, ocomponent
                  in zip(self._exponents_, other._exponents_)))


class InfinitePolynomial_sparser(CommutativeAlgebraElement):
    def __init__(self, parent, data):
        super(InfinitePolynomial_sparser, self).__init__(parent=parent)

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
        parent = self.parent()
        return iter((coefficient, InfinitePolynomial_sparser(parent, monomial))
                    for monomial, coefficient in iteritems(self._summands_))

    @cached_method
    def _sorted_monomials_and_coefficients_(self, reverse=False):
        def key(monomial_and_coefficient):
            monomial = monomial_and_coefficient[0]
            return self.parent()._sorting_key_monomial_(monomial)
        return sorted(iteritems(self._summands_), key=key, reverse=reverse)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
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

            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='degrevlex')
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
        return hash((self.parent(),
                     tuple(self._sorted_monomials_and_coefficients_(reverse=True))))

    def __bool__(self):
        return bool(self._summands_)

    __nonzero__ = __bool__

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
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
        if (isinstance(other, InfinitePolynomial_sparser)
            and not other._summands_):
                return not bool(self)
        try:
            return not bool(self - other)
        except (TypeError, ValueError):
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def _add_(self, other):
        r"""
        EXAMPLES::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
            sage: x[0] + y[0]
            x_0 + y_0
            sage: x[0] + y[0] + x[0]
            2*x_0 + y_0

        TESTS::

            sage: x[0] + 1
            x_0 + 1
        """
        summands = updated_by_adding_coefficients(
            self._summands_, other._summands_)
        return self.parent().element_class(self.parent(), summands)

    def _lmul_(self, other):
        r"""
        EXAMPLES::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
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
        EXAMPLES::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
            sage: x[0] - x[0]
            0
            sage: x[1] - y[2]
            x_1 - y_2
            sage: y[1] - x[2]
            -x_2 + y_1
        """
        return self + self.parent().coefficient_ring()(-1) * other

    def _mul_by_summand_(self, other_monomial, other_coefficient):
        summands = {monomial*other_monomial: coefficient*other_coefficient
                    for monomial, coefficient in iteritems(self._summands_)}
        return self.parent().element_class(self.parent(), summands)

    def _mul_(self, other):
        return sum((self._mul_by_summand_(monomial, coefficient)
                    for monomial, coefficient in iteritems(other._summands_)),
                   self.parent().zero())

    def degree(self):
        r"""
        EXAMPLES::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
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


class InfinitePolynomialGen_sparser(InfinitePolynomialGen_generic):

    def __init__(self, parent, name, index):
        self._index = index
        super(InfinitePolynomialGen_sparser, self).__init__(parent, name)

    def __hash__(self):
        return hash((self._parent, self._name, self._index))

    @cached_method
    def __getitem__(self, i):
        monomial = Monomial(
            tuple(({i: 1} if index == self._index else {})
                  for index in range(len(self._parent._names_))))
        return InfinitePolynomial_sparser(self._parent, monomial)


class InfinitePolynomialRing_sparser(Algebra, UniqueRepresentation):

    Element = InfinitePolynomial_sparser

    @staticmethod
    def __classcall__(cls, coefficient_ring, names,
                      order='lex',
                      category=None):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
            sage: P
            Infinite polynomial ring in x, y over Rational Field

        ::

            sage: InfinitePolynomialRing(QQ, names=('x', 'y'), order='deglex') is P
            True
            sage: InfinitePolynomialRing(QQ, ('x', 'y'), order='deglex') is P
            True
        """

        names = tuple(names)

        if category is None:
            from sage.categories.commutative_algebras import CommutativeAlgebras
            from sage.categories.rings import Rings
            category = CommutativeAlgebras(Rings())

        return super(InfinitePolynomialRing_sparser,
                     cls).__classcall__(cls, coefficient_ring, names,
                                        order=order,
                                        category=category)

    def __init__(self, coefficient_ring, names, order, category):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
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

        super(InfinitePolynomialRing_sparser, self).__init__(
            base_ring=coefficient_ring,
            category=category)

    def _repr_(self):
        return 'Infinite polynomial ring in {} over {}'.format(
            ', '.join(self._names_),
            self.coefficient_ring())

    @cached_method
    def __hash__(self):
        return hash((self.coefficient_ring(),
                     self._names_,
                     self.term_order()))

    def coefficient_ring(self):
        return self._coefficient_ring_

    def term_order(self):
        return self._order_

    @cached_method
    def gens(self):
        r"""
        EXAMPLES::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')  # indirect doctest
            sage: x
            x_*
            sage: x[3]
            x_3
            sage: y
            y_*
            sage: y[42]
            y_42
        """
        return tuple(InfinitePolynomialGen_sparser(self, name, index)
                     for index, name in enumerate(self._names_))

    def gen(self, n=0):
        return self.gens()[0]

    def _index_by_name_(self, name):
        try:
            return self._names_.index(name)
        except ValueError:
            raise ValueError("'{}' does not specify a generator of {}".format(
                name, self))

    def gen_by_name(self, name):
        return self.gen(self._index_by_name_(name))

    def ngens(self):
        return len(self.gens())

    def _sorting_key_monomial_(self, monomial):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x> = InfinitePolynomialRing(QQ, order='deglex')
            sage: monomial = next(iter(x[0]._summands_))
            sage: P._sorting_key_monomial_(monomial)
            (1, ((0, 1),))

        ::

            sage: P.<x> = InfinitePolynomialRing(QQ, order='degrevlex')
            sage: monomial = next(iter(x[0]._summands_))
            sage: P._sorting_key_monomial_(monomial)
            (1, ((0, -1),))
        """
        return getattr(monomial,
                       '_sorting_key_{}_'.format(self.term_order().name()))()

    def _monomial_one_(self):
        return Monomial(tuple({} for _ in self._names_))

    def _element_constructor_(self, data):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')

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

            sage: Q = InfinitePolynomialRing(QQ, names=('y', 'x'), order='deglex')
            sage: Q(x[0]), Q(y[0])
            (x_0, y_0)
            sage: Q(x[0] * y[0])
            y_0*x_0
            sage: Q(x[0] + y[0])
            y_0 + x_0
            sage: Q(x[42]^3 * y[24]^5)
            y_24^5*x_42^3

            sage: R = InfinitePolynomialRing(QQ, names=('x',), order='deglex')
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

            sage: Z.<z> = InfinitePolynomialRing(QQ, order='deglex')
            sage: R(z[3])
            x_3
        """
        from sage.rings.asymptotic.misc import combine_exceptions

        if isinstance(data, dict):
            return self.element_class(self, data)

        elif isinstance(data, InfinitePolynomial_sparser):
            if self.ngens() == 1 and data.parent().ngens() == 1:
                return self.element_class(self, data._summands_)
            else:
                def map_index(index):
                    try:
                        return self._index_by_name_(data.parent()._names_[index])
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

                summands = {rewire(monomial): coefficient
                            for monomial, coefficient
                            in iteritems(data._summands_)}
                return self.element_class(self, summands)

        elif data == 0:
            return self.element_class(self, {})
        else:
            return self.element_class(self, {self._monomial_one_(): data})


def InfinitePolynomialRing(*args, **kwds):
    return InfinitePolynomialRing_sparser(*args, **kwds)
