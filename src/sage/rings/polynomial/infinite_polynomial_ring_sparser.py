r"""
Infinite Multivariate Polynomial Ring with Sparse Exponents

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

from .infinite_polynomial_ring import InfinitePolynomialGen as InfinitePolynomialGen_generic


def monomial_factory(data):
    if isinstance(data, Monomial):
        return data
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
                        for index, exponent in iteritems(component))

    __str__ = __repr__

    def deg(self):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import Monomial
            sage: a = Monomial(({0: 2, 2: 3}, {0: 3, 3: 1})); a
            x0_0^2*x0_2^3*x1_0^3*x1_3
            sage: a.deg()
            9
        """
        return sum(exponent
                   for component in self._exponents_
                   for exponent in itervalues(component))

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

            sage: names = ('x', 'y')
            sage: x_0 = Monomial(({0: 1}, {}))
            sage: x_1 = Monomial(({1: 1}, {}))
            sage: y_0 = Monomial(({}, {0: 1}))
            sage: y_1 = Monomial(({}, {1: 1}))
            sage: x_0._sorting_key_lex_() < x_1._sorting_key_lex_()
            True
            sage: y_0._sorting_key_lex_() < y_1._sorting_key_lex_()
            True
            sage: x_0._sorting_key_lex_() < y_0._sorting_key_lex_()
            True
            sage: x_1._sorting_key_lex_() < y_0._sorting_key_lex_()
            True
            sage: x_0._sorting_key_lex_() < y_1._sorting_key_lex_()
            True
            sage: x_1._sorting_key_lex_() < y_1._sorting_key_lex_()
            True
            

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
            sage: a._sorting_key_lex_() < c._sorting_key_lex_()
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
                     for component in reversed(self._exponents_))

    def _sorting_key_revlex_(self):
        return tuple(tuple(-t for t in T) for T in self._sorting_key_lex_)

    def _sorting_key_deglex_(self):
        return (self.deg(), self._sorting_key_lex_())

    def _sorting_key_degrevlex_(self):
        return (self.deg(), self._sorting_key_revlex_())

    def __hash__(self):
        return hash(self._sorting_key_lex_())

class InfinitePolynomial_sparser(CommutativeAlgebraElement):
    def __init__(self, parent, data, _copy=True):
        super(InfinitePolynomial_sparser, self).__init__(parent=parent)

        coefficient_ring = parent.coefficient_ring()
        if isinstance(data, Monomial):
            self._summands_ = {data: coefficient_ring(1)}
        elif _copy:
            self._summands_ = {monomial_factory(monomial):
                               coefficient_ring(coefficient)
                               for monomial, coefficient in iteritems(data)}
        else:
            self._summands_ = data

    def __iter__(self):
        parent = self.parent()
        return iter((coefficient, InfinitePolynomial_sparser(parent, monomial))
                    for monomial, coefficient in iteritems(self._summands_))

    @cached_method
    def _sorted_monomials_and_coefficients_(self):
        def key(monomial_and_coefficient):
            monomial = monomial_and_coefficient[0]
            return self.parent()._sorting_key_monomial_(monomial)
        return sorted(iteritems(self._summands_), key=key, reverse=True)

    def _repr_(self):
        def summand(monomial, coefficient):
            factors = ('{}*'.format(coefficient) if coefficient != 1 else '',
                       monomial.__repr__(self.parent()._names_))
            s = '*'.join(f for f in factors if f)
            return s or '1'

        r = ' + '.join(summand(monomial, coefficient)
                       for monomial, coefficient
                       in self._sorted_monomials_and_coefficients_())
        return r or '0'

    def __hash__(self):
        return hash((self.parent(),
                     tuple(self._sorted_monomials_and_coefficients_())))


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

    def __init__(self, coefficient_ring, names, order='lex'):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_ring_sparser import InfinitePolynomialRing
            sage: P.<x, y> = InfinitePolynomialRing(QQ, order='deglex')
            sage: P
            Infinite polynomial ring in x, y over Rational Field
        """
        from sage.rings.polynomial.term_order import TermOrder
        from sage.symbolic.ring import isidentifier

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
            category=None)

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
        TESTS::

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
        """
        return getattr(monomial,
                       '_sorting_key_{}_'.format(self.term_order().name()))()

    def _element_constructor_(self, data):
        return self.element_class(self, data)


def InfinitePolynomialRing(*args, **kwds):
    return InfinitePolynomialRing_sparser(*args, **kwds)
