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

    def _sorting_key_deglex_(self):
        return (self.deg(), self._sorting_key_revlex_())


