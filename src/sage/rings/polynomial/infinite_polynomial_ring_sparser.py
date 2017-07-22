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
