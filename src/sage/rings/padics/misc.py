r"""
Miscellaneous Functions

This file contains two miscellaneous functions used by `p`-adics.

- ``min`` -- a version of ``min`` that returns `\infty` on empty input.
- ``max`` -- a version of ``max`` that returns `-\infty` on empty input.

AUTHORS:

- David Roe
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from six.moves.builtins import min as python_min
from six.moves.builtins import max as python_max
from sage.rings.infinity import infinity


def min(*L):
    r"""
    Return the minimum of the inputs, where the minimum of the empty
    list is `\infty`.

    EXAMPLES::

        sage: from sage.rings.padics.misc import min
        sage: min()
        +Infinity
        sage: min(2,3)
        2
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_min(L)
    except ValueError:
        return infinity


def max(*L):
    r"""
    Return the maximum of the inputs, where the maximum of the empty
    list is `-\infty`.

    EXAMPLES::

        sage: from sage.rings.padics.misc import max
        sage: max()
        -Infinity
        sage: max(2,3)
        3
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_max(L)
    except ValueError:
        return -infinity

def precprint(prec_type, prec_cap, p):
    """
    String describing the precision mode on a p-adic ring or field.

    EXAMPLES::

        sage: from sage.rings.padics.misc import precprint
        sage: precprint('capped-rel', 12, 2)
        'with capped relative precision 12'
        sage: precprint('capped-abs', 11, 3)
        'with capped absolute precision 11'
        sage: precprint('floating-point', 1234, 5)
        'with floating precision 1234'
        sage: precprint('fixed-mod', 1, 17)
        'of fixed modulus 17^1'
    """
    precD = {'capped-rel':'with capped relative precision %s'%prec_cap,
             'capped-abs':'with capped absolute precision %s'%prec_cap,
             'floating-point':'with floating precision %s'%prec_cap,
             'fixed-mod':'of fixed modulus %s^%s'%(p, prec_cap)}
    return precD[prec_type]
