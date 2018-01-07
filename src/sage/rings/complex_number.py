r"""
Deprecated module

TESTS::

    sage: from sage.rings.complex_number import create_ComplexNumber
    sage: create_ComplexNumber(0, 1)
    doctest:warning
    ...
    DeprecationWarning:
    Importing create_ComplexNumber from here is deprecated.
    If you need to use it, please import it directly from sage.rings.complex_mpfr
    See http://trac.sagemath.org/24483 for details.
    1.00000000000000*I
"""
#*****************************************************************************
#       Copyright (C) 2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function

from sage.misc.lazy_import import LazyImport
ComplexNumber = LazyImport('sage.rings.complex_mpfr', 'ComplexNumber', deprecation=24483)
create_ComplexNumber = LazyImport('sage.rings.complex_mpfr', 'create_ComplexNumber', deprecation=24483)
is_ComplexNumber = LazyImport('sage.rings.complex_mpfr', 'is_ComplexNumber', deprecation=24483)
set_global_complex_round_mode = LazyImport('sage.rings.complex_mpfr', 'set_global_complex_round_mode', deprecation=24483)
CCtoCDF = LazyImport('sage.rings.complex_mpfr', 'CCtoCDF', deprecation=24483)
RRtoCC = LazyImport('sage.rings.complex_mpfr', 'RRtoCC', deprecation=24483)
