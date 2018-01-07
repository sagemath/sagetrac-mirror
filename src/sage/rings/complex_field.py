r"""
Deprecated module

TESTS::

    sage: from sage.rings.complex_field import ComplexField
    sage: ComplexField()   # not deprecated
    doctest:warning
    ...
    DeprecationWarning:
    Importing ComplexField from here is deprecated. If you need to use it,
    please import it directly from sage.rings.complex_mpfr
    See http://trac.sagemath.org/24483 for details.
    Complex Field with 53 bits of precision
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

is_ComplexField = LazyImport('sage.rings.complex_mpfr', 'is_ComplexField', deprecation=24483)
ComplexField_class = LazyImport('sage.rings.complex_mpfr', 'ComplexField_class', deprecation=24483)
ComplexField = LazyImport('sage.rings.complex_mpfr', 'ComplexField', deprecation=24483)
