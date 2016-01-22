"""
Support Python's numbers abstract base class

.. SEEALSO:: :pep:`3141` for more information about :class:`numbers`.
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import numbers


def register_sage_classes():
    """
    Register all relevant Sage classes in the :class:`numbers`
    hierarchy.

    EXAMPLES::

        sage: import numbers
        sage: isinstance(5, numbers.Integral)
        True
        sage: isinstance(5, numbers.Number)
        True
        sage: isinstance(5/1, numbers.Integral)
        False
        sage: isinstance(22/7, numbers.Rational)
        True
        sage: isinstance(1.3, numbers.Real)
        True
        sage: isinstance(CC(1.3), numbers.Real)
        False
        sage: isinstance(CC(1.3 + I), numbers.Complex)
        True
        sage: isinstance(RDF(1.3), numbers.Real)
        True
        sage: isinstance(CDF(1.3, 4), numbers.Complex)
        True
        sage: isinstance(AA(sqrt(2)), numbers.Real)
        True
        sage: isinstance(QQbar(I), numbers.Complex)
        True

    This doesn't work with symbolic expressions at all::

        sage: isinstance(pi, numbers.Real)
        False
        sage: isinstance(I, numbers.Complex)
        False
        sage: isinstance(sqrt(2), numbers.Real)
        False
    """
    from sage.rings.integer import Integer
    from sage.rings.rational import Rational
    from sage.rings.real_mpfr import RealNumber
    from sage.rings.real_double import RealDoubleElement
    from sage.rings.complex_number import ComplexNumber
    from sage.rings.complex_double import ComplexDoubleElement
    from sage.rings.complex_mpc import MPComplexNumber
    from sage.rings.qqbar import AlgebraicReal, AlgebraicNumber

    numbers.Integral.register(Integer)
    numbers.Rational.register(Rational)
    numbers.Real.register(RealNumber)
    numbers.Real.register(RealDoubleElement)
    numbers.Real.register(AlgebraicReal)
    numbers.Complex.register(ComplexNumber)
    numbers.Complex.register(MPComplexNumber)
    numbers.Complex.register(ComplexDoubleElement)
    numbers.Complex.register(AlgebraicNumber)

register_sage_classes()
