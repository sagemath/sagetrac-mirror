r"""
The field of real number.

.. TODO::

    The real lazy field (:class:`sage.rings.real_lazy.RealLazyField_class`) should be
    considered as the set of *computable* numbers and not the set of real numbers. The
    precise interaction between the real lazy field and the ``RealFloatingPointField`` defined in
    this module has to be precised.
"""
#*****************************************************************************
#       Copyright (C) 2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import

from sage.misc.fast_methods import Singleton

from sage.structure.element import parent
from sage.structure.parent import Parent

from .infinity import Infinity

def create_RealFloatingPointField(prec=53, type=None, rnd="RNDN", sci_not=0):
    """
    Create a real field with given precision, type, rounding mode and
    scientific notation.

    Some options are ignored for certain types (RDF for example).

    INPUT:

    - ``prec`` -- a positive integer or infinity

    - ``type`` -- type of real field:

      - ``'RR'`` -- the field of real numbers
      - ``'RDF'`` -- the Sage real field corresponding to native doubles
      - ``'MPFR'`` or ``'FloatingPoint'`` -- floating point real numbers
        implemented using the MPFR
      - ``'MPIR'`` or ``'Interval'`` -- real fields implementing interval arithmetic
      - ``'ARB'`` or ``'Ball'`` -- real fields implementing ball arithmetic
      - ``'RLF'`` -- the real lazy field

    - ``rnd`` -- rounding mode:

      - ``'RNDN'`` -- round to nearest
      - ``'RNDZ'`` -- round toward zero
      - ``'RNDD'`` -- round down
      - ``'RNDU'`` -- round up

    - ``sci_not`` -- boolean, whether to use scientific notation for printing

    EXAMPLES::

        sage: from sage.rings.real_field import create_RealFloatingPointField
        sage: create_RealFloatingPointField(30)
        Real Field with 30 bits of precision
        sage: create_RealFloatingPointField(oo)
        Real Field
        sage: create_RealFloatingPointField(20, 'RDF') # ignores precision
        Real Double Field
        sage: create_RealFloatingPointField(60, 'Interval')
        Real Interval Field with 60 bits of precision
        sage: create_RealFloatingPointField(40, 'RLF') # ignores precision
        Real Lazy Field
    """
    if type is None:
        if prec == Infinity:
            type = "RR"
        else:
            type = "MPFR"

    if type == "RDF" or type == "Double":
        from .real_double import RDF
        return RDF
    elif type == "MPIR" or type == "Interval":
        from .real_mpfi import RealIntervalField
        return RealIntervalField(prec, sci_not)
    elif type == "ARB" or type == "Ball":
        from .real_arb import RealBallField
        return RealBallField(prec)
    elif type == "MPFR" or type == "FloatingPoint":
        from .real_mpfr import RealFloatingPointField
        return RealFloatingPointField(prec, sci_not, rnd)
    elif type == "RLF":
        from .real_lazy import RLF
        return RLF
    elif type == "RR":
        from .real_field import RealFloatingPointField
        return RealFloatingPointField()
    else:
        raise ValueError('unrecognized type {!r} for real field'.format(type))


class RealFloatingPointField(Singleton, Parent):
    r"""
    The field of real numbers.

    EXAMPLES::

        sage: R = QQ.completion(oo, oo)

        sage: 3 in R
        True
        sage: 1/2 in R
        True
        sage: AA(5).sqrt() in R
        True

    The behavior with respect to number fields is certainly not desirable::

        sage: K.<sqrt2> = QuadraticField(2, embedding=AA(2).sqrt())
        sage: sqrt2 in R  # good
        True
        sage: R(sqrt2)    # bad
        1.414213562373095?
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.rings.real_field import RealFloatingPointField
            sage: R = RealFloatingPointField()
            sage: TestSuite(R).run()
            sage: loads(dumps(R)) is R
            True
        """
        from sage.categories.fields import Fields
        from sage.rings.rational_field import QQ
        from sage.rings.qqbar import AA
        Parent.__init__(self,
                facade=[QQ, AA],
                category=Fields().Metric().Complete().Infinite())

    def metric(self):
        r"""
        EXAMPLES::

            sage: R = QQ.completion(oo, oo)
            sage: phi = R.metric()
            sage: phi(1, 3/4)
            1/4
            sage: phi(AA(2).sqrt(), 1)
            0.4142135623730951?
        """
        return lambda a,b: (a-b).abs()

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.rings.real_field import RealFloatingPointField
            sage: RealFloatingPointField()
            Real Field
        """
        return "Real Field"

    def _latex_(self):
        r"""
        TESTS::

            sage: R = QQ.completion(oo, oo)
            sage: latex(R)
            \Bold{R}
        """
        return "\\Bold{R}"

    def is_exact(self):
        r"""
        TESTS::

            sage: QQ.completion(oo, oo).is_exact()
            True
        """
        return True

    def construction(self):
        r"""
        Return the functorial construction of the real field.

        EXAMPLES::

            sage: R = QQ.completion(oo, oo)
            sage: functor, base = R.construction()
            sage: functor
            Completion[+Infinity, prec=+Infinity]
            sage: base
            Rational Field
            sage: functor(base) is R
            True
        """
        from sage.categories.pushout import CompletionFunctor
        from sage.rings.rational_field import QQ
        return (CompletionFunctor(Infinity, Infinity, {}), QQ)

    def _coerce_map_from_(self, S):
        r"""
        TESTS::

            sage: R = QQ.completion(oo, oo)
            sage: R.has_coerce_map_from(QQ)
            True
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.qqbar import AA
        if S is ZZ or S is QQ or S is AA:
            return True

        from sage.rings.number_field.number_field_base import is_NumberField
        if is_NumberField(S) and S.coerce_embedding() is not None:
            return AA.has_coerce_map_from(S.coerce_embedding().codomain())

    def __contains__(self, x):
        r"""
        TESTS::

            sage: R = QQ.completion(oo, oo)
            sage: 1 in R
            True
        """
        return self.has_coerce_map_from(parent(x))

    def some_elements(self):
        r"""
        Some elements in the field of real numbers.

        EXAMPLES::

            sage: R = QQ.completion(oo, oo)
            sage: R.some_elements()
            (0, 1, -1, 2, -2, 1/2, -1/2, 2, -2, 0, 1, 1, a, 2*a, -1/3*a + 2/3, 1/2*a)
        """
        from itertools import islice
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.qqbar import AA
        from sage.rings.number_field.number_field import NumberField
        from sage.rings.polynomial.polynomial_ring import polygen

        x = polygen(QQ)
        nf = NumberField(x**3 - 2, 'a', embedding=AA(2)**QQ((1,3)))

        subsets = [ZZ, QQ, AA, nf]
        l = []
        for k in subsets:
            l.extend(islice(k.some_elements(), 5))
        return tuple(l)
