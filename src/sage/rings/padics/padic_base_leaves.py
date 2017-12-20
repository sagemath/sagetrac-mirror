"""
`p`-Adic Base Leaves

Implementations of `\mathbb{Z}_p` and `\mathbb{Q}_p`

AUTHORS:

- David Roe
- Genya Zaytman: documentation
- David Harvey: doctests
- William Stein: doctest updates

EXAMPLES:

`p`-Adic rings and fields are examples of inexact structures, as the
reals are.  That means that elements cannot generally be stored
exactly: to do so would take an infinite amount of storage.  Instead,
we store an approximation to the elements with varying precision.

There are two types of precision for a `p`-adic element.  The first is
relative precision, which gives the number of known `p`-adic digits::

    sage: R = Qp(5, 20, 'capped-rel', 'series'); a = R(675); a
    2*5^2 + 5^4 + O(5^22)
    sage: a.precision_relative()
    20

The second type of precision is absolute precision, which gives the
power of `p` that this element is stored modulo::

    sage: a.precision_absolute()
    22

The number of times that `p` divides the element is called the
valuation, and can be accessed with the functions ``valuation()`` and
``ordp()``:

    sage: a.valuation()
    2

The following relationship holds:

``self.valuation() + self.precision_relative() == self.precision_absolute().``

    sage: a.valuation() + a.precision_relative() == a.precision_absolute()
    True

In the capped relative case, the relative precision of an element
is restricted to be at most a certain value, specified at the
creation of the field.  Individual elements also store their own
precision, so the effect of various arithmetic operations on
precision is tracked.  When you cast an exact element into a
capped relative field, it truncates it to the precision cap of the
field.::

    sage: R = Qp(5, 5); a = R(4006); a
    1 + 5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(17/3); b
    4 + 2*5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5)
    sage: c = R(4025); c
    5^2 + 2*5^3 + 5^4 + 5^5 + O(5^7)
    sage: a + b
    4*5 + 3*5^2 + 3*5^3 + 4*5^4 + O(5^5)
    sage: a + b + c
    4*5 + 4*5^2 + 5^4 + O(5^5)

::

    sage: R = Zp(5, 5, 'capped-rel', 'series'); a = R(4006); a
    1 + 5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(17/3); b
    4 + 2*5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5)
    sage: c = R(4025); c
    5^2 + 2*5^3 + 5^4 + 5^5 + O(5^7)
    sage: a + b
    4*5 + 3*5^2 + 3*5^3 + 4*5^4 + O(5^5)
    sage: a + b + c
    4*5 + 4*5^2 + 5^4 + O(5^5)

In the capped absolute type, instead of having a cap on the
relative precision of an element there is instead a cap on the
absolute precision.  Elements still store their own precisions,
and as with the capped relative case, exact elements are truncated
when cast into the ring.::

    sage: R = ZpCA(5, 5); a = R(4005); a
    5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(4025); b
    5^2 + 2*5^3 + 5^4 + O(5^5)
    sage: a * b
    5^3 + 2*5^4 + O(5^5)
    sage: (a * b) // 5^3
    1 + 2*5 + O(5^2)
    sage: type((a * b) // 5^3)
    <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
    sage: (a * b) / 5^3
    1 + 2*5 + O(5^2)
    sage: type((a * b) / 5^3)
    <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>

The fixed modulus type is the leanest of the p-adic rings: it is
basically just a wrapper around `\mathbb{Z} / p^n \mathbb{Z}`
providing a unified interface with the rest of the `p`-adics.  This is
the type you should use if your primary interest is in speed (though
it's not all that much faster than other `p`-adic types).  It does not
track precision of elements.::

    sage: R = ZpFM(5, 5); a = R(4005); a
    5 + 2*5^3 + 5^4 + O(5^5)
    sage: a // 5
    1 + 2*5^2 + 5^3 + O(5^5)

`p`-Adic rings and fields should be created using the creation
functions ``Zp`` and ``Qp`` as above.  This will ensure that there is
only one instance of `\mathbb{Z}_p` and `\mathbb{Q}_p` of a given
type, `p`, print mode and precision.  It also saves typing very long
class names.::

    sage: Qp(17,10)
    17-adic Field with capped relative precision 10
    sage: R = Qp(7, prec = 20, print_mode = 'val-unit'); S = Qp(7, prec = 20, print_mode = 'val-unit'); R is S
    True
    sage: Qp(2)
    2-adic Field with capped relative precision 20

Once one has a `p`-Adic ring or field, one can cast elements into it
in the standard way.  Integers, ints, longs, Rationals, other `p`-Adic
types, pari `p`-adics and elements of `\mathbb{Z} / p^n \mathbb{Z}`
can all be cast into a `p`-Adic field.::

    sage: R = Qp(5, 5, 'capped-rel','series'); a = R(16); a
    1 + 3*5 + O(5^5)
    sage: b = R(23/15); b
    5^-1 + 3 + 3*5 + 5^2 + 3*5^3 + O(5^4)
    sage: S = Zp(5, 5, 'fixed-mod','val-unit'); c = S(Mod(75,125)); c
    5^2 * 3 + O(5^5)
    sage: R(c)
    3*5^2 + O(5^5)

In the previous example, since fixed-mod elements don't keep track
of their precision, we assume that it has the full precision of
the ring.  This is why you have to cast manually here.

While you can cast explicitly as above, the chains of automatic
coercion are more restricted.  As always in Sage, the following
arrows are transitive and the diagram is commutative.::

    int -> long -> Integer -> Zp capped-rel -> Zp capped_abs -> IntegerMod
    Integer -> Zp fixed-mod -> IntegerMod
    Integer -> Zp capped-abs -> Qp capped-rel

In addition, there are arrows within each type.  For capped relative
and capped absolute rings and fields, these arrows go from lower
precision cap to higher precision cap.  This works since elements
track their own precision: choosing the parent with higher precision
cap means that precision is less likely to be truncated unnecessarily.
For fixed modulus parents, the arrow goes from higher precision cap to
lower.  The fact that elements do not track precision necessitates
this choice in order to not produce incorrect results.

TESTS::

    sage: R = Qp(5, 15, print_mode='bars', print_sep='&')
    sage: repr(R(2777))[3:]
    '0&0&0&0&0&0&0&0&0&0&4&2&1&0&2'
    sage: TestSuite(R).run()

    sage: R = Zp(5, 15, print_mode='bars', print_sep='&')
    sage: repr(R(2777))[3:]
    '0&0&0&0&0&0&0&0&0&0&4&2&1&0&2'
    sage: TestSuite(R).run()

    sage: R = ZpCA(5, 15, print_mode='bars', print_sep='&')
    sage: repr(R(2777))[3:]
    '0&0&0&0&0&0&0&0&0&0&4&2&1&0&2'
    sage: TestSuite(R).run()

"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.richcmp import op_LE

from .generic_nodes import pAdicFieldBaseGeneric, \
                          pAdicCappedRelativeFieldGeneric, \
                          pAdicRingBaseGeneric, \
                          pAdicCappedRelativeRingGeneric, \
                          pAdicFixedModRingGeneric, \
                          pAdicCappedAbsoluteRingGeneric, \
                          pAdicFloatingPointRingGeneric, \
                          pAdicFloatingPointFieldGeneric, \
                          pAdicGeneric
from .padic_capped_relative_element import pAdicCappedRelativeElement
from .padic_capped_absolute_element import pAdicCappedAbsoluteElement
from .padic_fixed_mod_element import pAdicFixedModElement
from .padic_floating_point_element import pAdicFloatingPointElement
from .padic_lattice_element import pAdicLatticeElement
from .lattice_precision import PrecisionLattice

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity


class pAdicRingCappedRelative(pAdicRingBaseGeneric, pAdicCappedRelativeRingGeneric):
    r"""
    An implementation of the `p`-adic integers with capped relative
    precision.
    """
    def __init__(self, p, prec, print_mode, names):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options.
        - ``names`` -- how to print the prime.

        EXAMPLES::

            sage: R = ZpCR(next_prime(10^60)) #indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicRingCappedRelative_with_category'>

        TESTS::

            sage: R = ZpCR(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^10)], max_runs = 2^12, skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = ZpCR(3, 1)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^3)])

            sage: R = ZpCR(3, 2)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^6)], skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = ZpCR(next_prime(10^60))
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^3)], max_runs = 2^5, skip='_test_log') # long time
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)]) # long time
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedRelativeElement)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = Zp(17)
            sage: K(1) + 1 #indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
            sage: K.has_coerce_map_from(Zp(17,10))
            True
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
        """
        #if isistance(R, pAdicRingLazy) and R.prime() == self.prime():
        #    return True
        if isinstance(R, pAdicRingCappedRelative) and R.prime() == self.prime():
            if R.precision_cap() < self.precision_cap():
                return True
            elif (R.precision_cap() == self.precision_cap() and
                  self._printer.richcmp_modes(R._printer, op_LE)):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: Zp(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Ring with capped relative precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

class pAdicRingCappedAbsolute(pAdicRingBaseGeneric, pAdicCappedAbsoluteRingGeneric):
    r"""
    An implementation of the `p`-adic integers with capped absolute precision.
    """
    def __init__(self, p, prec, print_mode, names):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options.
        - ``names`` -- how to print the prime.

        EXAMPLES::

            sage: R = ZpCA(next_prime(10^60)) #indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicRingCappedAbsolute_with_category'>

        TESTS::

            sage: R = ZpCA(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^10)], max_runs = 2^12, skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = ZpCA(3, 1)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^3)])

            sage: R = ZpCA(3, 2)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^6)], skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = ZpCA(next_prime(10^60))
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^3)], max_runs = 2^5, skip='_test_log') # long time
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedAbsoluteElement)

    def _coerce_map_from_(self, R):
        """
        Returns ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = ZpCA(17)
            sage: K(1) + 1 #indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
            sage: K.has_coerce_map_from(ZpCA(17,10))
            True
            sage: K.has_coerce_map_from(Zp(17,40))
            True
        """
        #if isistance(R, pAdicRingLazy) and R.prime() == self.prime():
        #    return True
        if isinstance(R, pAdicRingCappedRelative) and R.prime() == self.prime():
            return True
        if isinstance(R, pAdicRingCappedAbsolute) and R.prime() == self.prime():
            if R.precision_cap() < self.precision_cap():
                return True
            elif (R.precision_cap() == self.precision_cap() and
                  self._printer.richcmp_modes(R._printer, op_LE)):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: ZpCA(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Ring with capped absolute precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

class pAdicRingFloatingPoint(pAdicRingBaseGeneric, pAdicFloatingPointRingGeneric):
    r"""
    An implementation of the `p`-adic integers with floating point
    precision.
    """
    def __init__(self, p, prec, print_mode, names):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options.
        - ``names`` -- how to print the prime.

        EXAMPLES::

            sage: R = ZpFP(next_prime(10^60)) #indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicRingFloatingPoint_with_category'>

        TESTS::

            sage: R = ZpFP(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^10)], max_runs = 2^12, skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = ZpFP(3, 1)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^3)])

            sage: R = ZpFP(3, 2)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^6)], skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = ZpFP(next_prime(10^60))
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^3)], max_runs = 2^5, skip='_test_log') # long time
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicFloatingPointElement)

    def _coerce_map_from_(self, R):
        """
        Returns ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = ZpFP(17)
            sage: K(1) + 1 #indirect doctest
            2
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
            sage: K.has_coerce_map_from(Zp(17,10))
            False
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
        """
        if isinstance(R, pAdicRingFloatingPoint) and R.prime() == self.prime():
            if R.precision_cap() > self.precision_cap():
                return True
            elif R.precision_cap() == self.precision_cap() and self._printer.richcmp_modes(R._printer, op_LE):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: ZpFP(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Ring with floating precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

class pAdicRingFixedMod(pAdicRingBaseGeneric, pAdicFixedModRingGeneric):
    r"""
    An implementation of the `p`-adic integers using fixed modulus.
    """
    def __init__(self, p, prec, print_mode, names):
        """
        Initialization

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options.
        - ``names`` -- how to print the prime.

        EXAMPLES::

            sage: R = ZpFM(next_prime(10^60)) #indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicRingFixedMod_with_category'>

        TESTS::

            sage: R = ZpFM(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^10)], max_runs = 2^12, skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = ZpFM(3, 1)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^3)])

            sage: R = ZpFM(3, 2)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^6)], skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = ZpFM(next_prime(10^60))
            sage: TestSuite(R).run(skip='_test_log')
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^4)], max_runs = 2^6, skip='_test_log') # long time
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])

        Fraction fields work after :trac:`23510`::

            sage: R = ZpFM(5)
            sage: K = R.fraction_field(); K
            5-adic Field with floating precision 20
            sage: K(R(90))
            3*5 + 3*5^2
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicFixedModElement)

    def _coerce_map_from_(self, R):
        """
        Returns ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = ZpFM(17)
            sage: K(1) + 1 #indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Zp(7))
            False
            sage: K.has_coerce_map_from(ZpFM(17,40))
            True
            sage: K.has_coerce_map_from(ZpFM(17,10))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
        """
        #if isistance(R, pAdicRingLazy) and R.prime() == self.prime():
        #    return True
        if isinstance(R, pAdicRingFixedMod) and R.prime() == self.prime():
            if R.precision_cap() > self.precision_cap():
                return True
            elif (R.precision_cap() == self.precision_cap() and
                  self._printer.richcmp_modes(R._printer, op_LE)):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: ZpFM(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Ring of fixed modulus 7^20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

class pAdicFieldCappedRelative(pAdicFieldBaseGeneric, pAdicCappedRelativeFieldGeneric):
    r"""
    An implementation of `p`-adic fields with capped relative precision.

    EXAMPLES::

        sage: K = Qp(17, 1000000) #indirect doctest
        sage: K = Qp(101) #indirect doctest

    """

    def __init__(self, p, prec, print_mode, names):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options.
        - ``names`` -- how to print the prime.

        EXAMPLES::

            sage: K = Qp(next_prime(10^60)) # indirect doctest
            sage: type(K)
            <class 'sage.rings.padics.padic_base_leaves.pAdicFieldCappedRelative_with_category'>

        TESTS::

            sage: R = Qp(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^10)], max_runs = 2^12, skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = Qp(3, 1)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^6)], skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = Qp(3, 2)
            sage: TestSuite(R).run(elements=[R.random_element() for i in range(3^9)], skip="_test_metric") # long time
            sage: R._test_metric(elements=[R.random_element() for i in range(3^3)])

            sage: R = Qp(next_prime(10^60))
            sage: TestSuite(R).run(skip='_test_log')
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^3)], max_runs = 2^5, skip='_test_log') # long time
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])
        """
        pAdicFieldBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedRelativeElement)

    def _coerce_map_from_(self, R):
        """
        Returns ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = Qp(17)
            sage: K(1) + 1 #indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            True
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Qp(17,40))
            False
            sage: K.has_coerce_map_from(Qp(17,10))
            True
            sage: K.has_coerce_map_from(Zp(17,40))
            True

        """
        #if isinstance(R, pAdicRingLazy) or isinstance(R, pAdicFieldLazy) and R.prime() == self.prime():
        #    return True
        if isinstance(R, (pAdicRingCappedRelative, pAdicRingCappedAbsolute)) and R.prime() == self.prime():
            return True
        if isinstance(R, pAdicFieldCappedRelative) and R.prime() == self.prime():
            if R.precision_cap() < self.precision_cap():
                return True
            elif (R.precision_cap() == self.precision_cap() and
                  self._printer.richcmp_modes(R._printer, op_LE)):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: Qp(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Field with capped relative precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

    def random_element(self, algorithm='default'):
        r"""
        Returns a random element of ``self``, optionally using the ``algorithm``
        argument to decide how it generates the element. Algorithms currently
        implemented:

        - default: Choose an integer `k` using the standard
          distribution on the integers.  Then choose an integer `a`
          uniformly in the range `0 \le a < p^N` where `N` is the
          precision cap of ``self``.  Return ``self(p^k * a, absprec =
          k + self.precision_cap())``.

        EXAMPLES::

            sage: Qp(17,6).random_element()
            15*17^-8 + 10*17^-7 + 3*17^-6 + 2*17^-5 + 11*17^-4 + 6*17^-3 + O(17^-2)
        """
        if (algorithm == 'default'):
            k = ZZ.random_element()
            a = ZZ.random_element(self.prime()**self.precision_cap())
            return self(self.prime()**k * a, absprec = k + self.precision_cap())
        else:
            raise NotImplementedError("Don't know %s algorithm"%algorithm)

class pAdicFieldFloatingPoint(pAdicFieldBaseGeneric, pAdicFloatingPointFieldGeneric):
    r"""
    An implementation of the `p`-adic rationals with floating point
    precision.
    """
    def __init__(self, p, prec, print_mode, names):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options.
        - ``names`` -- how to print the prime.

        EXAMPLES::

            sage: R = QpFP(next_prime(10^60)) #indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicFieldFloatingPoint_with_category'>

        TESTS::

            sage: R = QpFP(2)
            sage: TestSuite(R).run()
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^10)], max_runs = 2^12, skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = QpFP(3, 1)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^3)])

            sage: R = QpFP(3, 2)
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(3^6)], skip='_test_metric') # long time
            sage: R._test_metric(elements = [R.random_element() for i in range(2^3)]) # long time

            sage: R = QpFP(next_prime(10^60))
            sage: TestSuite(R).run(skip='_test_log')
            sage: TestSuite(R).run(elements = [R.random_element() for i in range(2^3)], max_runs = 2^5, skip='_test_log') # long time
            sage: R._test_log(max_runs=2, elements=[R.random_element() for i in range(4)])
        """
        pAdicFieldBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicFloatingPointElement)

    def _coerce_map_from_(self, R):
        """
        Returns ``True`` if there is a coerce map from ``R`` to ``self``.

        EXAMPLES::

            sage: K = QpFP(17)
            sage: K(1) + 1 #indirect doctest
            2
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            True
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
            sage: K.has_coerce_map_from(Qp(17,10))
            False
            sage: K.has_coerce_map_from(ZpFP(17))
            True
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
        """
        if isinstance(R, (pAdicRingFixedMod, pAdicRingFloatingPoint, pAdicFieldFloatingPoint)) and R.prime() == self.prime():
            if R.precision_cap() > self.precision_cap():
                return True
            elif R.precision_cap() == self.precision_cap() and self._printer.richcmp_modes(R._printer, op_LE):
                return True

    def _convert_map_from_(self, R):
        """
        Finds conversion maps from R to this ring.

        EXAMPLES::

            sage: QpFP(7).convert_map_from(Zmod(343))
            Lifting morphism:
              From: Ring of integers modulo 343
              To:   7-adic Field with floating precision 20
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if isinstance(R, IntegerModRing_generic):
            N = R.cardinality()
            p = self.prime()
            n = N.exact_log(p)
            if N == p**n:
                from sage.rings.padics.padic_generic import ResidueLiftingMap
                return ResidueLiftingMap._create_(R, self)

    def _repr_(self, do_latex=False):
        r"""
        Print representation.

        EXAMPLES::

            sage: K = QpFP(17); K #indirect doctest
            17-adic Field with floating precision 20
            sage: latex(K)
            \QQ_{17}
        """
        if do_latex:
            return "\\QQ_{%s}" % self.prime()
        return "%s-adic Field with floating precision %s"%(self.prime(), self.precision_cap())



# Lattice precision
###################

# Maybe the next class should go to sage.rings.padics.generic_nodes but I 
# don't understand quite well the structure of all classes in this directory
class pAdicLatticeGeneric(pAdicGeneric):
    def __init__(self, p, prec, label=None, proof=False):
        if proof:
            raise NotImplementedError("p-adic with *proved* lattice precision not implemented yet")
        if label is None:
            self._label = None
        else:
            self._label = str(label)
        self._prec_cap = prec
        self._precision = PrecisionLattice(p, label)
        # We do not use the standard attribute element_class 
        # because we need to be careful with precision
        # Instead we implement _element_constructor_ (cf below)
        self._element_class = pAdicLatticeElement

    def _prec_type(self):
        """
        Return the precision handling type

        EXAMPLES::

            sage: ZpLP(5)._prec_type()
            'lattice'
        """
        return 'lattice'

    def precision(self):
        """
        Return the lattice precision object attached to this parent

        EXAMPLES::

            sage: R = ZpLP(5,label='precision')
            sage: R.precision()
            Precision Lattice on 0 object (label: precision)

            sage: x = R(1,10); y = R(1,5)
            sage: R.precision()
            Precision Lattice on 2 objects (label: precision)

        .. SEEALSO::

            :class:`sage.rings.padics.lattice_precision.PrecisionLattice`
        """

        return self._precision

    def label(self):
        """
        Return the label of this parent

        NOTE:

        Labels can be used to distinguish between parents with
        the same defining data.

        They are useful in the lattice precision framework in order
        to limit the size of the lattice modeling the precision (which
        is roughly the number of elements having this parent).

        Elements of a parent with some label do not coerse to aparent 
        with a different label. However conversions are allowed.

        EXAMPLES:

            sage: R = ZpLP(5)
            sage: R.label()  # no label by default

            sage: R = ZpLP(5, label='mylabel')
            sage: R.label()
            'mylabel'

        Labels are typically useful to isolate computations.
        For example, assume that we first want to do some calculations
        with matrices::

            sage: R = ZpLP(5, label='matrices')
            sage: M = random_matrix(R, 4, 4)
            sage: M.determinant()     # somewhat random

        Now, if we want to do another unrelated computations, we can
        use a different label::

            sage: R = ZpLP(5, label='polynomials')
            sage: S.<x> = PolynomialRing(R)
            sage: P = (x-1)*(x-2)*(x-3)*(x-4)*(x-5)

        If we have not used labels, the software would have modeled the
        precision on the matrices and on the polynomials using the same
        lattice (resulting then in manipulating a lattice of higher
        dimension, having then a significant impact on performances).
        """
        return self._label

    def _element_constructor_(self, x, prec=None):
        """
        Create an element of this parent

        INPUT:

        - ``x``: the datum from which the element is created

        - ``prec`` -- an integer or ``None`` (the default); the
          absolute precision of the created element

        NOTE:

        This function tries to be sharp on precision as much as
        possible.
        For instance, if the datum ``x`` is itself an element of the
        same parent, the software remembers that the created element
        is actually equal to ``x`` (at infinite precision)::

            sage: R = ZpLP(2, prec=50)
            sage: x = R(1,10); x
            1 + O(2^10)
            sage: y = R(x)   # indirect doctest
            sage: y
            1 + O(2^10)
            sage: x - y
            O(2^50)
        """
        # We first try the __copy__ method which is sharp on precision
        try:
            if prec is None:
                return x.__copy__(parent=self)
            else:
                return x.__copy__(parent=self).add_bigoh(prec)
        except (TypeError, ValueError, AttributeError):
            pass
        return self._element_class(self, x, prec)

    def convert_multiple(self, *elts):
        """
        Convert a list of elements to this parent

        NOTE:

        This function tries to be sharp on precision as much as
        possible.
        In particular, if the precision of the input elements are 
        handled by a lattice, diffused digits of precision are
        preserved during the conversion.

        EXAMPLES::

            sage: R = ZpLP(2)
            sage: x = R(1,10); y = R(1,5)
            sage: x,y = x+y, x-y

        Remark that the pair `(x,y)` has diffused digits of precision::

            sage: x
            2 + O(2^5)
            sage: y
            O(2^5)
            sage: x + y
            2 + O(2^11)

            sage: R.precision().number_of_diffused_digits([x,y])
            6

        As a consequence, if we convert ``x`` and ``y`` separately, we
        loose some precision::

            sage: R2 = ZpLP(2, label='copy')
            sage: x2 = R2(x); y2 = R2(y)
            sage: x2
            2 + O(2^5)
            sage: y2
            O(2^5)
            sage: x2 + y2
            2 + O(2^5)

            sage: R2.precision().number_of_diffused_digits([x2,y2])
            0

        On the other hand, this issue dissapears when we use multiple
        conversion::

            sage: x2,y2 = R2.convert_multiple(x,y)
            sage: x2 + y2
            2 + O(2^11)

            sage: R2.precision().number_of_diffused_digits([x2,y2])
            6
        """
        p = self.prime()

        # We sort elements by precision lattice
        elt_by_prec = { }
        elt_other = [ ]
        indices = { }
        for i in range(len(elts)):
            x = elts[i]; idx = id(x)
            if indices.has_key(idx):
                indices[idx].append(i)
            else:
                indices[idx] = [i]
            if isinstance(x, pAdicLatticeElement):
                prec = x.parent().precision()
                if prec.prime() != p:
                    raise TypeError("conversion between different p-adic rings not supported")
                if elt_by_prec.has_key(prec):
                    elt_by_prec[prec].append(x)
                else:
                    elt_by_prec[prec] = [x]
            else:
                elt_other.append(x)

        # We create the elements
        ans = len(elts)*[None]
        selfprec = self._precision
        # First the elements with precision lattice
        for (prec, L) in elt_by_prec.iteritems():
            if prec is selfprec:
                # Here, we use the __copy__ method in order
                # to be sharp on precision
                for x in L:
                    y = x.__copy__(parent=self)
                    for i in indices[id(x)]:
                        ans[i] = y
            else:
                lattice = prec.precision_lattice(L)
                for j in range(len(L)):
                    x = L[j]; dx = [ ]
                    for i in range(j):
                        dx.append([L[i], lattice[i,j]])
                    prec = lattice[j,j].valuation(p)
                    y = self._element_class(self, x.value(), prec, dx=dx, dx_mode='values', check=False, reduce=False)
                    for i in indices[id(x)]:
                        ans[i] = y
                    L[j] = y
        # Now the other elements
        for x in elt_other:
            y = self._element_class(self, x)
            for i in indices[id(x)]:
                ans[i] = y

        # We return the created elements
        return ans



class pAdicRingLattice(pAdicLatticeGeneric, pAdicRingBaseGeneric):
    """
    An implementation of the `p`-adic integers with lattice precision
    """
    def __init__(self, p, prec, print_mode, names, label=None, proof=False):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options
        - ``names`` -- how to print the prime
        - ``label`` -- the label of this ring

        EXAMPLES::

            sage: R = ZpLP(next_prime(10^60)) # indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicRingLattice_with_category'>

            sage: R = ZpLP(2,label='init') # indirect doctest
            sage: R
            2-adic Ring with lattice precision (label: init)

        .. SEEALSO::

            :meth:`label`
        """
        pAdicLatticeGeneric.__init__(self, p, prec, label, proof)
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, None)

    def _repr_(self, do_latex=False):
        """
        Return a representation of this parent

        EXAMPLES::

            sage: R = ZpLP(2); R   # indirect doctest
            2-adic Ring with lattice precision
            sage: latex(R)
            \mathbb Z_{2}
        """
        if do_latex:
            if self._label is not None:
                return "\\verb'%s' (\simeq \\mathbb Z_{%s})" % (self._label, self.prime())
            else:
                return "\\mathbb Z_{%s}" % self.prime()
        else:
            if self._label is not None:
                return "%s-adic Ring with lattice precision (label: %s)" % (self.prime(), self._label)
            else:
                return "%s-adic Ring with lattice precision" % self.prime()

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to this ring

        EXAMPLES::

        """
        if isinstance(R, pAdicRingLattice) and R.precision() is self.precision():
            return True

    def random_element(self, prec=None):
        """
        Return a random element of this ring

        INPUT:

        - ``prec`` -- an integer or ``None`` (the default): the
          absolute precision of the generated random element

        EXAMPLES::

            sage: R = ZpLP(2)
            sage: R.random_element()    # random
            2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^10 + 2^11 + 2^14 + 2^15 + 2^16 + 2^17 + 2^18 + 2^19 + O(2^20)

            sage: R.random_element(prec=10)    # random
            1 + 2^3 + 2^4 + 2^7 + O(2^10)
        """
        if prec is None:
            prec = self._prec_cap
        p = self.prime()
        x = ZZ.random_element(p**prec)
        return self._element_class(self, x, prec=prec)

    def integer_ring(self, print_mode=None):
        """
        Return the integer ring of this ring

        INPUT:

        - ``print_mode`` -- the priting mode of the returned ring

        EXAMPLES::

            sage: R = ZpLP(2); R
            2-adic Ring with lattice precision
            sage: R.integer_ring()
            2-adic Ring with lattice precision
            sage: R.integer_ring() is R
            True

            sage: R2 = R.integer_ring(print_mode='terse'); R2
            2-adic Ring with lattice precision
            sage: R2 is R
            False
            sage: x = R(121,10); x
            1 + 2^3 + 2^4 + 2^5 + 2^6 + O(2^10)
            sage: x2 = R2(x); x2
            121 + O(2^10)

        Labels are kept unchanged by this function::

            sage: R = ZpLP(2, label='test'); R
            2-adic Ring with lattice precision (label: test)
            sage: R.integer_ring()
            2-adic Ring with lattice precision (label: test)

        TESTS::

            sage: R.integer_ring(print_mode='series') is R
            True
        """
        if print_mode is None:
            return self
        from sage.rings.padics.factory import Zp
        return Zp(self.prime(), self.precision_cap(), 'lattice', print_mode=self._modified_print_mode(print_mode), 
                  names=self._uniformizer_print(), label=self._label)

    def fraction_field(self, print_mode=None):
        """
        Return the fraction field of this ring

        INPUT:

        - ``print_mode`` -- the priting mode of the returned ring

        EXAMPLES::

            sage: R = ZpLP(2); R
            2-adic Ring with lattice precision
            sage: K = R.fraction_field(); K
            2-adic Field with lattice precision

            sage: K2 = R.fraction_field(print_mode='terse'); K2
            2-adic Field with lattice precision
            sage: K2 is K
            False
            sage: x = R(121,10); x
            1 + 2^3 + 2^4 + 2^5 + 2^6 + O(2^10)
            sage: x2 = K2(x); x2
            121 + O(2^10)
 
        Labels are kept unchanged by this function::

            sage: R = ZpLP(2, label='test'); R
            2-adic Ring with lattice precision (label: test)
            sage: R.fraction_field()
            2-adic Field with lattice precision (label: test)
        """
        from sage.rings.padics.factory import Qp
        return Qp(self.prime(), self.precision_cap(), 'lattice', print_mode=self._modified_print_mode(print_mode), 
                  names=self._uniformizer_print(), label=self._label)


class pAdicFieldLattice(pAdicLatticeGeneric, pAdicFieldBaseGeneric):
    """
    An implementation of the `p`-adic numbers with lattice precision
    """
    def __init__(self, p, prec, print_mode, names, label=None, proof=False):
        """
        Initialization.

        INPUT:

        - ``p`` -- prime
        - ``prec`` -- precision cap
        - ``print_mode`` -- dictionary with print options
        - ``names`` -- how to print the prime
        - ``label`` -- the label of this ring

        EXAMPLES::

            sage: R = QpLP(next_prime(10^60)) # indirect doctest
            sage: type(R)
            <class 'sage.rings.padics.padic_base_leaves.pAdicFieldLattice_with_category'>

            sage: R = QpLP(2,label='init') # indirect doctest
            sage: R
            2-adic Field with lattice precision (label: init)

        .. SEEALSO::

            :meth:`label`
        """
        pAdicLatticeGeneric.__init__(self, p, prec, label, proof)
        pAdicFieldBaseGeneric.__init__(self, p, prec, print_mode, names, None)

    def _repr_(self, do_latex=False):
        """
        Return a representation of this parent

        EXAMPLES::

            sage: K = QpLP(2); K   # indirect doctest
            2-adic Field with lattice precision
            sage: latex(K)
            \mathbb Q_{2}
        """
        if do_latex:
            if self._label is not None:
                return "\\verb'%s' (\simeq \\mathbb Q_{%s})" % (self._label, self.prime())
            else:
                return "\\mathbb Q_{%s}" % self.prime()
        else:
            if self._label is not None:
                return "%s-adic Field with lattice precision (label: %s)" % (self.prime(), self._label)
            else:
                return "%s-adic Field with lattice precision" % self.prime()

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coerce map from ``R`` to this ring

        EXAMPLES::

        """

        if isinstance(R, (pAdicRingLattice, pAdicFieldLattice)) and R.precision() is self.precision():
            return True

    def random_element(self, prec=None): # integral=False
        """
        Return a random element of this ring

        INPUT:

        - ``prec`` -- an integer or ``None`` (the default): the
          absolute precision of the generated random element

        EXAMPLES::

            sage: K = QpLP(2)
            sage: K.random_element()    # random
            2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^10 + 2^11 + 2^14 + 2^15 + 2^16 + 2^17 + 2^18 + 2^19 + O(2^20)

            sage: K.random_element(prec=10)    # random
            1 + 2^3 + 2^4 + 2^7 + O(2^10)
        """
        # TODO: do not pick only among integers
        if prec is None:
            prec = self._prec_cap
        p = self.prime()
        x = ZZ.random_element(p**prec)
        return self._element_class(self, x, prec)

    def integer_ring(self, print_mode=None):
        """
        Return the integer ring of this field

        INPUT:

        - ``print_mode`` -- the priting mode of the returned ring

        EXAMPLES::

            sage: K = QpLP(2); K
            2-adic Field with lattice precision
            sage: R = K.integer_ring(); R
            2-adic Ring with lattice precision

            sage: R2 = K.integer_ring(print_mode='terse'); R2
            2-adic Ring with lattice precision
            sage: R2 is R
            False
            sage: x = R(121,10); x
            1 + 2^3 + 2^4 + 2^5 + 2^6 + O(2^10)
            sage: x2 = R2(x); x2
            121 + O(2^10)

        Labels are kept unchanged by this function::

            sage: K = QpLP(2, label='test'); K
            2-adic Field with lattice precision (label: test)
            sage: K.integer_ring()
            2-adic Ring with lattice precision (label: test)

        """
        from sage.rings.padics.factory import Zp
        return Zp(self.prime(), self.precision_cap(), 'lattice', print_mode=self._modified_print_mode(print_mode), 
                  names=self._uniformizer_print(), label=self._label)

    def fraction_field(self, print_mode=None):
        """
        Return the fraction field of this ring

        INPUT:

        - ``print_mode`` -- the priting mode of the returned ring

        EXAMPLES::

            sage: K = QpLP(2); K
            2-adic Field with lattice precision
            sage: K.fraction_field()
            2-adic Field with lattice precision
            sage: K.fraction_field() is K
            True

            sage: K2 = K.fraction_field(print_mode='terse'); K2
            2-adic Field with lattice precision
            sage: K2 is K
            False
            sage: x = K(121,10); x
            1 + 2^3 + 2^4 + 2^5 + 2^6 + O(2^10)
            sage: x2 = K2(x); x2
            121 + O(2^10)
 
        Labels are kept unchanged by this function::

            sage: K = QpLP(2, label='test'); K
            2-adic Field with lattice precision (label: test)
            sage: K.fraction_field()
            2-adic Field with lattice precision (label: test)

        TESTS::

            sage: K.fraction_field(print_mode='series') is K
            True
        """
        if print_mode is None:
            return self
        from sage.rings.padics.factory import Qp
        return Qp(self.prime(), self.precision_cap(), 'lattice', print_mode=self._modified_print_mode(print_mode), 
                  names=self._uniformizer_print(), label=self._label)
