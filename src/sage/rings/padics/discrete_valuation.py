r"""
Discrete valuations

This file defines abstract base classes for discrete (pseudo-)valuations.

AUTHORS:

- Julian Rueth (2013-03-16): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.morphism import Morphism
from sage.misc.abstract_method import abstract_method
from sage.structure.unique_representation import UniqueRepresentation

from sage.sets.set import Set
from sage.rings.all import QQ
from sage.rings.infinity import infinity
discrete_valuation_codomain = Set(QQ).union(Set([infinity]))

class DiscretePseudoValuation(Morphism):
    r"""
    Abstract base class for discrete pseudo-valuations, i.e., discrete
    valuations which might send more that just zero to infinity.

    INPUT:

    - ``domain`` -- an integral domain

    - ``check`` -- a boolean (default: ``True``), if ``True``, we check that
      ``domain`` is an integral domain.

    EXAMPLES::

        sage: ZZ.valuation(2) # indirect doctest
        2-adic valuation

    """
    @staticmethod
    def __classcall__(cls, domain, check=True):
        r"""
        Normalize such that ``check`` is not part of the key for the
        :class:`sage.structure.factory.UniqueRepresentation`.

        TESTS::

            sage: TrivialDiscretePseudoValuation(ZZ, check=True) is TrivialDiscretePseudoValuation(ZZ, check=False)
            True
            sage: TrivialDiscreteValuation(ZZ, check=True) is TrivialDiscreteValuation(ZZ, check=False)
            True

        """
        if check:
            from sage.categories.integral_domains import IntegralDomains
            if not domain in IntegralDomains():
                raise ValueError("domain of a discrete valuation must be an integral domain")
        return super(DiscretePseudoValuation, cls).__classcall__(cls, domain)

    def __init__(self, domain):
        r"""
        Initialization.

        TESTS::

            sage: from sage.rings.padics.discrete_valuation import DiscretePseudoValuation
            sage: isinstance(ZZ.valuation(2), DiscretePseudoValuation)
            True

        """
        from sage.categories.homset import Hom
        from sage.categories.sets_cat import Sets
        Morphism.__init__(self, Hom(domain, discrete_valuation_codomain))

    def is_equivalent(self, f, g):
        r"""
        Return whether ``f`` and ``g`` are equivalent modulo this valuation,
        i.e., whether their difference has positive valuation.

        INPUT:

        - ``f``, ``g`` -- elements in the domain of this valuation

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: v.is_equivalent(2,4)
            True
            sage: v.is_equivalent(2,3)
            False

        """
        return self(f-g)>0

    def _test_add(self, **options):
        r"""
        Helper method to test that this valuation behaves correctly with
        respect to additions.

        TESTS::

            sage: v = ZZ.valuation(2)
            sage: v._test_add()

        """
        tester = self._tester(**options)
        S = tester.some_elements(self.domain().some_elements())
        from sage.combinat.cartesian_product import CartesianProduct
        for x,y in tester.some_elements(CartesianProduct(S,S)):
            tester.assertGreaterEqual(self(x+y),min(self(x),self(y)))
            if self(x) != self(y):
                tester.assertEqual(self(x+y),min(self(x),self(y)))

    def _test_infty(self, **options):
        r"""
        Helper method to test that the consistency of the elements which are
        mapped to `\infty`.

        TESTS::

            sage: from sage.rings.padics.discrete_valuation import DiscretePseudoValuation
            sage: v = ZZ.valuation(2)
            sage: DiscretePseudoValuation._test_infty(v) # _test_infty is overwritten by DiscreteValuation

        """
        tester = self._tester(**options)
        for x in tester.some_elements(self.domain().some_elements()):
            if x.is_zero():
                tester.assertEqual(self(x), infinity)

    def _test_mul(self, **options):
        r"""
        Helper method to test that this valuation behaves correctly with
        respect to multiplications.

        TESTS::

            sage: v = ZZ.valuation(2)
            sage: v._test_mul()

        """
        tester = self._tester(**options)
        S = self.domain().some_elements()
        from sage.combinat.cartesian_product import CartesianProduct
        for x,y in tester.some_elements(CartesianProduct(S,S)):
            tester.assertEqual(self(x*y),self(x)+self(y))

    @abstract_method
    def value_group(self):
        r"""
        Return the :class:`discrete_value_group.DiscreteValueGroup` of this
        valuation.

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: v.value_group()
            DiscreteValueGroup(1)

        For a trivial discrete pseudo-valuation, this is not defined::

            sage: v = TrivialDiscretePseudoValuation(ZZ)
            sage: v.value_group()
            Traceback (most recent call last):
            ...
            NotImplementedError: the value group of a trivial discrete pseudo-valuation is not defined.

        """
        pass

    @abstract_method
    def uniformizer(self):
        r"""
        Return a uniformizer for this valuation.

        OUTPUT:

        An element in the domain of this valuation with positive valuation
        whose valuation generates the :meth:`value_group` of this valuation.

        EXAMPLES::

            sage: ZZ.valuation(2).uniformizer()
            2

        For trivial valuations, this is undefined::

            sage: TrivialDiscreteValuation(ZZ).uniformizer()
            Traceback (most recent call last):
            ...
            NotImplementedError: the uniformizer of a trivial discrete valuation is not defined.
            sage: TrivialDiscretePseudoValuation(ZZ).uniformizer()
            Traceback (most recent call last):
            ...
            NotImplementedError: the uniformizer of a trivial discrete valuation is not defined.

        """
        pass

class TrivialDiscretePseudoValuation(DiscretePseudoValuation, UniqueRepresentation):
    r"""
    The trivial discrete pseudo-valuation on ``R``, i.e., the valuation which
    maps every element of ``R`` to `\infty`.

    INPUT:

    - ``domain`` -- an integral domain

    - ``check`` -- an boolean (default: ``True``), if ``True``, we check that
      ``domain`` is an integral domain.

    EXAMPLES::

        sage: v = TrivialDiscretePseudoValuation(ZZ); v
        Trivial pseudo-valuation on Integer Ring

    TESTS::

        sage: TestSuite(v).run()

    """
    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: TrivialDiscretePseudoValuation(ZZ)
            Trivial pseudo-valuation on Integer Ring

        """
        return "Trivial pseudo-valuation on {0}".format(self.domain())

    def _call_(self, x):
        r"""
        Evaluation this valuation at ``x``.

        INPUT:

        - ``x`` -- an element of the domain of this valuation.

        EXAMPLES::

            sage: v = TrivialDiscretePseudoValuation(ZZ)
            sage: v(1)
            +Infinity

        """
        return infinity

    def value_group(self):
        r"""
        The value group of this trivial valuation is undefined, hence this
        raises a ``NotImplementedError``.

        EXAMPLES::

            sage: v = TrivialDiscretePseudoValuation(ZZ)
            sage: v.value_group()
            Traceback (most recent call last):
            ...
            NotImplementedError: the value group of a trivial discrete pseudo-valuation is not defined.

        """
        raise NotImplementedError("the value group of a trivial discrete pseudo-valuation is not defined.")

    def uniformizer(self):
        r"""
        The uniformizer of this trivial valuation is undefined, hence this
        raise a ``NotImplementedError``.

        EXAMPLES::

            sage: v = TrivialDiscretePseudoValuation(ZZ)
            sage: v.uniformizer()
            Traceback (most recent call last):
            ...
            NotImplementedError: the uniformizer of a trivial discrete valuation is not defined.

        """
        raise NotImplementedError("the uniformizer of a trivial discrete valuation is not defined.")

    def __reduce__(self):
        r"""
        Return the parameters passed to the constructor; used for pickling.

        EXAMPLES::

            sage: v = TrivialDiscreteValuation(ZZ)
            sage: v.__reduce__()
            (<class 'sage.rings.padics.discrete_valuation.TrivialDiscreteValuation'>, (Integer Ring,))

        """
        return TrivialDiscretePseudoValuation, (self.domain(),)

class DiscreteValuation(DiscretePseudoValuation):
    r"""
    Abstract base class for discrete valuations.

    INPUT:

    - ``domain`` -- an integral domain

    - ``check`` -- an boolean (default: ``True``), if ``True``, we check that
      ``domain`` is an integral domain.

    EXAMPLES::

        sage: Zp(2).valuation() # indirect doctest
        2-adic valuation

    """
    def _test_infty(self, **options):
        r"""
        Helper method to test that the consistency of the elements which are
        mapped to `\infty`.

        TESTS::

            sage: v = ZZ.valuation(2)
            sage: v._test_infty()

        """
        tester = self._tester(**options)
        for x in tester.some_elements(self.domain().some_elements()):
            tester.assertEqual(x.is_zero(), self(x) == infinity)

    def value_group(self):
        r"""
        Return the value group of this valuation, i.e., the subgroup of \QQ
        generated by the valuation of the :meth:`uniformizer`.

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: v.value_group()
            DiscreteValueGroup(1)

        """
        from discrete_value_group import DiscreteValueGroup
        return DiscreteValueGroup(self(self.uniformizer()))

class TrivialDiscreteValuation(DiscreteValuation, UniqueRepresentation):
    r"""
    The trivial discrete valuation on ``R``, i.e., the valuation that sends
    every non-zero element to zero, and zero to `\infty`.

    INPUT:

    - ``domain`` -- an integral domain

    - ``check`` -- an boolean (default: ``True``), if ``True``, we check that
      ``domain`` is an integral domain.

    EXAMPLES::

        sage: v = TrivialDiscreteValuation(ZZ); v
        Trivial valuation on Integer Ring

    TESTS::

        sage: TestSuite(v).run()

    """
    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: TrivialDiscreteValuation(ZZ)
            Trivial valuation on Integer Ring

        """
        return "Trivial valuation on {0}".format(self.domain())

    def _call_(self, x):
        r"""
        Evaluate this valuation at ``x``.

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        EXAMPLES::

            sage: v = TrivialDiscreteValuation(ZZ)
            sage: v(0)
            +Infinity
            sage: v(1)
            0
            sage: v(2)
            0
        """
        if x.is_zero():
            return infinity
        from sage.rings.all import ZZ
        return ZZ(0)

    def __reduce__(self):
        r"""
        Return the parameters passed to the constructor; used for pickling.

        EXAMPLES::

            sage: v = TrivialDiscreteValuation(ZZ)
            sage: v.__reduce__()
            (<class 'sage.rings.padics.discrete_valuation.TrivialDiscreteValuation'>, (Integer Ring,))

        """
        return TrivialDiscreteValuation, (self.domain(),)

    def uniformizer(self):
        r"""
        The uniformizer of this trivial valuation is undefined, hence this
        raise a ``NotImplementedError``.

        EXAMPLES::

            sage: v = TrivialDiscreteValuation(ZZ)
            sage: v.uniformizer()
            Traceback (most recent call last):
            ...
            NotImplementedError: the uniformizer of a trivial discrete valuation is not defined.

        """
        raise NotImplementedError("the uniformizer of a trivial discrete valuation is not defined.")

    def value_group(self):
        r"""
        Return the value group of this valuation, i.e., the subgroup of \QQ
        generated by `0`.

        EXAMPLES::

            sage: v = TrivialDiscreteValuation(ZZ)
            sage: v.value_group()
            DiscreteValueGroup(0)

        For the trivial ring, the value group would be undefined. However, the
        trivial ring is not considered a integral domain::

            sage: v = TrivialDiscreteValuation(ZZ.quo(1))
            Traceback (most recent call last):
            ...
            ValueError: domain of a discrete valuation must be an integral domain

        """
        from discrete_value_group import DiscreteValueGroup
        return DiscreteValueGroup(0)
