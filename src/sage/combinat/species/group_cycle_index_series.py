# -*- coding: utf-8 -*-
"""
Group Cycle Indices

This file implements the group cycle indices of Henderson and Gainer-Dewar.

For a group `\\Gamma` and a ring `R`, a `\\Gamma`-cycle index over `R`
is an element of the free module over the ring of cycle index series
over `R` with basis `\\Gamma`.

In other words, a `\\Gamma`-cycle index series over `R` is a formal sum

.. MATH::

    F = \\sum_{\gamma \\in \Gamma} F[\gamma] \cdot \gamma,

where each coefficient `F[\\gamma]` is a cycle index series over `R`.
(By functoriality, if `F` is the `\\Gamma`-cycle index of a
`\\Gamma`-species, it must be a class function of `\\Gamma`.)

These objects are of interest because they can be used to enumerate `\\Gamma`-species;
they serve the same role in that theory as ordinary cycle indices do for classical
species.

Just like classical species, `\\Gamma`-species may be combined in ways
which correspond naturally to algebraic operations on their
`\\Gamma`-cycle indices.  Many standard operations on
`\\Gamma`-species are available, including addition, multiplication,
and differentiation.  The associated operations `+`, `\\cdot`, and
`\\prime` on `\\Gamma`-cycle indices are defined componentwise and are
implemented by this class.

Additionally, there is a natural way to compose `\\Gamma`-species,
corresponding to the construction of
:class:`~sage.combinat.species.composition_species.CompositionSpecies`.
The `\\Gamma`-cycle index of such a composition can be computed using
a variant of plethysm, but this depends on having access to all the
terms of both component `\\Gamma`-cycle indices; it cannot be computed
componentwise.  More details of this operation are available in the
documentation for :meth:`~sage.combinat.species.group_cycle_index_series.GroupCycleIndexSeries.composition`, which implements it in Sage.

AUTHORS:

- Andrew Gainer-Dewar (2013): initial version

EXAMPLES::

    sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
    sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(4))
    sage: loads(dumps(GCISR))
    Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field

.. SEEALSO::

    :mod:`Examples of Group Cycle Indices <~sage.combinat.species.group_cycle_index_series_library>`

.. TODO::

    Implement (optional?) optimizations using assumption that a
    `\\Gamma`-cycle index is a class function of `\\Gamma`.

"""
#*****************************************************************************
#       Copyright (C) 2013 Andrew Gainer-Dewar <andrew.gainer.dewar@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_function
from sage.rings.all import QQ, NN
from sage.combinat.free_module import CombinatorialFreeModule,CombinatorialFreeModuleElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.species.series_stream import (LinearCombinationStream, TermStream,
                                                 SumGeneratorStream, PowerStream)
from sage.combinat.species.generating_series import CycleIndexSeries, CycleIndexSeriesRing

class GroupCycleIndexSeriesRing(CombinatorialFreeModule, UniqueRepresentation):
    """
    Return the ring of group cycle index series.

    EXAMPLES::

        sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
        sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(2)); GCISR
        Ring of (Symmetric group of order 2! as a permutation group)-Cycle Index Series over Rational Field

    TESTS:

    We test to make sure the test suite runs::

        sage: TestSuite(GCISR).run()
    """
    def __init__(self, G, R=QQ):
        """
        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(4)); GCISR
            Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field
            sage: GCISR == loads(dumps(GCISR))
            True
        """
        from sage.categories.algebras_with_basis import AlgebrasWithBasis
        
        self._coeff_ring = R
        CISR = CycleIndexSeriesRing(R)
        self._cisr = CISR
        self._group = G

        CombinatorialFreeModule.__init__(self, CISR, G, element_class=GroupCycleIndexSeries,
                                         category = AlgebrasWithBasis(CISR), prefix='G')
    
    def product_on_basis(self, left, right):
        """
        Return the product of two basis elements ``left`` and ``right`` of ``self``.
        
        Multiplication of `\Gamma`-cycle indices is defined componentwise.
        That is, if `F` and `G` are two `\Gamma`-cycle indices, then `(F \cdot G) [\gamma] = F [\gamma] \cdot G [\gamma]`,
        where the multiplication on the right-hand side is ordinary multiplication of cycle indices.

        This is handled in Sage by defining multiplication on the basis of monomials induced
        by elements of `\Gamma`.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(4))
            sage: e = SymmetricGroup(4).identity()
            sage: t = SymmetricGroup(4)([4,3,2,1])
            sage: GCISR.product_on_basis(t,t) == GCISR(t)
            True
        
        """
        if left == right:
            return self.monomial(left)
        else:
            return self(0)

    def one(self):
        """
        Return the multiplicative identity element of this algebra.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(2))
            sage: GCISR.one()
            G[()] + G[(1,2)]

        """
        basis = self.basis()
        return self.sum(basis[k] for k in basis.keys())

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: GroupCycleIndexSeriesRing(SymmetricGroup(4))
            Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field
        """
        return "Ring of (%s)-Cycle Index Series over %s" %(self._group, self._coeff_ring)

    def group(self):
        """
        Returns the group of this :class:`GroupCycleIndexSeriesRing`.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: G = GroupCycleIndexSeriesRing(SymmetricGroup(4))
            sage: G.group()
            Symmetric group of order 4! as a permutation group
        """
        return self._group

    def recursive_series(self):
        """
        Returns a group cycle index series ring that can de defined
        recursively.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: S2 = SymmetricGroup(2)
            sage: GCISR = GroupCycleIndexSeriesRing(S2)
            sage: from sage.combinat.species.group_cycle_index_series_library import LinearOrderWithReversalGroupCycleIndex
            sage: L = LinearOrderWithReversalGroupCycleIndex()
            sage: X = GCISR(species.SingletonSpecies().cycle_index_series())
            sage: A = GCISR.recursive_series()
            sage: A.define(X*L(A))
            sage: A.quotient().isotype_generating_series().counts(8)
            [0, 1, 1, 2, 4, 10, 26, 76]
        """
        terms = [(g, self._cisr()) for g in self.basis().keys()]
        return self.sum_of_terms(terms, distinct=True)

class GroupCycleIndexSeries(CombinatorialFreeModuleElement):
    """
    EXAMPLES::

        sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
        sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(4))
        sage: GCISR.an_element()
        G[()] + 2*G[(3,4)] + 3*G[(2,3)] + G[(1,2,3,4)]
    """
    
    def quotient(self):
        r"""
        Return the quotient of this group cycle index.

        This is defined to be the ordinary cycle index `F / \Gamma` obtained from a
        `\Gamma`-cycle index `F` by:

        .. MATH::
            F / \Gamma = \\frac{1}{\\lvert \Gamma \\rvert} \sum_{\gamma \in \Gamma} F [\gamma].

        It is shown in [AGdiss]_ that, if `F` is the `\Gamma`-cycle
        index of a `\Gamma`-species, `F / \Gamma` is the ordinary
        cycle index of orbits of structures under the action of
        `\Gamma`.  This corresponds to the notion of quotient species
        introduced in [Bousquet]_.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: S4 = SymmetricGroup(4)
            sage: GCISR = GroupCycleIndexSeriesRing(S4)
            sage: GCISR.an_element()
            G[()] + 2*G[(3,4)] + 3*G[(2,3)] + G[(1,2,3,4)]
            sage: GCISR.an_element().quotient().coefficients(4)
            [7/24*p[], 0, 0, 0]

        REFERENCES:

        .. [AGdiss] Andrew Gainer. "`\Gamma`-species, quotients, and graph enumeration". Ph.D. diss. Brandeis University, 2012.
        .. [Bousquet] Michel Bousquet. "Especes de structures et applications du denombrement de cartes et de cactus planaires". Ph.D. diss. Universite du Quebec a Montreal, 1999.
           http://lacim.uqam.ca/publications_pdf/24.pdf

        """

        return 1/self.parent().group().cardinality() * sum(self.coefficients())

    def composition(self, y, check=False):
        r"""
        Return the group-cycle-index plethysm of ``self`` with ``y``.
        
        Plethysm of group cycle index series is defined by a sort of
        'mixing' operation in [Hend]_:

        .. MATH::
            (F \circ G) [\gamma] (p_{1}, p_{2}, p_{3}, \dots) =
            F [\gamma] \left( G [\gamma] (p_{1}, p_{2}, p_{3}, \dots),
            G [\gamma^{2}] (p_{2}, p_{4}, p_{6}, \dots), \dots \\right).

        It is shown in [Hend]_ that this operation on $\Gamma$-cycle
        indices corresponds to the 'composition' operation on
        $\Gamma$-species, a natural analogue of the
        :class:`~sage.combinat.species.composition_species.CompositionSpecies`
        operation on ordinary species.

        It is required that each of the components of `y` has zero
        constant term.  However, testing this may not be possible if
        `y` is of an exotic class.  Set ``check`` to ``False`` to
        suppress any testing.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: from sage.combinat.species.group_cycle_index_series_library import CyclicOrderWithReversalGroupCycleIndex
            sage: S2 = SymmetricGroup(2)
            sage: GCISR = GroupCycleIndexSeriesRing(S2)
            sage: E = sage.combinat.species.set_species.SetSpecies().cycle_index_series()
            sage: C = CyclicOrderWithReversalGroupCycleIndex()
            sage: example = (E*GCISR.one())(C)
            sage: example.quotient().generating_series().counts(8)
            [1, 1, 2, 5, 17, 73, 398, 2636]
            sage: x = oeis("A007868")         # optional -- internet
            sage: x.first_terms(8)            # optional -- internet
            (1, 1, 2, 5, 17, 73, 398, 2636)
            sage: x.name()                    # optional -- internet
            'Number of inverse pairs of elements in symmetric group S_n, or pairs of total orders on n nodes (average of A000085 and A000142).'            
            sage: example[S2.identity()].coefficients(4)
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
            sage: example[S2.gen()].coefficients(4)
            [p[], p[1], p[1, 1] + p[2], 2/3*p[1, 1, 1] + 2*p[2, 1] + 1/3*p[3]]

            sage: example[S2.gen()].generating_series().counts(8)
            [1, 1, 2, 4, 10, 26, 76, 232]
            sage: x = oeis("A000085")         # optional -- internet
            sage: x.first_terms(8)            # optional -- internet
            (1, 1, 2, 4, 10, 26, 76, 232)
            sage: x.name()                    # optional -- internet
            'Number of self-inverse permutations on n letters, also known as involutions; number of Young tableaux with n cells.'

        REFERENCES:

        .. [Hend] Anthony Henderson. "Species over a finite field". J. Algebraic Combin., 21(2):147-161, 2005.

        """
        assert self.parent() == y.parent()

        if check:
            for ycis in y.coefficients():
                assert ycis.coefficient(0) == 0

        # We define y_powers in this scope so that they can be reused
        # across all of the component series
        y_powers = {g: PowerStream(y[g]._stream) for g in self.parent().group()}

        # Map each stream to a composition stream
        def component_builder(g, series):
            return (g, series._new(self.CompositionStream, series._stream, y_powers, g))
        return self.map_item(component_builder)

    class CompositionStream(SumGeneratorStream):
        def __init__(self, stream, y_powers, group_element, **kwds):
            """
            A class for the coefficients of the composition (plethysm)
            of the group cycle index stream ``stream`` and the group
            cycle index stream corresponding to ``y_powers``.

            This stream is a :class:`SumGeneratorStream` whose
            component streams are given by :meth:`stream_iterator`.

            EXAMPLES::

                sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
                sage: from sage.combinat.species.group_cycle_index_series_library import CyclicOrderWithReversalGroupCycleIndex
                sage: from sage.combinat.species.series_stream import PowerStream
                sage: S2 = SymmetricGroup(2)
                sage: e = S2.gen()
                sage: GCISR = GroupCycleIndexSeriesRing(S2)
                sage: E = sage.combinat.species.set_species.SetSpecies().cycle_index_series()
                sage: C = CyclicOrderWithReversalGroupCycleIndex()
                sage: EG = E*GCISR.one()
                sage: y_powers = {g: PowerStream(C[g]._stream) for g in S2}
                sage: s = EG.CompositionStream(EG[e]._stream, y_powers, e, base_ring=E.base_ring())
                sage: [s[i] for i in range(4)]
                [p[], p[1], p[1, 1] + p[2], 2/3*p[1, 1, 1] + 2*p[2, 1] + 1/3*p[3]]
            """
            self._stream = stream
            self._y_powers = y_powers
            self._g = group_element
            super(GroupCycleIndexSeries.CompositionStream, self).__init__(self.stream_iterator(), **kwds)

        @staticmethod
        def example():
            """
            Returns an instance of
            :class:`GroupCycleIndexSeries.CompositionStream`.

            EXAMPLES::

                sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeries
                sage: e = GroupCycleIndexSeries.CompositionStream.example(); e
                <class 'sage.combinat.species.group_cycle_index_series.CompositionStream'>
            """
            from sage.groups.all import SymmetricGroup
            from sage.combinat.species.group_cycle_index_series_library import CyclicOrderWithReversalGroupCycleIndex
            from sage.combinat.all import species
            S2 = SymmetricGroup(2)
            GCISR = GroupCycleIndexSeriesRing(S2)
            E = species.SetSpecies().cycle_index_series()
            C = CyclicOrderWithReversalGroupCycleIndex()
            example = (E*GCISR.one())(C)
            return example[S2.gen()]._stream


        def stream_iterator(self):
            """
            A generator for the composition of each of the terms in
            ``self._stream`` with the group cycle index stream
            corresponding to ``self._y_powers``.

            EXAMPLES::

                sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeries
                sage: e = GroupCycleIndexSeries.CompositionStream.example()
                sage: outer = e._stream
                sage: inner = e._y_powers[e._g][0]
                sage: it = e.stream_iterator()
                sage: t0 = it.next()
                sage: outer[0]
                p[]
                sage: [t0[i] for i in range(4)]
                [p[], 0, 0, 0]
                sage: t1 = it.next() # outer[1](y)
                sage: outer[1]
                p[1]
                sage: [t1[i] for i in range(4)]
                [0, p[1], 1/2*p[1, 1] + 1/2*p[2], p[2, 1]]
                sage: [inner[i] for i in range(4)]
                [0, p[1], 1/2*p[1, 1] + 1/2*p[2], p[2, 1]]
            """
            for i in NN:
                term = self._stream[i]
                composed_terms = [(self.monomial_composer(p), c) for (p, c) in term if c != 0]
                yield LinearCombinationStream(composed_terms,
                                              base_ring=self._base_ring)

        def monomial_composer(self, partition):
            """
            Returns a stream which corresponds to composing a single
            monomial ``partition`` with the group cycle index series
            $y$.

            EXAMPLES::

                sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeries
                sage: e = GroupCycleIndexSeries.CompositionStream.example()
                sage: inner = e._y_powers[e._g][0]
                sage: p = Partition([2])
                sage: s = e.monomial_composer(p)
                sage: [s[i] for i in range(7)]
                [0, 0, p[2], 0, 1/2*p[2, 2] + 1/2*p[4], 0, 1/3*p[2, 2, 2] + 2/3*p[6]]
                sage: [e._y_powers[e._g**2][0][i] for i in range(4)]
                [0, p[1], 1/2*p[1, 1] + 1/2*p[2], 1/3*p[1, 1, 1] + 2/3*p[3]]
            """
            from sage.misc.misc_c import prod
            if not partition:
                return TermStream(n=0, value=self._base_ring.one(),
                                  base_ring=self._base_ring)
            exp = partition.to_exp_dict()
            streams = []
            for k, multiplicity in exp.items():
                s = self._y_powers[self._g**k][multiplicity - 1]   # s = (y[g^k])^(multiplicity)
                streams.append(CycleIndexSeries.StretchStream(s, k, base_ring=self._base_ring))
            return prod(streams)

    __call__ = composition

    def derivative(self):
        """
        Return the cycle-index derivative of ``self``.
        
        Differentiation of group cycle index series is defined termwise:

        .. MATH::
            (F')[\gamma] = (F [\gamma])'

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: S4 = SymmetricGroup(4)
            sage: GCISR = GroupCycleIndexSeriesRing(S4)
            sage: E = GCISR(species.SetSpecies().cycle_index_series())
            sage: E.derivative()[S4.an_element()].coefficients(6) == E[S4.an_element()].coefficients(6)
            True

        """
        return self.map_coefficients(lambda x: x.derivative())

    def restricted(self, min=None, max=None):
        """
        Return the restriction of ``self`` with coefficients starting at degree
        ``min`` and going up to, but not including, degree ``max``.
        If ``min`` is not specified, it is assumed to be zero. If ``max`` is not
        specified, it assumed to be infinity.

        This method simply calls :meth:`~sage.combinat.species.series.LazyPowerSeries.restricted` on each term
        of ``self``.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series_library import LinearOrderWithReversalGroupCycleIndex
            sage: L = LinearOrderWithReversalGroupCycleIndex()
            sage: e,t = L.parent().basis().keys()
            sage: L.restricted(min=2,max=4)[t].coefficients(5)
            [0, 0, p[2], p[2, 1], 0]

        """
        return self.map_coefficients(lambda x: x.restricted(min=min, max=max))

    def define(self, x):
        """
        Set ``self`` equal to ``x``, possibly recursively.

        EXAMPLES:

        Consider the `\mathfrak{S}_{2}`-species `\mathcal{A}` of rooted ordered trees,
        where the action of the nontrivial element of `\mathfrak{S}_{2}` inverts all the
        orderings on subtrees.

        Then we have the recursive functional equation

        .. MATH::
            \mathcal{A} = X \cdot \mathcal{L} (\mathcal{A})

        where `X` is the trivial singleton `\mathfrak{S}_{2}`-species and `\mathcal{L}` is the
        `\mathfrak{S}_{2}`-species of linear orderings with the reversal action
        described previously.

        To enumerate this species using Sage, we first need to set up an environment::
        
            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: S2 = SymmetricGroup(2)
            sage: GCISR = GroupCycleIndexSeriesRing(S2)
            sage: from sage.combinat.species.group_cycle_index_series_library import LinearOrderWithReversalGroupCycleIndex
            sage: L = LinearOrderWithReversalGroupCycleIndex()
            sage: X = GCISR(species.SingletonSpecies().cycle_index_series())

        We can then use ``define`` to set up `\mathcal{A}`::
        
            sage: A = GCISR.recursive_series()
            sage: A.define(X*L(A))
            sage: A.quotient().isotype_generating_series().counts(8)
            [0, 1, 1, 2, 4, 10, 26, 76]

        (Compare :oeis:`A007123`.)
        
        """        
        keys = self.parent().basis().keys()
        for key in keys:
            self[key].define(x[key])
