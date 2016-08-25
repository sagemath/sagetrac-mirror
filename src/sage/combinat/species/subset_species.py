"""
Subset Species
"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .species import GenericCombinatorialSpecies
from .set_species import SetSpecies
from .generating_series import _integers_from, factorial_stream
from .structure import GenericSpeciesStructure
from sage.rings.all import ZZ
from sage.misc.cachefunc import cached_function
from sage.combinat.species.misc import accept_size
from sage.structure.unique_representation import UniqueRepresentation

class SubsetSpeciesStructure(GenericSpeciesStructure):
    def __repr__(self):
        """
        EXAMPLES::

            sage: S = species.SubsetSpecies()
            sage: a = S.structures(["a","b","c"]).random_element(); a
            {}
        """
        s = GenericSpeciesStructure.__repr__(self)
        return "{"+s[1:-1]+"}"

    def canonical_label(self):
        """
        EXAMPLES::

            sage: P = species.SubsetSpecies()
            sage: S = P.structures(["a", "b", "c"])
            sage: [s.canonical_label() for s in S]
            [{}, {'a'}, {'a'}, {'a'}, {'a', 'b'}, {'a', 'b'}, {'a', 'b'}, {'a', 'b', 'c'}]
        """
        rng = range(1, len(self._list)+1)
        return self.__class__(self.parent(), self._labels, rng)


    def label_subset(self):
        """
        Returns a subset of the labels that "appear" in this
        structure.

        EXAMPLES::

            sage: P = species.SubsetSpecies()
            sage: S = P.structures(["a", "b", "c"])
            sage: [s.label_subset() for s in S]
            [[], ['a'], ['b'], ['c'], ['a', 'b'], ['a', 'c'], ['b', 'c'], ['a', 'b', 'c']]
        """
        return [self._relabel(i) for i in self._list]

    def transport(self, perm):
        """
        Returns the transport of this subset along the permutation perm.

        EXAMPLES::

            sage: F = species.SubsetSpecies()
            sage: a = F.structures(["a", "b", "c"])[5]; a
            {'a', 'c'}
            sage: p = PermutationGroupElement((1,2))
            sage: a.transport(p)
            {'b', 'c'}
            sage: p = PermutationGroupElement((1,3))
            sage: a.transport(p)
            {'a', 'c'}
        """
        l = sorted([perm(i) for i in self._list])
        return SubsetSpeciesStructure(self.parent(), self._labels, l)

    def automorphism_group(self):
        """
        Returns the group of permutations whose action on this subset leave
        it fixed.

        EXAMPLES::

            sage: F = species.SubsetSpecies()
            sage: a = F.structures([1,2,3,4])[6]; a
            {1, 3}
            sage: a.automorphism_group()
            Permutation Group with generators [(2,4), (1,3)]

        ::

            sage: [a.transport(g) for g in a.automorphism_group()]
            [{1, 3}, {1, 3}, {1, 3}, {1, 3}]
        """
        from sage.groups.all import SymmetricGroup, PermutationGroup
        from sage.misc.all import uniq
        a = SymmetricGroup(self._list)
        b = SymmetricGroup(self.complement()._list)
        return PermutationGroup(a.gens()+b.gens())

    def complement(self):
        """
        EXAMPLES::

            sage: F = species.SubsetSpecies()
            sage: a = F.structures(["a", "b", "c"])[5]; a
            {'a', 'c'}
            sage: a.complement()
            {'b'}
        """
        new_list = [i for i in range(1, len(self._labels)+1) if i not in self._list]
        return SubsetSpeciesStructure(self.parent(), self._labels, new_list)

class SubsetSpecies(GenericCombinatorialSpecies, UniqueRepresentation):
    @staticmethod
    @accept_size
    def __classcall__(cls, *args, **kwds):
        """
        EXAMPLES::

            sage: S = species.SubsetSpecies(); S
            Subset species
        """
        return super(SubsetSpecies, cls).__classcall__(cls, *args, **kwds)

    def __init__(self, min=None, max=None, weight=None):
        """
        Returns the species of subsets.

        EXAMPLES::

            sage: S = species.SubsetSpecies()
            sage: S.generating_series().coefficients(5)
            [1, 2, 2, 4/3, 2/3]
            sage: S.isotype_generating_series().coefficients(5)
            [1, 2, 3, 4, 5]

            sage: S = species.SubsetSpecies()
            sage: c = S.generating_series().coefficients(3)
            sage: S._check()
            True
            sage: S == loads(dumps(S))
            True
        """
        GenericCombinatorialSpecies.__init__(self, min=None, max=None, weight=None)
        self._name = "Subset species"

    _default_structure_class = SubsetSpeciesStructure

    def _structures(self, structure_class, labels):
        """
        EXAMPLES::

            sage: S = species.SubsetSpecies()
            sage: S.structures([1,2]).list()
            [{}, {1}, {2}, {1, 2}]
            sage: S.structures(['a','b']).list()
            [{}, {'a'}, {'b'}, {'a', 'b'}]
        """
        from sage.combinat.combination import Combinations
        for c in Combinations(range(1, len(labels)+1)):
            yield structure_class(self, labels, c)

    def _isotypes(self, structure_class, labels):
        """
        EXAMPLES::

            sage: S = species.SubsetSpecies()
            sage: S.isotypes([1,2]).list()
            [{}, {1}, {1, 2}]
            sage: S.isotypes(['a','b']).list()
            [{}, {'a'}, {'a', 'b'}]
        """
        for i in range(len(labels)+1):
            yield structure_class(self, labels, range(1, i+1))

    def _gs_iterator(self, base_ring):
        """
        The generating series for the species of subsets is
        `e^{2x}`.

        EXAMPLES::

            sage: S = species.SubsetSpecies()
            sage: S.generating_series().coefficients(5)
            [1, 2, 2, 4/3, 2/3]
        """
        for n in _integers_from(0):
            yield  base_ring(2)**n/base_ring(factorial_stream[n])

    def _itgs_iterator(self, base_ring):
        """
        The generating series for the species of subsets is
        `e^{2x}`.

        EXAMPLES::

            sage: S = species.SubsetSpecies()
            sage: S.isotype_generating_series().coefficients(5)
            [1, 2, 3, 4, 5]
        """
        for n in _integers_from(1):
            yield base_ring(n)

    def _cis(self, series_ring, base_ring):
        r"""
        The cycle index series for the species of subsets satisfies

        .. math::

             Z_{\mathfrak{p}} = Z_{\mathcal{E}} \cdot Z_{\mathcal{E}}


        EXAMPLES::

            sage: S = species.SubsetSpecies()
            sage: S.cycle_index_series().coefficients(5)
            [p[],
             2*p[1],
             2*p[1, 1] + p[2],
             4/3*p[1, 1, 1] + 2*p[2, 1] + 2/3*p[3],
             2/3*p[1, 1, 1, 1] + 2*p[2, 1, 1] + 1/2*p[2, 2] + 4/3*p[3, 1] + 1/2*p[4]]
        """
        ciset = SetSpecies().cycle_index_series(base_ring)
        res = ciset**2
        if self.is_weighted():
            res *= self._weight
        return res

#Backward compatibility
SubsetSpecies_class = SubsetSpecies
