# -*- coding: utf-8 -*-
"""
BWSym

AUTHOR:

 - Ali Chouria
 - Olivier Mallet
 - Jean-Baptiste Priez
"""
#*****************************************************************************
#       Copyright (C) 2013 Ali Chouria <ali.chouria@univ-rouen.fr>,
#                          Olivier Mallet <olivier.mallet@univ-rouen.fr>,
#                          Jean-Baptiste Priez <jbp@kerios.fr>,
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
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.lazy_import import LazyImport
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.combinat.set_partition_into_lists import SetPartitionIntoList, SetPartitionIntoLists

from tools.generic_basis import GenericBasis


class BWSym(UniqueRepresentation, Parent):
    r'''

    '''

    def __init__(self, R):
        '''
        '''
        from sage.categories.rings import Rings
        assert(R in Rings()), '%s must be a ring' % R
        Parent.__init__(
            self, base=R,
            category=GradedHopfAlgebras(R).WithRealizations()
        )

    def _repr_(self):
        return "The combinatorial Hopf algebra of BWSym Functions" + \
               " over the %s" % self.base_ring()

    def a_realization(self):
        return self.Fundamental()

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [GradedHopfAlgebras(R).Realizations(),
                    GradedHopfAlgebrasWithBasis(R).Realizations()]

        class ParentMethods:

            def build_morphisms(self):
                '''
                Define morphisms associated to the current basis
                '''

            def __init_extra__(self):
                self.build_morphisms()

        class Base(GenericBasis):
            _prefix = "** TO DEFINE **"
            _basis_indices = SetPartitionIntoLists()

    class Fundamental(Bases.Base):
        _prefix = "Phi"

        def product_on_basis(self, pf1, pf2):
            """

            """
            return self.monomial(self.basis().keys()(
                list(pf1) + map(
                    lambda li: map(
                        lambda i: i + pf1.size(),
                        list(li)
                    ), pf2)
            ))

        def coproduct_on_basis(self, pf):
            """
            TESTS::

                sage: from sage.combinat.set_partition_into_lists import SetPartitionIntoList
                sage: from sage.combinat.cha.bwsym import BWSym
                sage: Phi = BWSym(QQ).Phi()
                sage: spl = SetPartitionIntoList([[1,2],[3]])
                sage: Phi(spl).coproduct()
                Phi[] # Phi[(1, 2), (3,)] + Phi[(1,)] # Phi[(1, 2)] + Phi[(1, 2)] # Phi[(1,)] + Phi[(1, 2), (3,)] # Phi[]
                sage: spl = SetPartitionIntoList([[1,3],[2,4,6],[6]])
                sage: Phi(spl).coproduct()
                Phi[] # Phi[(1, 3), (5,), (2, 4, 6)] + Phi[(1,)] # Phi[(1, 3), (2, 4, 5)] + Phi[(1, 2)] # Phi[(1, 2, 3), (4,)] + Phi[(1, 2), (3,)] # Phi[(1, 2, 3)] + Phi[(1, 2, 3)] # Phi[(1, 2), (3,)] + Phi[(1, 2, 3), (4,)] # Phi[(1, 2)] + Phi[(1, 3), (2, 4, 5)] # Phi[(1,)] + Phi[(1, 3), (5,), (2, 4, 6)] # Phi[]
            """
            from sage.combinat.permutation import to_standard
            from sage.sets.set import Set

            def std(lli):
                acc_pos = []
                tmp = []
                for li in lli:
                    tmp += list(li)
                    acc_pos.append(len(tmp))
                tmp = to_standard(tmp)
                nlli = []
                old_pos = 0
                for pos in acc_pos:
                    nlli.append(tmp[old_pos:pos])
                    old_pos = pos
                return self.basis().keys()(nlli)
            S = Set(list(pf))
            return self.tensor_square().sum_of_monomials(
                (std(p), std(S.difference(p))) for p in S.subsets()
            )

    Phi = Fundamental
