# -*- coding: utf-8 -*-
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>,
#                          Rémi Maurice <maurice@univ-mlv.fr>.
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
from sage.misc.abstract_method import abstract_method
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.combinat.tableau import StandardTableaux

from tools.generic_basis import GenericBasis
from fsym import FreeSymmetricFunctions
from all import FQSym


class FreeSymmetricFunctionsDual(UniqueRepresentation, Parent):

    def __init__(self, R):
        '''

        '''
        from sage.categories.rings import Rings
        assert(R in Rings()), '%s must be a ring' % R
        Parent.__init__(
            self, base=R,
            category=GradedHopfAlgebras(R).Quotients().WithRealizations()
        )

    def dual(self):
        return FreeSymmetricFunctions(self.base())

    def ambient(self):
        return FQSym(self.base())

    def _repr_(self):
        return "The combinatorial Hopf algebra of Free Symmetric Functions" + \
            " over the %s" % self.base_ring()

    def a_realization(self):
        return self.Fundamental()

    Gs = Fundamental = LazyImport(
        "sage.combinat.cha._fsym.fundamental_dual_basis", "FundamentalDual")

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [ModulesWithBasis(R),
                    GradedHopfAlgebras(R).Realizations(),
                    GradedHopfAlgebras(R).Quotients()]

        class ParentMethods:

            def __init_extra__(self):
                self.build_morphisms()

            def build_morphisms(self):
                '''
                Define morphisms associated to the current basis
                '''

            @abstract_method
            def ambient(self):
                pass

            def lift(self, Gs_elem):
                return self.ambient()(Gs_elem)

            def retract(self, G_elem):
                return self(G_elem)

        class Base(GenericBasis):
            _prefix = "** TO DEFINE **"
            _basis_indices = StandardTableaux()
