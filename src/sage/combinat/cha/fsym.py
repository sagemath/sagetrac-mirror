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
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.graded_hopf_algebras_with_basis import \
    GradedHopfAlgebrasWithBasis
from sage.combinat.tableau import StandardTableaux

from tools.generic_basis import GenericBasis


class FreeSymmetricFunctions(UniqueRepresentation, Parent):

    def __init__(self, R):
        '''
        TESTS::

            sage: FreeSymmetricFunctions(QQ)
            The combinatorial Hopf algebra of Free Symmetric Functions over the Rational Field
        '''
        from sage.categories.rings import Rings
        assert(R in Rings()), '%s must be a ring' % R
        Parent.__init__(
            self, base=R,
            category=GradedHopfAlgebras(R).WithRealizations()
        )

    def _repr_(self):
        return "The combinatorial Hopf algebra of Free Symmetric" + \
            " Functions over the %s" % self.base_ring()

    @lazy_attribute
    def dual(self):
        FSymDual = LazyImport(
            "sage.combinat.cha.fsym_dual", "FreeSymmetricFunctionsDual")
        return FSymDual(self.base())

    Fs = Fundamental = LazyImport(
        "sage.combinat.cha._fsym.fundamental_basis", "Fundamental")
    Es = Elementary = LazyImport(
        "sage.combinat.cha._fsym.elementary_basis", "Elementary")

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [GradedHopfAlgebrasWithBasis(R),
                    GradedHopfAlgebras(R).Realizations()]

        class ParentMethods:

            def __init_extra__(self):
                self.build_morphisms()

            def build_morphisms(self):
                '''
                Define morphisms associated to the current basis
                '''

        class Base(GenericBasis):
            _prefix = "** TO DEFINE **"
            _basis_indices = StandardTableaux()

            def _latex_term(self, m):
                from sage.misc.latex import latex
                prefix = self.print_options()['latex_prefix']
                if len(m) == 0:
                    return "1"
                return prefix + \
                    "_{\\vcenter{\\hbox{\\scalebox{.5}\n{" + latex(m) + "}}}}"
