# -*- coding: utf-8 -*-
r"""
The fundamental basis of FSym dual Hopf algebra.
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>.
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
from sage.combinat.cha.fsym_dual import FreeSymmetricFunctionsDual


class FundamentalDual(FreeSymmetricFunctionsDual.Bases.Base):
    '''
    TESTS::

        sage: Gs = FSymDual(QQ).Gs(); Gs
        The combinatorial Hopf algebra of Free Symmetric Functions over the Rational Field in the realization Fundamental
        sage: G = FQSym(QQ).G(); G
        The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field in the realization FundamentalDual
        sage: Gs([[1]])**5 == Gs(G.sum_of_monomials(
        ....:     [sig for sig in Permutations(5)]))
        True
    '''
    _prefix = "Gs"

    def ambient(self):
        return self.realization_of().ambient().G()

    def build_morphisms(self):
        G = self.ambient()
        # G -> Gs 
        G._module_morphism(
            lambda sigma: self(sigma.robinson_schensted()[0]),
            codomain=self
        ).register_as_coercion()
        # Gs -> G 
        self.module_morphism(
            lambda stdtab: G(stdtab.reading_word_permutation()),
            codomain=G
        ).register_as_coercion()
