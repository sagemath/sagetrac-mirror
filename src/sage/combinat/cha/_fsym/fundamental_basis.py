# -*- coding: utf-8 -*-
r"""
The fundamental basis of FSym Hopf algebra.
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
from sage.combinat.cha.fsym import FreeSymmetricFunctions
from sage.combinat.cha.all import FQSym


class Fundamental(FreeSymmetricFunctions.Bases.Base):
    '''
    TODO:: tests
    '''
    _prefix = "Fs"

    def plactic_class(self, tab):
        from sage.combinat.tools import transitive_ideal
        from sage.combinat.permutation import Permutation

        def knuth_relations(w):
            li = []
            # bca -> bac with a < b <= c
            b = w[0]; c = w[1]; a = w[2]
            if a < b and b <= c:
                li.append([b, a, c])
            # bca <- bac with a < b <= c
            b = w[0]; a = w[1]; c = w[2]
            if a < b and b <= c:
                li.append([b, c, a])
            # acb -> cab with a <= b < c
            a = w[0]; c = w[1]; b = w[2]
            if a <= b and b < c:
                li.append([c, a, b])
            # acb <- cab with a <= b < c
            c = w[0]; a = w[1]; b = w[2];
            if a <= b and b < c:
                li.append([a, c, b])
            return li

        sigma = tab.reading_word_permutation()

        def parse(sigma):
            plactic_class = []
            for i in range(len(sigma) - 2):
                for trans in knuth_relations(sigma[i:i + 3]):
                    plactic_class.append(
                        Permutation(sigma[:i] + trans + sigma[i + 3:])
                    )
            return plactic_class
        return transitive_ideal(parse, sigma)

    def build_morphisms(self):
        # Fs -> F
        F = FQSym(self.base()).F()
        self.module_morphism(
            on_basis=lambda tab: F.sum_of_monomials(self.plactic_class(tab)),
            codomain=F
        ).register_as_coercion()

    def product_on_basis(self, stdtab1, stdtab2):
        '''
        TESTS::

            sage: Fs = FSym(QQ).Fs()
            sage: t = StandardTableau([[1, 2], [3]]);t
            [[1, 2], [3]]
            sage: t2 = StandardTableau([[1,2,3]]);t2
            [[1, 2, 3]]
            sage: Fs(t) * Fs(t2)
            Fs[[1, 2, 4, 5, 6], [3]] + Fs[[1, 2, 5, 6], [3], [4]] + Fs[[1, 2, 5, 6], [3, 4]] + Fs[[1, 2, 6], [3, 5], [4]]
            sage: t3 = StandardTableau([[1, 3, 4], [2]])
            sage: t4 = StandardTableau([[1, 3], [2]])
            sage: Fs(t3) * Fs(t4)
            Fs[[1, 3, 4], [2, 5], [6, 7]] + Fs[[1, 3, 4], [2, 5, 7], [6]] + Fs[[1, 3, 4], [2, 7], [5], [6]] + Fs[[1, 3, 4, 5], [2, 6, 7]] + Fs[[1, 3, 4, 5], [2, 7], [6]] + Fs[[1, 3, 4, 5, 7], [2], [6]] + Fs[[1, 3, 4, 5, 7], [2, 6]] + Fs[[1, 3, 4, 7], [2], [5], [6]] + Fs[[1, 3, 4, 7], [2, 5], [6]]

        TODO:: use taquin...
        '''
        from sage.misc.misc import uniq
        from sage.combinat.rsk import RSK
        F = FQSym(self.base()).F()
        return self.sum_of_monomials(uniq(map(
            lambda sig: RSK(sig)[0],
            (F(self(stdtab1)) * F(self(stdtab2))).\
                monomial_coefficients().keys()
        )))
