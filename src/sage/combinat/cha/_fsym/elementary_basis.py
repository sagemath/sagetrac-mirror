# -*- coding: utf-8 -*-
r"""
The elementary basis of FQSym Hopf algebra.
"""
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
from sage.combinat.cha.fsym import FreeSymmetricFunctions


class Elementary( FreeSymmetricFunctions.Bases.Base ):
    '''

    '''
    _prefix = "Es"

    def product_on_basis(self, t1, t2):
        '''
        TESTS::

            sage: Es = FSym(QQ).Es()
            sage: Tableau.set_pretty_repr('compact')
            sage: t1 = Tableau([[1, 2], [3]]); t1
            |3|
            |1|2|
            sage: t2 = Tableau([[1,3],[2],[4]]);t2
            |4|
            |2|
            |1|3|
            sage: Es(t1)*Es(t2)
            Es
              |7|
              |3|5|
              |1|2|4|6|
        '''
        if len(t1) == 0:
            return self.monomial(t2)
        n = t1.size()

        def new_tab():
            new_tab = []
            for i in range(min(len(t1), len(t2))):
                new_tab.append(t1[i] + [j + n for j in t2[i]])
            new_tab.extend(t1[len(t2):] if len(t1) > len(t2) else \
                map(lambda line: [j + n for j in line], t2[len(t1):]))
            return new_tab
        return self.monomial(self.basis().keys()(new_tab()))

    def build_morphisms(self):
        '''
        TODO:: morphism from Es <-> Fs???
        '''
