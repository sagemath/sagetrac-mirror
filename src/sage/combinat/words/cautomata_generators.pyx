# coding=utf8
"""
DetAutomaton generators

AUTHORS:

- Paul Mercat (2018)- I2M AMU Aix-Marseille Universite - initial version
- Dominique Benielli (2018) Labex Archimede - I2M -
  AMU Aix-Marseille Universite - Integration in -SageMath

"""

# *****************************************************************************
#       Copyright (C) 2014 Paul Mercat <paul.mercat@univ-amu.fr>
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
# *****************************************************************************
from __future__ import print_function
from libc.stdlib cimport malloc, free
from sage.graphs.digraph import DiGraph

cimport sage.combinat.words.cautomata
from cysignals.signals cimport sig_on, sig_off, sig_check
from cpython cimport bool as c_bool 
from sage.combinat.words.cautomata import DetAutomaton
from sage.misc.prandom import randint, random


class DetAutomatonGenerators(object):
    def Random(self, n=None, A=None, density_edges=None, 
               density_finals=None, verb=False):
        if density_edges is None:
            density_edges = random()
        if density_finals is None:
            density_finals = random()
        if n is None:
            n = randint(2, 1000)
        if A is None:
            if random() < .5:
                A = list('abcdefghijklmnopqrstuvwxyz')[:randint(0, 25)]
            else:
                A = range(1, randint(1, 1000))
        if verb:
            print("Random automaton with %s states, density of leaving edges %s, density of final states %s and alphabet %s"%(n, density_edges, density_finals, A))
        L = []
        for i in range(n):
            for j in range(len(A)):
                if random() < density_edges:
                    L.append((i, randint(0, n-1), A[j]))
        if verb:
            print(L)
        F = []
        for i in range(n):
            if random() < density_finals:
                F.append(i)
        if verb:
            print("final states %s" % F)
        return DetAutomaton(L, i=randint(0, n-1), final_states=F)

 