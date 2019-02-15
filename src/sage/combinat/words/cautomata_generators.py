# coding=utf8
r"""
DeterministicAutomaton generators

AUTHORS:

- Paul Mercat (2018)- I2M AMU Aix-Marseille Universite - initial version
- Dominique Benielli (2018) Labex Archimede - I2M -
  AMU Aix-Marseille Universite - Integration in -SageMath

"""

# *****************************************************************************
#       Copyright (C) 2018 Paul Mercat <paul.mercat@univ-amu.fr>
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

from sage.combinat.words.cautomata import DeterministicAutomaton
from sage.misc.prandom import randint, random


class DeterministicAutomatonGenerators(object):
    """
    This class permits to generate various usefull DetAutomata.

    EXAMPLES::

        #. DeterministicAutomaton recognizing every words over a given alphabet
        sage: dag.AnyWord(['a','b'])
        DeterministicAutomaton with 1 state and an alphabet of 2 letters

        #. DeterministicAutomaton recognizing every letters of a given alphabet
        sage: dag.AnyLetter(['a','b'])
        DeterministicAutomaton with 2 states and an alphabet of 2 letters

        #. DeterministicAutomaton whose language is a single word
        sage: dag.Word(['a','b','a'])
        DeterministicAutomaton with 4 states and an alphabet of 2 letters

        #. DeterministicAutomaton with a given alphabet, recognizing the empty word
        sage: dag.EmptyWord(['a','b'])
        DeterministicAutomaton with 1 state and an alphabet of 2 letters

        #. random DeterministicAutomaton
        sage: dag.Random()      # random
        DeterministicAutomaton with 244 states and an alphabet of 183 letters
    """

    def AnyLetter(self, A, A2=None):
        """
        Generate a DeterministicAutomaton recognizing every letter of the alphabet A.

        INPUT:

        - ``A`` -- the result recognize every letter of this alphabet

        - ``A2`` (default: ``None``) -- alphabet of the result (must contain A)

        OUTPUT:

        A :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: dag.AnyLetter(['a', 'b'])
            DeterministicAutomaton with 2 states and an alphabet of 2 letters

            sage: dag.AnyLetter(['a', 'b'], ['a', 0, 'b', 1])
            DeterministicAutomaton with 2 states and an alphabet of 4 letters

        TESTS::

            sage: dag.AnyLetter(['a', 'b'], [0,1])
            Traceback (most recent call last):
            ...
            ValueError: A label of a transition is not in the alphabet [0, 1]
        """
        return DeterministicAutomaton([(0, 1, i) for i in A], A=A2, i=0, final_states=[1])

    def AnyWord(self, A, A2=None):
        """
        Generate a DeterministicAutomaton recognizing every words over the alphabet A.

        INPUT:

        - ``A`` -- the result recognize every word over this alphabet

        - ``A2`` (default: ``None``) -- alphabet of the result (must contain A)

        OUTPUT:

        A :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: dag.AnyWord(['a', 'b'])
            DeterministicAutomaton with 1 state and an alphabet of 2 letters

            sage: dag.AnyWord(['a', 'b'], ['a', 0, 'b', 1])
            DeterministicAutomaton with 1 state and an alphabet of 4 letters

        TESTS::

            sage: dag.AnyWord(['a', 'b'], [0,1])
            Traceback (most recent call last):
            ...
            ValueError: A label of a transition is not in the alphabet [0, 1]
        """
        return DeterministicAutomaton([(0, 0, i) for i in A], A=A2, i=0, final_states=[0])

    def Empty(self, A):
        """
        Generate a DeterministicAutomaton recognizing the empty language over the alphabet A.

        INPUT:

        - ``A`` -- alphabet of the result

        OUTPUT:

        A :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: dag.Empty(['a', 'b'])
            DeterministicAutomaton with 0 state and an alphabet of 2 letters
        """
        return DeterministicAutomaton([], A=A)

    def EmptyWord(self, A):
        """
        Generate a DeterministicAutomaton recognizing the empty word and having the alphabet A.

        INPUT:

        - ``A`` -- alphabet of the result

        OUTPUT:

        A :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: dag.EmptyWord(['a', 'b'])
            DeterministicAutomaton with 1 state and an alphabet of 2 letters
        """
        return DeterministicAutomaton([], S=[0], i=0, final_states=[0], A=A)

    def Word(self, w, A=None):
        """
        Generate a DeterministicAutomaton recognizing the word w.

        INPUT:

        - ``w`` -- the result recognize this word

        - ``A`` (default: ``None``) -- alphabet of the result (must contain the letters of w)

        OUTPUT:

        A :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: dag.Word(['a', 'b'])
            DeterministicAutomaton with 3 states and an alphabet of 2 letters

            sage: dag.Word(['a', 'b'], ['a', 0, 'b', 1])
            DeterministicAutomaton with 3 states and an alphabet of 4 letters

        TESTS::

            sage: dag.Word(['a', 'b'], [0,1])
            Traceback (most recent call last):
            ...
            ValueError: A label of a transition is not in the alphabet [0, 1]
        """
        return DeterministicAutomaton(
            [(i, i+1, j) for i, j in enumerate(w)], A=A, i=0, final_states=[len(w)])

    def Random(self, n=None, A=None, density_edges=None, 
               density_finals=None, verb=False):
        """
        Generate a random DeterministicAutomaton.

        INPUT:

        - ``n`` - int (default: ``None``) -- the number of states 

        - ``A`` (default: ``None``) -- alphabet of the result

        - ``density_edges`` (default: ``None``) -- the density of the transitions among all possible transitions

        - ``density_finals`` (default: ``None``) -- the density of final states among all states

        - ``verb`` - bool (default: ``False``) -- print informations for debugging

        OUTPUT:

        A :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: dag.Random()      # random
            DeterministicAutomaton with 401 states and an alphabet of 836 letters

            sage: dag.Random(3, ['a','b'])
            DeterministicAutomaton with 3 states and an alphabet of 2 letters

        """
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
        return DeterministicAutomaton(L, A=A, S=range(n), i=randint(0,n-1), final_states=F)


# Easy access to the automaton generators:
dag = DeterministicAutomatonGenerators()
