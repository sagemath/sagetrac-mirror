r"""
Motzkin Words

AUTHORS:

- Judith Braunsteiner (2017): Initial version
"""
#*****************************************************************************
#       Copyright (C) 2017 Judith Braunsteiner <judith.braunsteiner@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function, absolute_import

from sage.structure.list_clone import ClonableArray
from sage.combinat.backtrack import GenericBacktracker
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.rings.all import ZZ, QQ, Integer, infinity
from sage.arith.all import bernoulli, binomial, factorial

from six.moves import range
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.words.word import Word
from sage.structure.global_options import GlobalOptions
from sage.interfaces.all import maxima
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.libs.all import pari
from sage.misc.prandom import randint
from sage.misc.all import prod
from sage.structure.sage_object import SageObject
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.enumerated_sets import EnumeratedSets



class MotzkinWord(Element):
    r"""
    A Motzkin word.

    A Motzkin word is a sequence of up, down and horizontal steps such that
    every down stap has a corresponding up step preciding it.
    A Motzkin word of length `n` is a list with `k1` entries 1, `k2` entries
    -1 and `k3` entries 0 such that `k1 + k2 + k3 = n` and the first `i`
    entries als ways have at least as many 1s among them as -1.
    
    A Motzkin word may also be thought of as a lattice path in the
    `\mathbb{Z}^2` grid, starting at the origin `(0,0)`, and with steps in the
    North-East `NE = (1,1)`, South-East `SE = (1,-1)` and East `E = (1,0)`
    directions such that is does not pass the horizontal axis
    (x-axis, `y = 0`). A North-East step is represented by a 1 in the list, a
    South-Est step is represented by a -1 and an East step y a 0.

    A path representing a Motzkin word is called a Motzkin path.
  
    EXAMPLES::

    """

    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, mw):
        """
        Create a MotzkinWord.

        EXAMPLES::

            sage: MotzkinWord([1,1,0,-1,1,-1,0])
            [1,1,0,-1,1,-1,0]
        """
        mw = list(mw)
        n=len(mw)
        M = MotzkinWords(n)
        for i in range(1,len(mw)+1):
            if sum(mw[0:i]) < 0:
                raise ValueError("A Motzkin word is not allowed to go beneath the starting level!")
            if ((mw[i-1]) not in [1, -1, 0]):
                raise ValueError("In a Motzkin word only steps in {1,-1,0} are allowed!")
        return M(mw)

        asm = matrix(asm)
        if not asm.is_square():
            raise ValueError("The alternating sign matrices must be square")
        P = AlternatingSignMatrices(asm.nrows())
        if asm not in P:
            raise ValueError("Invalid alternating sign matrix")
        return P(asm)

    def __init__(self, parent, mw):
        """
        Initialize ``self``.
        """
        self._word=mw
        Element.__init__(self, parent)

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return "MotzkinWord {}".format(self._word)

    def number_of_up_steps(self):
        r"""
        Return the number of up steps in ``self``.

        EXAMPLES::

            sage: MotzkinWord([1, 0, 1, -1]).number_of_up_steps()
            2
            sage: MotzkinWord([1, -1, 1, 1, 0]).number_of_up_steps()
            3

        TESTS::

            sage: MotzkinWord([]).number_of_up_steps()
            0
        """
        return len([x for x in self._word if x == 1])

    def number_of_down_steps(self):
        r"""
        Return the number of down steps in ``self``.

        EXAMPLES::

            sage: MotzkinWord([1, 0, 1, -1]).number_of_down_steps()
            1
            sage: MotzkinWord([1, -1, 1, 0, -1]).number_of_down_steps()
            2

        TESTS::

            sage: MotzkinWord([]).number_of_down_steps()
            0
        """
        return len([x for x in self._word if x == -1])

    def number_of_horizontal_steps(self):
        r"""
        Return the number of horizontal steps in ``self``.

        EXAMPLES::

            sage: MotzkinWord([1, 0, 1, -1, 0]).number_of_horizontal_steps()
            2
            sage: MotzkinWord([0, 1, 0, 0, -1, 0]).number_of_horizontal_steps()
            4

        TESTS::

            sage: MotzkinWord([]).number_of_horizontal_steps()
            0
        """
        return len([x for x in self._word if x == 0])


    def ends_at_horizontal_axis(self):
        r"""
        Return ``True`` if ``self`` ends at the horizontal axis.

        A Dyck word `d` ends at the horizontal axis if `d` contains as many
        up steps as down steps.

        EXAMPLES::

            sage: DyckWord([1, 0, -1, 0, 1, -1]).is_complete()
            True
            sage: DyckWord([1, -1, 0, 1, -1, 0]).is_complete()
            True
            sage: DyckWord([1, 0, 1, -1, 0]).is_complete()
            False

        TESTS::

            sage: DyckWord([]).is_complete()
            True
        """
        return self.number_of_up_steps() == self.number_of_down_steps()

    def height(self):
        r"""
        Return the height of ``self``.

        We view the Motzkin word as a Motzkin path from `(0, 0)` to
        `(n, k1-k2)` in the first quadrant by letting ``1``'s represent
        steps in the direction `(1, 1)`, ``-1``'s represent steps in the
        direction `(1, -1)`,  and ``0``'s represent steps in the direction
        `(1, 0)`.

        The height is the maximum `y`-coordinate reached.

        EXAMPLES::

            sage: DyckWord([]).height()
            0
            sage: DyckWord([1,0]).height()
            1
            sage: DyckWord([1, 1, -1, 0, 1]).height()
            2
            sage: DyckWord([1, 1, 0, 1, -1]).height()
            3
            sage: DyckWord([1, 1, 0, -1, 0]).height()
            2
            sage: DyckWord([0, 0]).height()
            0
            sage: DyckWord([1, 1, 0, 0, 1, -1, -1, 0, 0, 0]).height()
            3
        """
        height = 0
        for i in range(0,len(self._word)):
            if sum(self._word[0:i+1]) > height:
                height=sum(self._word[0:i+1])
        return height 

    def has_horizontal_step_on_top(self):
        r"""
        Returns True if the word has a horizontal step on the highest
        height of ``self``.

        We view the Motzkin word as a Motzkin path from `(0, 0)` to
        `(n, k1-k2)` in the first quadrant by letting ``1``'s represent
        steps in the direction `(1, 1)`, ``-1``'s represent steps in the
        direction `(1, -1)`,  and ``0``'s represent steps in the direction
        `(1, 0)`.

        The height is the maximum `y`-coordinate reached.

        EXAMPLES::

            sage: DyckWord([]).height()
            False
            sage: DyckWord([1,0]).height()
            True
            sage: DyckWord([1, 1, -1, 0, 1]).height()
            False
            sage: DyckWord([1, 1, 0, 1, -1]).height()
            False
            sage: DyckWord([1, 1, 0, -1, 0]).height()
            True
            sage: DyckWord([0, 0]).height()
            True
            sage: DyckWord([1, 1, 0, 0, 1, -1, -1, 0, 0, 0]).height()
            False
        """
        horizontaltop = False
        height = 0
        for i in range(0,len(self._word)):
            if sum(self._word[0:i+1]) > height:
                height=sum(self._word[0:i+1])
                horizontaltop=False
            if (i<len(self._word)) and (sum(self._word[0:i])==height and sum(self._word[0:(i+1)])==height):
                horizontaltop=True
        return horizontaltop 


    def heights(self):
        r"""
        Return the heights of ``self``.

        We view the Motzkin word as a Motzkin path from `(0, 0)` to
        `(n, k1-k2)` in the first quadrant by letting ``1``'s represent
        steps in the direction `(1, 1)`, ``-1``'s represent steps in the
        direction `(1, -1)`,  and ``0``'s represent steps in the direction
        `(1, 0)`.

        The heights is the sequence of the `y`-coordinates of all
        `n+1` lattice points along the path.

        EXAMPLES::

            sage: DyckWord([]).heights()
            (0)
            sage: DyckWord([1,0]).heights()
            (0, 1, 1)
            sage: DyckWord([1, 1, -1, -1]).heights()
            (0, 1, 2, 1, 0)
            sage: DyckWord([0, 1, 1, -1, 1, 0]).heights()
            (0, 0, 1, 2, 1, 2, 2)
            sage: DyckWord([1, 1, -1, 0, 1, 0]).heights()
            (0, 1, 2, 1, 1, 2, 2)
            sage: DyckWord([1, 0, 1, 0]).heights()
            (0, 1, 1, 2, 2)
            sage: DyckWord([1, 1, 0, -1, 1, -1, 1, 0, -1, 0]).heights()
            (0, 1, 2, 2, 1, 2, 1, 2, 2, 1, 1)
        """
        heights=[0]*(len(self._word)+1)
        for i in range(0,len(self._word)):
            heights[i + 1] = sum(self._word[0:i+1])
        return tuple(heights)

    def pretty_print(self):
        r"""
        A path representation of the Motzkin word consisting of steps
        ``/``, ``\`` and ``_``.

        EXAMPLES::

            sage: print(MotzkinWord([1, -1, 1, -1, 0, 0]).pretty_print())
            /\/\__
            sage: print(MotzkinWord([1, 1, 0, 0, -1,-1]).pretty_print())
              __
             /  \
            /    \
            sage: print(MotzkinWord([1,1,0,-1,1,1,-1,0,-1,0,-1]).pretty_print())
              _  /\_
             / \/   \_
            /         \
        """
        space = ' '
        up    = '/'
        down  = '\\'
        hztl  = '_'

        if self.has_horizontal_step_on_top():
            res = [([space]*len(self._word)) for _ in range(self.height()+1)]
        else:
            res = [([space]*len(self._word)) for _ in range(self.height())]
        h = 1
        i = 0
        for p in self._word:
            if p == 1:
                res[-h][i] = up
                h += 1
            elif p == -1:
                h -= 1
                res[-h][i] = down
            else:
                res[-h][i] = hztl
            i=i+1
        return "\n".join("".join(l) for l in res)



class MotzkinWords(UniqueRepresentation, Parent):
    r"""
    Motzkin words.

    A Motzkin word is a sequence of up down and horizontal steps such that
    every down stap has a corresponding up step preciding it.
    A Motzkin word of length `n` is a list with `k1` entries 1, `k2` entries
    -1 and `k3` entries 0 such that `k1 + k2 + k3 = n` and the first `i`
    entries als ways have at least as many 1s among them as -1s.

    EXAMPLES::

    """

    def __init__(self,n):
        """
        Intialize ``self``.

        EXAMPLES::

            sage: TestSuite(MotzkinWords).run()
        """
        self._n = n
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = MotzkinWords(4); M
            Motzkinwords of length 4
        """
        return "Motzkinwords of length {}".format(self._n)


#    def _element_constructor_(self, mword):
#        """
#        Construct an element of ``self``.

#        EXAMPLES::

#            sage: D = MotzkinWords()
#            sage: elt = D([1, 1, 0, 1, 0, 0]); elt
#            [1, 1, 0, 1, 0, 0]
#            sage: elt.parent() is D
#            True
#        """
#        if isinstance(mword, MotzkinWord) and mword.parent() is self:
#            return mword
#        return self.element_class(self, list(mword))

    def __iter__(self):
        """
        Iterate over ``self``.
        """
        n = self._n
        for k1 in range(0, n+1):
            for k2 in range(0, min([n-k1,k1])+1):
                if n == 0:
                    yield self.element_class(self, [])
                elif (k1==n):
                    yield self.element_class(self, [1]*k1)
                elif n==1 and (k1 == 0 and k2==0):
                    yield self.element_class(self, [0])
                else:
                    for w in MotzkinWordBacktracker(k1, k2, n-k2-k1):
                        yield self.element_class(self, w)

    Element=MotzkinWord

class MotzkinWordBacktracker(GenericBacktracker):
    r"""
    This class is an iterator for all Motzkin words
    with `k1` up steps, `k2` down steps and `k3` horizontal steps using
    the backtracker class. It is used by the :class:`MotzkinWords` class.

    This is not really meant to be called directly.
    """

    def __init__(self, k1, k2, k3):
        r"""
        init for Generic Backtracker
        """
        self.k1 = ZZ(k1)
        self.k2 = ZZ(k2)
        self.k3 = ZZ(k3)
        self.n = k1 + k2 + k3
        GenericBacktracker.__init__(self, [], (0, 0, 0))

    def _rec(self, path, state):
        r"""
        _rec for Generic Backtracker
        """
        l1, l2, l3 = state
        if (l1+l2+l3) < (self.n - 1):
            if l1 < self.k1:
                yield path + [1], (l1+1,l2,l3), False
            if (l2 < self.k2) and (l2<l1):
                yield path + [-1], (l1,l2+1,l3), False
            if l3 < self.k3:
                yield path + [0], (l1,l2,l3+1), False
        else:
            if l1 < self.k1:
                yield path + [1], None, True
            elif l2 < self.k2:
                yield path + [-1], None, True
            else:
                yield path + [0], None, True
