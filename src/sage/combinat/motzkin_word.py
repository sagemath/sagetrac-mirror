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
from sage.structure.list_clone import ClonableArray
from .backtrack import GenericBacktracker
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

#    def check(self):
#        """
#        TODO CHANGE THIS DOCUMENTTION
#        Check if ``self`` is a valid 6 vertex configuration.

#        EXAMPLES::
#
#            sage: M = SixVertexModel(3, boundary_conditions='ice')
#            sage: M[0].check()
#        """
#        if self not in self.parent():
#            raise ValueError("invalid configuration")



#    __metaclass__ = InheritComparisonClasscallMetaclass

#    @staticmethod
#    def __classcall_private__(cls, mword):
#        """
#        Create a MotzkinWord.

#        EXAMPLES::

#            sage: MotzkinWord([1,1,0,-1,1,-1,0])
#            [1,1,0,-1,1,-1,0]
#        """
#        mword = list(mword)
#        for i in range(1,len(mword)+1):
#            if sum(mword[0:i]) < 0:
#                raise ValueError("A Motzkin word is not allowed to go beneath the starting level!")
#            if ((mword[i-1]) not in [1, -1, 0]):
#                raise ValueError("In a Motzkin word only steps in {1,-1,0} are allowed!")
#        return mword

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
        return len([x for x in self if x == -1])

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
        return len([x for x in self if x == 0])


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

        .. SEEALSO:: :meth:`heights`

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
        for i in range(1,self._n+1):
            if sum(self[1:i]) > height:
                height=sum(self[1:i])
        return height 
#        height = 0
#        height_max = 0
#        for letter in self:
#            if letter == 1:
#                height += 1
#                height_max = max(height, height_max)
#            elif letter == close_symbol:
#                height -= 1
#        return height_max

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
        for i in range(1,self._n + 1):
            if sum(self[1:i]) > height:
                height=sum(self[1:i])
                horizontaltop=False
            if (i<self._n) and (sum(self[1:i])==height and sum(self[1:(i+1)])==height):
                horizontaltop=True
        return height 


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
        heights=[0]*(len(self)+1)
        for i in range(1,self._n+1):
            heights[i + 1] = sum(self[1:i])
#        height = 0
#        heights = [0] * (len(self) + 1)
#        for i, letter in enumerate(self):
#            if letter == open_symbol:
#                height += 1
#            elif letter == close_symbol:
#                height -= 1
#            heights[i + 1] = height
        return tuple(heights)


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

#    __metaclass__ = InheritComparisonClasscallMetaclass

#    @staticmethod
#    def __classcall_private__(cls, n):
#        """
#        Choose the correct parent based upon input.
#
#        EXAMPLES::
#
#            sage: DW1 = DyckWords(3,3)
#            sage: DW2 = DyckWords(3)
#            sage: DW1 is DW2
#            True
#        """
#        return MotzkinWords()

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
        return "Motzkinwords of length %s" % self._n


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

        EXAMPLES::

            sage: it = MotzkinWords.__iter__()
            sage: [next(it) for x in range(10)]

        """
        n = 0
        yield self.element_class(self, [])
        while True:
            for k1 in range(1, n+1):
                for k2 in range(1, min([n+1-k1,k1])):
                    for x in MotzkinWords_size(k1, k2, n-k1-k2):
                        yield self.element_class(self, list(x))
            n += 1

    Element=MotzkinWord

class MotzkinWordBacktracker(GenericBacktracker):
    r"""
    This class is an iterator for all Dyck words
    with `k1` up steps, `k2` down steps and `k3` horizontal steps using
    the backtracker class. It is used by the :class:`MotzkinWords_size` class.

    This is not really meant to be called directly.
    """

    def __init__(self, k1, k2, k3):
        r"""
        TESTS::

        """
        GenericBacktracker.__init__(self, [], (0, 0))
        # note that the comments in this class think of our objects as
        # Dyck paths, not words; having k1 opening parens and k2 closing
        # parens corresponds to paths of length k1 + k2 ending at height
        # k1 - k2.
        k1 = ZZ(k1)
        k2 = ZZ(k2)
        k3 = ZZ(k3)
        self.n = k1 + k2 + k3
        self.endht = k1 - k2

    def _rec(self, path, state):
        r"""
        TESTS::
        !!!SOMETHING LIKE THAT BUT DIFFERENT (TO-DO)
            sage: from sage.combinat.dyck_word import DyckWordBacktracker
            sage: dwb = DyckWordBacktracker(3, 3)
            sage: list(dwb._rec([1,1,0],(3, 2)))
            [([1, 1, 0, 0], (4, 1), False), ([1, 1, 0, 1], (4, 3), False)]
            sage: list(dwb._rec([1,1,0,0],(4, 0)))
            [([1, 1, 0, 0, 1], (5, 1), False)]
            sage: list(DyckWordBacktracker(4, 4)._rec([1,1,1,1],(4, 4)))
            [([1, 1, 1, 1, 0], (5, 3), False)]
        """
        len, ht = state
        if len < self.n - 1:
            # if length is less than n-1, new path won't have length n, so
            # don't yield it, and keep building paths

            # if the path isn't too low and is not touching the x-axis, we can
            # yield a path with a downstep at the end
            if ht > (self.endht - (self.n - len)) and ht > 0:
                yield path + [-1], (len + 1, ht - 1), False

            # if the path isn't too high, it can also take an upstep
            if ht < (self.endht + (self.n - len)):
                yield path + [1], (len + 1, ht + 1), False

            # there can always be a horizontal step
            yield path + [0], (len + 1, ht), False

        else:
            # length is n - 1, so add a single step (up, down or horizontal,
            # according to current height and endht), don't try to
            # construct more paths, and yield the path
            if ht < self.endht:
                yield path + [1], None, True
            elif ht > self.endht:
                yield path + [-1], None, True
            else:
                yield path + [0], None, True


class MotzkinWords_size(MotzkinWords):
    """
    Motzkin words with `k1` up steps, `k2` down steps  and `k_3` horizontal
    steps
    """
    def __init__(self, k1, k2, k3):
        r"""
        TESTS:

        """
        self.k1 = ZZ(k1)
        self.k2 = ZZ(k2)
        self.k3 = ZZ(k3)
        MotzkinWords.__init__(self, k1+k2+k3)

    def __repr__(self):
        r"""
        TESTS::

            sage: MotzkinWords(4,2,3)
            Motzkin words with 4 upsteps, 2 downsteps and 3 horziontalsteps
        """
        return "Motzkin words with %s up steps, %s down steps and %s horziontal steps" % (self.k1, self.k2, self.k3)


    def __iter__(self):
        r"""
        Return an iterator for Motzkin words with `k1` up steps, `k2` down
        steps  and `k_3` horizontal steps.

        EXAMPLES::

        !!!SOMETHING LIKE THAT BUT DIFFERENT (TO-DO)
            sage: list(DyckWords(0))
            [[]]
            sage: list(DyckWords(1))
            [[1, 0]]
            sage: list(DyckWords(2))
            [[1, 0, 1, 0], [1, 1, 0, 0]]
            sage: len(DyckWords(5))
            42
            sage: list(DyckWords(3,2))
            [[1, 0, 1, 0, 1],
             [1, 0, 1, 1, 0],
             [1, 1, 0, 0, 1],
             [1, 1, 0, 1, 0],
             [1, 1, 1, 0, 0]]
        """
        if self.k1 == 0:
            yield self.element_class(self, [])
        elif (self.k2 == 0 and self.k3==0):
            yield self.element_class(self, [1]*self.k1)
        elif (self.k1 == 0 and self.k2==0):
            yield self.element_class(self, [0]*self.k1)
        else:
            for w in MotzkinWordBacktracker(self.k1, self.k2, self.k3):
                yield self.element_class(self, w)

    Element=MotzkinWord

#    def __contains__(self, x):
#        r"""
#        IRGENDSOWAS TODO!!
#        EXAMPLES::

#            sage: [1, 0, 0, 1] in DyckWords(2,2)
#            False
#            sage: [1, 0, 1, 0] in DyckWords(2,2)
#            True
#            sage: [1, 0, 1, 0, 1] in DyckWords(3,2)
#            True
#            sage: [1, 0, 1, 1, 0] in DyckWords(3,2)
#            True
#            sage: [1, 0, 1, 1] in DyckWords(3,1)
#            True
#        """
#        return is_a(x, self.k1, self.k2, self.k3)
