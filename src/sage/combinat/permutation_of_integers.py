r"""
Permutations of integers which move finitely many elements.

This is implemented as a permutation of an interval of integers::

    sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
    sage: p = PermutationOfIntegers(-2, 4, [-1, 0, 3, 2, -2, 1, 4]); p
    Permutation on [-2,4]: [-1, 0, 3, 2, -2, 1, 4]

"""
from copy import deepcopy
from sage.structure.sage_object import SageObject

class PermutationOfIntegers(SageObject):
    def __init__(self, low, high, lis):
        r"""
        Initialize ``self``.
        """
        self._low = low
        self._high = high
        self._values = deepcopy(lis)

    def _repr_(self):
        return "Permutation on [{},{}]: {}".format(self.low(), self.high(), self.values())

    def low(self):
        r"""
        The minimum of the domain of ``self``.
        
        EXAMPLES::

            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
            sage: p = PermutationOfIntegers(-2, 4, [-1, 0, 3, 2, -2, 1, 4])
            sage: p.low()
            -2
        """
        return self._low

    def high(self):
        r"""
        The maximum of the domain of ``self``.
        
        EXAMPLES::

            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
            sage: p = PermutationOfIntegers(-2, 4, [-1, 0, 3, 2, -2, 1, 4])
            sage: p.high()
            4
        """
        return self._high

    def values(self):
        r"""
        The list of values of ``self``.
        
        EXAMPLES::

            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
            sage: p = PermutationOfIntegers(-2, 4, [-1, 0, 3, 2, -2, 1, 4])
            sage: p.values()
            [-1, 0, 3, 2, -2, 1, 4]
        """
        return self._values

    def value(self, i):
        r"""
        The value of ``self`` at `i`.
        
        EXAMPLES::

            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
            sage: p = PermutationOfIntegers(-2, 4, [-1, 0, 3, 2, -2, 1, 4])
            sage: [[i, p.value(i)] for i in range (-4, 6)]
            [[-4, -4], [-3, -3], [-2, -1], [-1, 0], [0, 3], [1, 2], [2, -2], [3, 1], [4, 4], [5, 5]]
        """
        if i < self._low or i > self._high:
            return i
        return self._values[i-self._low]

    def func(self):
        r"""
        The function from integers to integers defined by ``self``.
        
        EXAMPLES::

            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
            sage: p = PermutationOfIntegers(-2, 4, [-1, 0, 3, 2, -2, 1, 4])
            sage: f = p.func()
            sage: [f(-5), f(3), f(-1)]
            [-5, 1, 0]
        """
        return lambda i: self.value(i)

    def left_action_product(self, w):
        r"""
        Compute ``w`` followed by ``self``.

        EXAMPLES::

            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
            sage: p = PermutationOfIntegers(-2, 4, [-1, 0, 3, 2, -2, 1, 4])
            sage: p.left_action_product(p)
            Permutation on [-2,4]: [0, 3, 1, -2, -1, 2, 4]
        """
        low = min(self.low(), w.low())
        high = max(self.high(), w.high())
        value_list = [self.value(w.value(i)) for i in range(low, high+1)]
        return PermutationOfIntegers(low, high, value_list)

    def nonpos_to_pos_set(w):
        r"""
        The set of nonpositive integers that `w` sends to positive.

        EXAMPLES::

            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
            sage: p = PermutationOfIntegers(-1,3,[1,3,-1,0,2])
            sage: p.nonpos_to_pos_set()
            [-1, 0]
        """
        if w.low() >= 1 or w.high() <= 0:
            return []
        return [i for i in range(w.low(), 1) if w.value(i)>0]

    def pos_to_nonpos_set(w):
        r"""
        The set of positive integers that `w` sends to nonpositive.

        EXAMPLES::

            sage: from sage.combinat.permutation_of_integers import PermutationOfIntegers
            sage: p = PermutationOfIntegers(-1,3,[1,3,-1,0,2])
            sage: p.pos_to_nonpos_set()
            [1, 2]
        """
        if w.low() >= 1 or w.high() <= 0:
            return []
        return [i for i in range(1,w.high()+1) if w.value(i)<=0]

