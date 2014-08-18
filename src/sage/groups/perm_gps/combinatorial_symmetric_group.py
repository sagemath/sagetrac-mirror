from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.symmetrica.all import charvalue
from sage.combinat.partition import Partitions, Partition
from sage.combinat.permutation import Permutation
from itertools import imap

class SymmetricGroupCombinatorial(SymmetricGroup):
    """
    The combinatorial version of the symmetric group.
    """
    def conjugacy_classes(self):
        return CombinatorialSymmetricGroupConjugacyClasses(self.degree())

    def irreducible_character(self, la):
        if sum(la)<>self._deg:
            raise ValueError("Input should be a partition of %s"%(self._deg))
        return lambda mu: charvalue(la,mu)

class CombinatorialSymmetricGroupConjugacyClass(ElementWrapper):
    """
    The combinatorial version of a conjugacy class in a symmetric group.
    """
    def __init__(self, parent, la):
        la = Partition(la)
        ElementWrapper.__init__(self, parent, la)
        
    def _repr_(self):
        """
        TESTS::

            sage: CombinatorialSymmetricGroupConjugacyClass(CombinatorialSymmetricGroupConjugacyClasses(3), [2,1])
            Permutations with cycle type [2,1]
        """
        return "Permutations with cycle type {}".format(self.value)

    def size(self):
        """
        Return the size of self.

        TESTS::

            sage: CombinatorialSymmetricGroupConjugacyClasses(3)([2,1]).size()
            3
        """
        return self.value.size()
    
    def partition(self):
        """
        Return the cycle type of elements of self.

        TESTS::

            sage: CombinatorialSymmetricGroupConjugacyClasses(3)([2,1]).partition()
            [2, 1]
        """
        return self.value
    
    def representative(self):
        """
        Return an element of a symmetric group which lies in self.

        TESTS::

            sage: CombinatorialSymmetricGroupConjugacyClasses(3)([2,2,1,1]).representative()
            (1,2)(3,4)
        """
        L = list()
        S = 0
        for i in range(len(self.value)):
            T = S + self.value[i]
            L.append(tuple(range(S+1, T+1)))
            S = T
        return SymmetricGroup(self.size())(L)

class CombinatorialSymmetricGroupConjugacyClasses(Parent, UniqueRepresentation):
    def __init__(self, n):
        Parent.__init__(self, category = FiniteEnumeratedSets())
        self._n = n
        
    def _repr_(self):
        return "Conjugacy classes of the symmetric group of order {}".format(self._n)
    
    def __iter__(self):
        for la in Partitions(self._n):
            yield self(la)
            
    Element = CombinatorialSymmetricGroupConjugacyClass
