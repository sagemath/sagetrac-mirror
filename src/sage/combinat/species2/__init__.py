from sage.combinat.permutation import Permutation
from sage.sets.positive_integers import PositiveIntegers


def partition_to_permutation(pi):
    it = iter(PositiveIntegers())
    return Permutation([tuple([it.next() for _ in range(part)]) for part in pi])
