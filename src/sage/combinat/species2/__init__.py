from sage.structure.parent import Parent
from sage.categories.species import Species
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.partition import Partitions
from sage.combinat.permutation import Permutation
from sage.sets.positive_integers import PositiveIntegers
from sage.sets.set import Set
from sage.sets.disjoint_set import DisjointSet
from sage.structure.element_wrapper import ElementWrapper
from sage.misc.cachefunc import cached_method


def partition_to_permutation(pi):
    it = iter(PositiveIntegers())
    return Permutation([tuple([it.next() for _ in range(i)])
                        for i in pi])

class SpeciesStructures(Parent):

    def __init__(self, U):
        Parent.__init__(self, category=FiniteEnumeratedSets())
        assert(U in FiniteEnumeratedSets())

    def ambient(self):
        return self._ambient_

class DefaultSpeciesIsoTypes(Parent):

    def __init__(self, n):
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._FU_ = self.ambient().structures(Set(range(1,n+1)))
        self._n_ = n
        self._isotypes_ = DisjointSet(self._FU_)
        self.compute_equivalence_class()

    def compute_equivalence_class(self):
        for s in self._FU_:
            for pi in Partitions(self._n_):
                sigma = partition_to_permutation(pi)
                self._isotypes_.union(s, self.ambient().transport(sigma)(s))

    def ambient(self):
        return self._ambient_

    def __iter__(self):
        for ts in self._isotypes_.root_to_elements_dict().values():
            yield ElementWrapper(self, ts)

    @cached_method
    def cardinality(self):
        return self._isotypes_.number_of_subsets()

class SpeciesDesign(Parent):

    def __init__(self):
        Parent.__init__(self, category=Species())

    def structures(self, U):
        FU = self.Structures(U)
        FU._ambient_ = self
        return FU

    Structures = NotImplemented

    def isomorphism_types(self, n):
        return self.SpeciesIsoTypes(n)

    SpeciesIsoTypes = DefaultSpeciesIsoTypes
