from sage.categories.combinatorial_structures import CombinatorialStructures
from sage.categories.category import Category
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.species import Species
from sage.combinat.partition import Partitions
from sage.misc.lazy_import import LazyImport
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.positive_integers import PositiveIntegers
from sage.sets.set import Set
from sage.sets.disjoint_set import DisjointSet
from sage.structure.element_wrapper import ElementWrapper
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.structures import Structures as Structs
Permutation = LazyImport('sage.combinat.permutation', 'Permutation')


def partition_to_permutation(pi):
    it = iter(PositiveIntegers())
    return Permutation([tuple([it.next() for _ in range(part)])
                        for part in pi])

class ClassOfIsoTypes(Structs):
    """
    An abstract class for all parents modelling the set of isormophism types of a species
    """
    def __init__(self, species):
        Structs.__init__(self)
        self._species_ = species

    Element = ElementWrapper

    def _repr_(self):
        return "Isomorphism types of %s"%repr(self._species_)

    def generating_series(self):
        return self._species_.isomorphism_type_generating_series()

    class GradedComponent(Structs.GradedComponent):
        # In this default implementation, the iso types are computed
        # by computing the orbits in the list of structures for the
        # action of the symmetric group. This is done lazily by
        # compute_equivalence_class the first time __iter__ is called,
        # and using DisjointSet.
        #
        # Whenever possible, a better implementation should be provided

        def __init__(self, ambient, n):
            Structs.GradedComponent.__init__(self, ambient, n)
            # This could invoke directly self.species().structures()
            FU = self.ambient()._species_.structures(Set(range(1, self.grading()+1)))
            self._isotypes_ = DisjointSet(FU)
            self._is_computed_ = False

        def compute_equivalence_class(self):
            F = self.ambient()._species_
            FU = F.structures(Set(range(1, self.grading()+1)))
            for s in FU:
                for pi in Partitions(self.grading()):
                    sigma = partition_to_permutation(pi)
                    t = F.transport(sigma)(s)
                    self._isotypes_.union(s, t)

        def __iter__(self):
            if not self._is_computed_:
                self.compute_equivalence_class()
            for ts in self._isotypes_.root_to_elements_dict().values():
                yield ElementWrapper(self.ambient(), ts)

        @cached_method
        def cardinality(self):
            F = self.ambient()._species_
            return F.isomorphism_type_generating_series().coefficient(self.grading())


class SpeciesDesign(UniqueRepresentation, Parent):
    """
    An abstract class for a species.

    Would be naturally called Species, but we currently use that name
    for the category (Species is spelled the same way in singular or plural) 
    """

    def __init__(self):
        Parent.__init__(self, category=Category.join([Species(), CombinatorialStructures()]))
        self._structures = {}

    # At some point there was a (good?) reason for not using a @cached_method?
    # We may want to have a weak cache or ...
    def structures(self, U):
        if U in self._structures:
            return self._structures[U]
        FU = self.Structures(self, U)
        self._structures[U] = FU
        FU._ambient_ = self     # could this be handled by the constructor of Structures?
        return FU

    def graded_component(self, k):
        return self.structures(Set(range(1, k+1)))

    def _an_element_(self):
        return self.first()

    def _element_constructor_(self, *args, **options):
        """
        Redefinition of that method to be coherent with the *_classcall_* of
        *Structure*.
        """
        return self.element_class(parent=self, *args, **options)

    class Structures(UniqueRepresentation, Parent):
        """
        An abstract class for all parents modelling the set of structures of a species
        """
        # TODO (to be confirmed): rename ambient into species
        def __init__(self, ambient, U):
            """

            :param U: a finite set


            """
            Parent.__init__(self, category=FiniteEnumeratedSets())
            self._finite_set_ = U
            self._ambient_ = ambient

        def ambient(self):
            return self._ambient_

        def finite_set(self):
            """
            Return the set of labels

            Check Bergeron what the right name for this method could be. label_set?
            """
            return self._finite_set_

        def cardinality(self):
            # self.species instead of self.ambient
            return self.ambient().exponential_generating_series().coefficient(self.grading())

        def grading(self):
            return self._finite_set_.cardinality()

        def _repr_(self):
            # self.species instead of self.ambient
            return repr(self.ambient()) + "-structures on " + repr(self.finite_set())

        @lazy_attribute
        def _element_constructor_(self):
            return self.ambient()._element_constructor_

    IsomorphismTypes = ClassOfIsoTypes
