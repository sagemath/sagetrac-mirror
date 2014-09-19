from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.categories.species import Species as SpeciesCat
from sage.combinat.structures import Structures
from sage.sets.set import Set
from sage.structure.parent import Parent

class Species(Structures):
    
    def __init__(self, category=SpeciesCat()):
        Parent.__init__(self, category=category)

    def graded_component(self, k):
        return self.structures(range(1, k+1))

    def structures(self, U):
        return self.Structures(self, Set(U))

    def transport(self, sigma):
        #assert(sigma in Bijections())
        return SetMorphism(Hom(
                self.structures(sigma.domain()),
                self.structures(sigma.codomain())
            ), lambda obj: obj.transport(sigma)
        )

    class Structures(Structures.GradedComponent):

        def __init__(self, ambient, U, category=SpeciesCat.Structures()):
            """
            TESTS::

                sage: from sage.categories.examples.combinatorial_structures_compositions import Compositions
                sage: Compositions(4)
                Compositions of integers of degree 4
                sage: Compositions().graded_component(4)
                Compositions of integers of degree 4
                sage: Compositions(3).list()
                [[3], [1, 2], [2, 1], [1, 1, 1]]
            """
            Parent.__init__(self, category=category)
            self._ambient_ = ambient
            self._underlying_set_ = U

        def underlying_set(self):
            return self._underlying_set_

        def grading(self):
            return len(self.underlying_set())

        def __contains__(self, obj):
            """
            TESTS::

            """
            if not isinstance(obj, self.ambient().element_class):
                try:
                    obj = self._element_constructor_(obj)
                except: return False
            return obj in self.ambient() and \
                   obj.underlying_set() == self.underlying_set()
