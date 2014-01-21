"""
Finiteness and Infiniteness Axioms
"""

from axiom import Axiom


class Finite(Axiom):
    
    class ParentMethods:
        
        def is_finite(self):
            """
            Returns ``True`` since ``self`` is finite.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.is_finite()
                True
            """
            return True
    

class Infinite(Axiom): 

    class ParentMethods:

        def is_finite(self):
            """
            Returns ``False`` since self is not finite.

            EXAMPLES::

                sage: C = InfiniteEnumeratedSets().example()
                sage: C.is_finite()
                False

            TESTS::

                sage: C.is_finite.im_func is sage.categories.sets_cat.Sets.Infinite.ParentMethods.is_finite.im_func
                True
            """
            return False

        def cardinality(self):
            """
            Counts the elements of the enumerated set.

            EXAMPLES::

                sage: NN = InfiniteEnumeratedSets().example()
                sage: NN.cardinality()
                +Infinity
            """
            from sage.rings.infinity import infinity
            return infinity


class FiniteDimensional(Axiom):
    pass
