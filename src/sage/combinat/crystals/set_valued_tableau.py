from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.graphs.all import DiGraph
from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType

from sage.misc.lazy_attribute import lazy_attribute

class HeckeMonoidFactorizationCrystal(UniqueRepresentation, Parent):
    r"""
    The crystal on Hecke monoid factorizations.

    INPUT:

    - ``w`` -- an element in a symmetric group

    - ``n`` -- the number of factors in the factorization

    EXAMPLES::
        
        sage: W = WeylGroup(['A',3,1], prefix='s')
        sage: w = W.from_reduced_word([2,3,2,1])
        sage: B = crystals.HeckeMonoidFactorizationCrystal(w,3); B
        Crystal on Hecke monoid factorizations associated to s2*s3*s2*s1    
    """
    '''
    ## YET TO BE IMPLEMENTED ##
    @staticmethod
    def __classcall_private__(cls, w, n, x = None, k = None):
        r"""
        Classcall to mend the input.

        TESTS::

            sage: A = crystals.AffineFactorization([[3,1],[1]], 4, k=3); A
            Crystal on affine factorizations of type A3 associated to s3*s2*s1
            sage: AC = crystals.AffineFactorization([Core([4,1],4),Core([1],4)], 4, k=3)
            sage: AC is A
            True
        """ 
        pass
        '''

    def __init__(self, w, n):
        r"""
        EXAMPLES::
            
            sage: W = WeylGroup(['A',3], prefix='s')
            sage: w = W.from_reduced_word([2,3,2,1])
            sage: B = crystals.HeckeMonoidFactorizationCrystal(w,4)
            sage: B.w
            s2*s3*s2*s1
            sage: B.k
            3
            sage: B.n
            4

        TESTS::

            sage: W = WeylGroup(['A',3], prefix='s')
            sage: w = W.from_reduced_word([2,3,2,1])
            sage: B = crystals.HeckeMonoidFactorizationCrystal(w,4)
            sage: TestSuite(B).run()
        """
        Parent.__init__(self, category = ClassicalCrystals())
        # n is the number of parenthesis pairs
        self.n = n
        # k is the highest letter that appears in the factorization
        self.k = w.parent().n-1
        # w is a given word in the symmetric group
        self.w = w
        cartan_type = CartanType(['A',n-1])
        self._cartan_type = cartan_type
        ## YET TO BE IMPLEMENTED ##
        ## self.module_generators (this enumerates all highest weight vectors!)

    def __repr__(self):
        """
        EXAMPLES::

            sage: 
        """
        return "Crystal of 0-Hecke monoid factorizations of type A_{} associated to {}".format(self.n-1, self.w)

    # temporary workaround while an_element is overriden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_

    class Element(ElementWrapper):
            

        def e(self, i):
            r"""
            Returns the action of `e_i` on ``self``.

            EXAMPLES::

                sage: 
            """

        def f(self, i):
            r"""
            Returns the action of `f_i` on ``self``.

            EXAMPLES::

                sage: 
            """
            if i not in self.index_set():
                raise ValueError("i must be in the index set")
            
