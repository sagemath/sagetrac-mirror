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
            

        def _get_signs(self, i):
            """
            auxiliary function for `e_i` and `f_i` methods.
    
            Assign each column of self a +1, -1 or 0 according to
            +1 if there is an unmatched `i+1` aka left paren '('
            -1 if there is an unmatched `i` aka right paren ')'
            0 if all the `i`s and `i+1`s are matched
            
            Return list of +1, -1, 0 with length equal to number of columns of self.
            """
            st = SemistandardSetValuedTableau(self)
            signs = []
            for col in st.conjugate():
                word = sum(col, [])
                if i in word and i+1 in word:
                    signs += [0]
                elif i in word: #i in word, i+1 not in word
                    signs += [-1]
                elif i+1 in word: # i not in word, i+1 is
                    signs += [+1]
                else: # neither i nor i+1 in word
                    signs += [0]
            return signs

        def _bracket(self, i, right=True):
            """
            auxiliary function for `e_i` and `f_i` methods.

            Return index of column in self with rightmost `i` to be changed to `i+1` 
            or leftmost `i+1` to be changed to `i`.

            If right is True (default), then return index of rightmost `i`.
            If right is False, then return index of leftmost `i+1`.

            If no `i` can be changed to `i+1` or vice versa, return -1.
            """
            x = _get_signs(self, i) # x is a list of +1, -1, 0
            if not right:
                x = [-j for j in x][::-1]
            count = 0
            index = -1
            for j in range(len(x)):
                if x[j] == -1:
                    if count == 0:
                        index = j
                    else:
                        count -= 1
                if x[j] == +1:
                    count += 1
            if right:
                return index
            else:
                return -1 if index < 0 else len(x)-1-index

        def e(self, i):
            r"""
            Returns the action of `e_i` on ``self``.

            EXAMPLES::

                sage: 
            """
            if i not in self.index_set():
                raise ValueError("i must be in the index set")
            pass

        def f(self, i):
            r"""
            Returns the action of `f_i` on ``self``.

            EXAMPLES::

                sage: 
            """
            if i not in self.index_set():
                raise ValueError("i must be in the index set")
            pass
            
