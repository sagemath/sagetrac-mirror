from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.graphs.all import DiGraph
from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType

from sage.misc.lazy_attribute import lazy_attribute

class SemistandardSetValuedTableaux(UniqueRepresentation, Parent):
    r"""
    The crystal on semistandard set-valued tableaux.

    INPUT:

    - ```` -- an element in a symmetric group

    - ``n`` -- the number of factors in the factorization

    EXAMPLES::
        
        sage: W = WeylGroup(['A',3,1], prefix='s')
        sage: w = W.from_reduced_word([2,3,2,1])
        sage: B = crystals.HeckeMonoidFactorizationCrystal(w,3); B
        Crystal on Hecke monoid factorizations associated to s2*s3*s2*s1    
    """
   
    @staticmethod
    def __classcall_private__(cls, shape, n):
        r"""
        Classcall to mend the input.

        TESTS::

            sage: B = crystals.SemistandardSetValuedTableaux([2,1],3); B
            Crystal of set-valued tableaux of type A_2 of shape [2,1]
        """ 
        
        

    def __init__(self, sh, n):
        r"""
        Initialize crystal of semistandard set-valued tableaux of a fixed shape and given maximum entry. 
        
        EXAMPLES::
            
            sage: B = crystals.SemistandardSetValuedTableaux([2,1],3)
            sage: B.n
            3
            sage: B.shape
            [2,1]

        TESTS::

            sage: W = WeylGroup(['A',3], prefix='s')
            sage: w = W.from_reduced_word([2,3,2,1])
            sage: B = crystals.HeckeMonoidFactorizationCrystal(w,4)
            sage: TestSuite(B).run()
        """
        Parent.__init__(self, category = ClassicalCrystals())
        # n is the max entry
        self.n = n
        cartan_type = CartanType(['A',n-1])
        self._cartan_type = cartan_type
        # (this enumerates all highest weight vectors!)
        self.module_generators = [] 

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: B = crystals.SemistandardSetValuedTableaux([2,1],3); B
            Crystal of set-valued tableaux of type A_2 of shape [2,1]
        """
        return "Crystal of set-valued tableaux of type A_{} of shape {}".format(self.n-1, self.shape)

    # temporary workaround while an_element is overriden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_


    class Element(ElementWrapper):
    
        def __init__(self,parent,tab):
            r"""
            Initialize self as a crystal element.

            EXAMPLES::

            sage: SVT = crystals.SemistandardSetValuedTableau([2,1],3)
            sage: T = SVT([ [[1,2],[2,3]],[[3]] ]); T
            [[[1, 2], [2,3]], [[3]] ]
            """
            self.n = parent.n
            self.value = SemistandardSetValuedTableau(tab)
            if hf.excess != excess:
                raise ValueError("The HeckeFactorization word must have correct excess")
            ElementWrapper.__init__(self, parent, self.value)

        def _get_signs(self, i):
            """
            Auxiliary function for `e_i` and `f_i` methods.
    
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
            col = _bracket(self,i,right=False)
            if col == -1:
                return None
            tab = [[[entry] for entry in r] for r in self]
            t = self.conjugate()
            column = t[col]
            row = min([ j for j in range(len(column)) if i+1 in column[j] ])
            # checks that there is a cell to the left and that the cell contains i and i+1
            if col>=1 and all(x in tab[row][col-1] for x in [i,i+1]):
                tab[row][col-1].remove(i+1)
                tab[row][col] = sorted(tab[row][col]+[i])
            else:
                tab[row][col].remove(i+1)
                tab[row][col] = sorted(tab[row][col]+[i])
            return self.parent()(tab)
            
        def f(self, i):
            r"""
            Returns the action of `f_i` on ``self``.

            EXAMPLES::

                sage: 
            """
            if i not in self.index_set():
                raise ValueError("i must be in the index set")
            col = _bracket(self,i,right=True)
            if col == -1:
                return None
            tab = [[[entry] for entry in r] for r in self]
            t = self.conjugate()
            column = t[col]
            row = min([ j for j in range(len(column)) if i in column[j] ])
            # checks that there is a cell to the right and that the cell contains i and i+1
            if col<len(self[row])-1 and all(x in tab[row][col+1] for x in [i,i+1]):
                tab[row][col+1].remove(i)
                tab[row][col] = sorted(tab[row][col]+[i+1])
            else:
                tab[row][col].remove(i)
                tab[row][col] = sorted(tab[row][col]+[i+1])
            return self.parent()(tab)

