


# Will eventually need to figure out import statements
# ...making tentative list below

from six import add_metaclass

from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.rings.all import ZZ
from sage.structure.list_clone import ClonableArray
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass

from sage.rings.integer import Integer
from sage.combinat.combinatorial_map import combinatorial_map

@add_metaclass(InheritComparisonClasscallMetaclass)
class IncreasingLabeling(ClonableArray):
    @staticmethod
    def __classcall_private__(cls, poset, labels):
        if isinstance(IncreasingLabeling, cls):
            return labels
        return LinearExtensionsOfPoset(poset)(labels)
    
    def check(self):
        # A function to check that we have an increasing labeling
        # I think we want this check in IncreasingLabelings classes, also
        P = self.parent().poset()
        is_inc_label = True
        for elem in P:
            for lc in P.lower_covers(elem):
                if restrict[lc]>=restrict[elem]:
                    is_inc_label = False
        if is_inc_label is False:
            raise ValueError("%s is not an increasing labeling of %s"%(self, P))

class IncreasingLabelings(UniqueRepresentation, Parent):
#    def __classcall_private__(self, n=None, **kwargs)
    # Class call will redirect to the right parent based on input
    # Ie, we have one parent for increasing labelings with given function
    # one parent for largest max entry, and this is the 'parent' of the parents
    # that decides where things go

    def __classcall_private__(self, poset,n=None,restrict=None):
        if n is None:
            if restrict is None:
                return IncreasingLabelings_all(poset)
            else:
                return IncreasingLabelings_restrict(poset,restrict=restrict)     
        else:
            if restrict is None:
                if isinstance(n, (int,Integer)):
                    return IncreasingLabelings_n(poset,n=n)
        if is_infinite:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=FiniteEnumeratedSets())

    Element = IncreasingLabeling

#    def _element_constructor_(self, lst, check=True):
#        if isinstance(lst, IncreasingLabeling):
#            lst = list(lst)
#        if not isinstance(lst, (list, tuple)):
#            raise TypeError("input should be a list or tuple")
#        lst = [self._poset(_) for _ in lst]
#        if self._is_facade:
#            return lst
#        else:
#            return self.element_class(self, lst, check)       

#    def __contains__(self, x):    



    
class IncreasingLabelings_all(IncreasingLabelings):

    @staticmethod
    def __classcall_private__(cls,poset):
        self.poset=poset
        
    def __init__(self, is_infinite=True):
            IncreasingLabelings.__init__(self, True)
    
class IncreasingLabelings_restrict(IncreasingLabelings):
# Class of increasing labelings with function for range
# of possible values 'restrict'

    def __init__(self, poset, restrict):
        """
        Initialize ``self``.

        """
        self.restrict=restrict
        IncreasingLabelings.__init__(self, poset, False)
        
        
#    def __iter__(self):
#        for il in inc_labels(self._poset, self.restrict)
#            yield self.element_class(self, il)
        

    
class IncreasingLabelings_n(IncreasingLabelings):
# class of increasing labelings with max entry n.
    
    def __classcall_private__(cls, n):
        """
        Standardize input to ensure a unique representation.

        EXAMPLES::

            sage: C = Compositions(5)
            sage: C2 = Compositions(int(5))
            sage: C3 = Compositions(ZZ(5))
            sage: C is C2
            True
            sage: C is C3
            True
        """
# Change doc test
        return super(IncreasingLabelings_n, cls).__classcall__(cls, Integer(n))

    def __init__(self, n):
        self.n = n

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(Compositions(3))
            'Compositions of 3'
        """
        return "Increasing labelings of %s with largest possible label $s. "%(self.n, self.poset)
# Fix doctests.    


    
def inc_labels(poset,restrict,i=None,f=None,labels=None):
    if labels is None:
        labels=[]
    if f is None:
        f={}
    if i is None:
        i=0
    if i<poset.cardinality():
        for j in restrict[i]:
            assign_val=True
            for k in poset.lower_covers(i):
                if j<=f.get(k):
                    assign_val=False
                    break
            if assign_val is True:
                f.update({i:j})
                labels=inc_labels(poset,restrict,i+1,f,labels)
                f.pop(i)
    else:
        if f not in labels:
            labels.append(copy(f))
    return labels


