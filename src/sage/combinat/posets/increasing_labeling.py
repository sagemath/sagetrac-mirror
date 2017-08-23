


# Will eventually need to figure out import statements
# ...making tentative list below

from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.rings.all import ZZ

from sage.rings.integer import Integer
from sage.combinat.combinatorial_map import combinatorial_map


class IncreasingLabeling(ClonableArray):
    
    def __classcall_private__()
    # Include any pre-parsing. Want base object to be list with respect to linear extension.
    # But want user to be able to give us a function or dictionary on poset elements
    # Without specifying a linear extension
    
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
            
   STUFFF

class IncreasingLabelings(UniqueRepresentation, Parent):
    def __classcall_private__()
    # Class call will redirect to the right parent based on input
    # Ie, we have one parent for increasing labelings with given function
    # one parent for largest max entry, and this is the 'parent' of the parents
    # that decides where things go
    
    def __init__(self, poset, is_infinite=False):
        """
        Initialize ``self``.

        """
        self._poset = poset

        if is_infinite:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=FiniteEnumeratedSets())

    def _element_constructor_(self, lst, check=True):
        if isinstance(lst, IncreasingLabeling):
            lst = list(lst)
        if not isinstance(lst, (list, tuple)):
            raise TypeError("input should be a list or tuple")
        lst = [self._poset(_) for _ in lst]
        if self._is_facade:
            return lst
        else:
            return self.element_class(self, lst, check)
        

    def __contains__(self, x):    

    Element = IncreasingLabeling

    
class IncreasingLabelings_all(IncreasingLabelings):
    def __classcall_private__(self, n=None, **kwargs):
        
    def __init__(self, poset, labels, is_infinite=True):
    
    
class IncreasingLabelings_restrict(IncreasingLabelings):
# Class of increasing labelings with function for range
# of possible values 'restrict'
        self._labels = labels

    def __init__(self, restrict):
        """
        Initialize ``self``.

        """
        self.restrict=restrict
        IncreasingLabelings.__init__(self, False)
        
        
    def __iter__(self):
        for il in inc_labels(self._poset, self.restrict)
            yield self.element_class(self, il)
        

    
class IncreasingLabelings_q(IncreasingLabelings):
# class of increasing labelings with max entry q.
    
    def __classcall_private__(cls, n):
        
    def __init__(self, n):
    


    
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


