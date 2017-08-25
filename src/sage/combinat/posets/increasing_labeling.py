


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
from sage.sets.family import Family
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass

from sage.rings.integer import Integer
from sage.combinat.combinatorial_map import combinatorial_map

@add_metaclass(InheritComparisonClasscallMetaclass)
class IncreasingLabeling(ClonableArray):
    
    def check(self):
        P = self.parent()._poset
        is_inc_label = True
        for elem in P:
            for lc in P.lower_covers(elem):
                if self[lc]>=self[elem]:
                    is_inc_label = False
        if is_inc_label is False:
            raise ValueError("%s is not an increasing labeling of %s"%(self, P))

class IncreasingLabelings(UniqueRepresentation, Parent):
#    def __classcall_private__(self, n=None, **kwargs)
    # Class call will redirect to the right parent based on input
    # Ie, we have one parent for increasing labelings with given function
    # one parent for largest max entry, and this is the 'parent' of the parents
    # that decides where things go
    @staticmethod
    def __classcall_private__(self, poset=None, n=None, restrict=None):
#        self._poset=poset
        if n is None:
            if restrict is None:
                return IncreasingLabelings_all(poset=poset)
            else:
                return IncreasingLabelings_restrict(poset=poset,restrict=restrict)     
        else:
            if restrict is None:
                if isinstance(n, (int,Integer)):
                    return IncreasingLabelings_n(poset=poset,n=n)
        raise ValueError("At least one of restrict and n must be None")

    def __init__(self, is_infinite=False):
        if is_infinite:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=FiniteEnumeratedSets())

#    def _element_constructor_(self, label, check=True):
#        if isinstance(label, IncreasingLabeling):
#            label = Family(label)
#        if not isinstance(label, (list, dict)):
#            raise TypeError("input should be a list")
#        label2 = Family({self._poset(x):label[x] for x in label})
#        return self.element_class(self, label2) 

    Element = IncreasingLabeling

    
class IncreasingLabelings_all(IncreasingLabelings):
    def __init__(self, poset):
        IncreasingLabelings.__init__(self, True)
        self._poset=poset
        

    def _repr_(self):
        return "The set of all increasing labelings of %s"%(self._poset)

    def _element_constructor_(self, label, check=True):
        if isinstance(label, IncreasingLabeling):
            label = Family(label)
        if not isinstance(label, (list, dict)):
            raise TypeError("input should be a list")
        label2 = Family({self._poset(x):label[x] for x in label})
        return self.element_class(self, label2) 

class IncreasingLabelings_n(IncreasingLabelings):
    def __init__(self,poset, n):
        IncreasingLabelings.__init__(self, False)
        self._n = n
        self._poset=poset

    def _repr_(self):
        return "The set of all increasing labelings of %s with largest possible part %s"%(self._poset, self._n)

    def _element_constructor_(self, label, check=True):
        if isinstance(label, IncreasingLabeling):
            label = Family(label)
        if not isinstance(label, (list, dict)):
            raise TypeError("input should be a list")
        label2 = Family({self._poset(x):label[x] for x in label})
        for p in label2:
            if p>self._n:
                raise ValueError("%s exceeds largest allowable label %s"%(p,self._n))
        return self.element_class(self, label2,)
            
class IncreasingLabelings_restrict(IncreasingLabelings):
    def __init__(self,poset, restrict):
       IncreasingLabelings.__init__(self, False)
       self._restrict = restrict
       self._poset=poset 

    def _repr_(self):
        return "The set of all increasing labelings of %s with parts restricted by %s"%(self._poset, self._restrict)

    def _element_constructor_(self, label, check=True):
        if isinstance(label, IncreasingLabeling):
            label = Family(label)
        if not isinstance(label, (list, dict)):
            raise TypeError("input should be a list")
        label2 = Family({self._poset(x):label[x] for x in label})
        for p in self._poset:
            if label2[p] not in restrict[p]:
                raise ValueError("%s is not a valid label for %s"%(label2[p],p))
        return self.element_class(self, label2,)
 
def inc_labels_iter(poset,restrict,i=None,f=None,labels=None):
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
                for il in inc_labels_iter(poset,restrict,i+1,f,labels):
                    yield il
                f.pop(i)
    else:
        if f not in labels:
            labels.append(copy(f))
            yield f
