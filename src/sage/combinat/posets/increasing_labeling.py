


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
    @staticmethod
    def __classcall_private__(self, poset=None, n=None, restrict=None):
        self.poset=poset
        if n is None:
            if restrict is None:
                return IncreasingLabelings_all(poset)
            else:
                return IncreasingLabelings_restrict(poset,restrict=restrict)     
        else:
            if restrict is None:
                if isinstance(n, (int,Integer)):
                    return IncreasingLabelings_n(n=n)
        raise ValueError("At least one of restrict and n must be None")

    def __init__(self, is_infinite=False):
        if is_infinite:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=FiniteEnumeratedSets())
        

    Element = IncreasingLabeling

    
class IncreasingLabelings_all(IncreasingLabelings):

    def __init__(self):
        IncreasingLabelings.__init__(self, True)

class IncreasingLabelings_n(IncreasingLabelings):

    def __init__(self, n):
        IncreasingLabelings.__init__(self, False)
        self.n = n
            
class IncreasingLabelings_restrict(IncreasingLabelings):

   def __init__(self, restrict):
       IncreasingLabelings.__init__(self, False)
       self.restrict = restrict
 
