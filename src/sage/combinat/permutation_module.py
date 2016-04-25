from sage.combinat.free_module import CombinatorialFreeModule

#from sage.combinat.partitions import Partition
from sage.combinat.tableau import Tableau,Tableaux,StandardTableaux
from sage.combinat.tabloid import Tabloid,Tabloids



class PermutationModule(CombinatorialFreeModule):

    def __init__(self, R, shape):
        """
        TESTS::

        """
#        shape = sage.combinat.partitions.Partition(shape)
#   this is a cludge, will -  ideally we will override other methos
#  to give right representation name
#  and have elements have a convienient method that is shortcut for
# current use of map_support 
        tabloids = Tabloids(shape)
        super(PermutationModule,self).__init__(R,tabloids)
  
#        self._shape = Partition(shape)
#        self._n = self._shape.size()
#        self._group = Permutations(self._n)
#        self._tabloids = Tabloids(self._shape)
#        CombinatorialFreeModule.__init__(self, R, tabloids)


class SpechtModule(PermutationModule):
    
    def __init__(self,R,shape):
    
    
        tabloids= StandardTabloids(shape)
        super(SpechtModule,self).__init__(R,tabloids)
