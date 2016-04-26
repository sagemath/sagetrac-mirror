from sage.combinat.free_module import CombinatorialFreeModule

#from sage.combinat.partitions import Partition
#from sage.combinat.tableau import Tableau,Tableaux,StandardTableaux




class PermutationModule(CombinatorialFreeModule):
    def __init__(self, R, shape):
        from sage.combinat.tabloid import Tabloid,Tabloids
        """
        TESTS::

        """
#        shape = sage.combinat.partitions.Partition(shape)
#   this is a cludge, will -  ideally we will override other methos
#  to give right representation name
#  and have elements have a convienient method that is shortcut for
# current use of map_support and lambda function

        tabloids = Tabloids(shape)
        super(PermutationModule,self).__init__(R,tabloids)
  
#        self._shape = Partition(shape)
#        self._n = self._shape.size()
#        self._group = Permutations(self._n)
#        self._tabloids = Tabloids(self._shape)
#        CombinatorialFreeModule.__init__(self, R, tabloids)


class SpechtModule(CombinatorialFreeModule):
    
    def __init__(self,R,shape):
        from sage.combinat.tableau import StandardTableaux
        from sage.sets.set import Set
        SYT= StandardTableaux(shape)
        basis=Set([])
        for T in SYT:
            basis+=Set([T.to_polytabloid()])
            
        super(SpechtModule,self).__init__(R,basis) 
        
        
        
