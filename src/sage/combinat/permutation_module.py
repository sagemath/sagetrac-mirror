from sage.combinat.free_module import CombinatorialFreeModule

#from sage.combinat.partitions import Partition
#from sage.combinat.tableau import Tableau,Tableaux,StandardTableaux

class PermutationModule(CombinatorialFreeModule):
    def __init__(self, R, shape):
        from sage.combinat.partition import Partition, Partitions
        from sage.combinat.tabloid import Tabloids
        """
        TESTS::

        """
        self.R=R
        
        if Partitions().is_parent_of(shape)==True:
            self.shape=shape
        else:
            self.shape=Partition(shape)            
#   this is a kludge, will -  ideally we will override other methos
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

    def to_specht_submodule(self,shape=0):
        if shape==0:
            shape=self.shape
#WARNING:  this is not fully implement, and results may not be consistent
        return self.submodule(SpechtModule(self.R,shape).basis())    
        
        
        
class SpechtModule(CombinatorialFreeModule):
    # this could greatly benefit by defining classes for the set of polytabloids and polytabloid elements
    def __init__(self,R,shape):
        from sage.combinat.tableau import StandardTableaux
        from sage.sets.set import Set
        SYT= StandardTableaux(shape)
        basis=Set([])
        for T in SYT:
            basis+=Set([T.to_polytabloid()])
            
        super(SpechtModule,self).__init__(R,basis) 
