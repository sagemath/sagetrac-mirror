from sage.structure.sage_object import SageObject
#from sage.misc.all import cached_method
from sage.rings.integer import Integer
from sage.combinat.permutation import Permutation
from sage.combinat.root_system.weyl_group import WeylGroup, WeylGroup_gens, WeylGroupElement
from sage.graphs.all import DiGraph
from sage.combinat.cluster_algebra_quiver.double_bruhat_digraph import DoubleBruhatDigraph
from sage.matrix.matrix import Matrix

class QuiverConstructionType(SageObject):

    def __init__(self,data):
        
        from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        from sage.combinat.cluster_algebra_quiver.cluster_seed import ClusterSeed
        
        self._description = None
        self._details = ""
        
        # The standard construction types
        if isinstance(data,ClusterSeed):
            self._description = "a cluster seed"
        elif isinstance(data,ClusterQuiver):
            self._description = "a quiver"
        elif isinstance(data,DiGraph):
            self._description = "a digraph"
        elif isinstance(data,Matrix):
            self._description = "a skew-symmetrizable matrix"
            
        # More exotic construction types
        elif type(data) in [list,tuple]:
            if data[0] == 'DB':
                
                # defaults to type A and interprets the list as a partition
                if len(data[1]) == 2:
                    try: 
                        u = Permutation(data[1][0])
                        v = Permutation(data[1][1])
                    except Exception:
                        print "Invalid input for construction type 'DB'"
                    CartanType = ['A',len(u)]
                    
                elif len(data[1]) == 3:
                    if isinstance(data[1][0],WeylGroup_gens):
                        CartanType = data[1][0].cartan_type()
                        
                    # In this case, the list is interpreted as a word in the generators of the Weyl group
                    elif type(data[1][0]) in [tuple,list]:
                        
                        if len(data[1][0]) == 2 and isinstance(data[1][0][0], str) and type(data[1][0][1]) in [int,Integer]:
                            CartanType = data[1][0]
                        else:
                            _construction_type_error(data)
                        W = WeylGroup(CartanType)
                    if isinstance(data[1][1],list) and isinstance(data[1][2],list):
                        S = W.simple_reflections()
                        u = 1
                        v = 1
                        for i in range(len(data[1][1])):
                            u = u*S[data[1][1][i]]
                            v = v*S[data[1][2][i]]
                    elif isinstance(data[1][1], WeylGroupElement) and isinstance(data[1][2], WeylGroupElement):
                        u = data[1][1]
                        v = data[1][2]
                    else:
                        try:
                            u = WeylGroupElement(WeylGroup(CartanType),data[1][1])
                            v = WeylGroupElement(WeylGroup(CartanType),data[1][2])
                        except Exception:
                            print "Invalid input for construction type 'DB'"
                        
            
                self._construction_digraph, self._frozen = DoubleBruhatDigraph(CartanType,u,v)
                
                # Modifies the original construction so that the digraph has skew-symmetric edges
                newEdges = []
                for edge in self._construction_digraph.edges():
                    if edge[0] in self._frozen or edge[1] in self._frozen:
                        newEdges.append((edge[0],edge[1],(1,-1)))
                    else:
                        newEdges.append(edge)
                self._digraph = DiGraph(newEdges)
                u = u.reduced_word()
                v = v.reduced_word()
                self._description = "a double Bruhat cell in a group of type " + str(CartanType)
                self._details = " using the reduced words " + str(u) + " and " +  str(v)
                
                
        if not self._description:
            self._description = str(type(data))
        
    def _repr_(self):
        """
        Return the string representation of ``self``.
        """
        
        return self._description
    
    def full_description(self):
        
        print "Quiver constructed from " + self._description + self._details  + "." 
    
def _construction_type_error(data):
    return ValueError(data)
    
    