from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.combinat.permutation import Permutation
from sage.combinat.root_system.weyl_group import WeylGroup, WeylGroup_gens, WeylGroupElement
from sage.combinat.cluster_algebra_quiver.double_bruhat_digraph import DoubleBruhatDigraph

class QuiverConstructionType(SageObject):
    def __init__(self,data):
        
        self._description = str(type(data))
        self._details = ""
        if type(data) in [list,tuple]:
            if data[0] == 'DB':
                if len(data[1]) == 2:
                    try: 
                        u = Permutation(data[1][0])
                        v = Permutation(data[1][1])
                    except Exception:
                        print "Invalid input for construction type 'DB'"
                    CartanType = ['A',len(u)]
                elif len(data[1]) == 3:
                    if isinstance(data[1][0],WeylGroup_gens):
                        CartanType = data[1][0].cartan_type
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
                        else :
                            try:
                                u = WeylGroupElement(WeylGroup(CartanType),data[1][1])
                                v = WeylGroupElement(WeylGroup(CartanType),data[1][2])
                            except Exception:
                                print "Invalid input for construction type 'DB'"
                        
            
                self._digraph, self._frozen = DoubleBruhatDigraph(CartanType,u,v)
                u = u.reduced_word()
                v = v.reduced_word()
                self._description = "a double Bruhat cell in a group of type " + str(CartanType)
                self._details = " using the reduced words " + str(u) + " and " +  str(v)
        
        
    def _repr_(self):
        """
        Return the string representation of ``self``.
        """
        
        return self._description
    
    def full_description(self):
        
        print "Quiver constructed from " + self._description + self._details  + "." 
    
def _construction_type_error(data):
    return ValueError(data)
    
    