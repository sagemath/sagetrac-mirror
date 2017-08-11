r"""
Quiver Construction Types

This file contains the QuiverConstructionType class, which records the method used to 
construct a quiver. This also provides infrastructure to implement more exotic quiver
constructions, such as the construction of a quiver from a double Bruhat cell.

AUTHORS:

- Ben Strasser (2017, initial version)
"""
from sage.structure.sage_object import SageObject
#from sage.misc.all import cached_method
from sage.rings.integer import Integer
from sage.combinat.permutation import Permutation
from sage.combinat.root_system.weyl_group import WeylGroup, WeylGroup_gens, WeylGroupElement
from sage.graphs.all import DiGraph
from sage.combinat.cluster_algebra_quiver.double_bruhat_digraph import DoubleBruhatDigraph
from sage.matrix.matrix import Matrix
from six.moves import range

class QuiverConstructionType(SageObject):
    """
    The *construction type* for a *quiver*. In the case of standard construction types,
    this class simply records the method used to construct the quiver and stores the 
    string representation of that quiver. For other construction types, this class 
    constructs the digraph and the list of frozen vertices. It also records the string
    representation of the quiver, with details about the construction method used.

    TESTS::
    sage: W = WeylGroup(['B',5])
    sage: s1,s2,s3,s4,s5 = W.simple_reflections()
    sage: ClusterQuiver(['DB',[W, s1*s2*s3,s1*s2*s3*s4*s5]])._construction_type
    a double Bruhat cell in a group of type ['B', 5]
    sage: ClusterQuiver(['DB',[W, s1*s2*s3,s1*s2*s3*s4*s5]])._construction_type.full_description()
    Quiver constructed from a double Bruhat cell in a group of type ['B', 5] using the reduced words [1, 2, 3] and [1, 2, 3, 4, 5].

    sage: ClusterQuiver(['A',3])._construction_type
    a digraph

    sage: S = ClusterSeed(['D',4])
    sage: S._construction_type
    a digraph
    sage: ClusterSeed(S)._construction_type
    a digraph
    sage: M = S.b_matrix()
    sage: ClusterSeed(M)._construction_type
    a skew-symmetrizable matrix
    sage: S._construction_type = None
    sage: ClusterQuiver(S)._construction_type
    a cluster seed
    sage: Q = ClusterQuiver(S)
    sage: Q._construction_type = None
    sage: ClusterQuiver(Q)._construction_type
    a quiver
    """

    def __init__(self,data):
        
        
        from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        from sage.combinat.cluster_algebra_quiver.cluster_seed import ClusterSeed
        
        self._description = None
        self._details = ""
        self._digraph = False
        
        # The standard construction types

        if isinstance(data,ClusterSeed):
            self._description = "Seed"
        elif isinstance(data,ClusterQuiver):
            self._description = "Quiver"
        elif isinstance(data,DiGraph):
            self._description = "Digraph"
        elif isinstance(data,Matrix):
            self._description = "Matrix"
            
        # More exotic construction types
        elif type(data) in [list,tuple]:
            if data[0] == 'DB':
                self._description = "DB"
                # defaults to assume input does not contain a two reduced words in 'list form'
                word = False
                
                # defaults to type A and interprets the list as a partition
                if len(data[1]) == 2:
                    try: 
                        u = Permutation(data[1][0])
                        v = Permutation(data[1][1])
                    except Exception:
                        print("Invalid input for construction type 'DB'")
                    CartanType = ['A',len(u)]
                    
                # The first element of data[1] should identify the relevant Weyl group
                elif len(data[1]) == 3:
                    if isinstance(data[1][0],WeylGroup_gens):
                        CartanType = data[1][0].cartan_type()
                        
                    elif type(data[1][0]) in [tuple,list]:
                        
                        if len(data[1][0]) == 2 and isinstance(data[1][0][0], str) and type(data[1][0][1]) in [int,Integer]:
                            CartanType = data[1][0]
                        else:
                            _construction_type_error(data)
                        W = WeylGroup(CartanType)

                    # In this case, the list is interpreted as a word in the generators of the Weyl group
                    if isinstance(data[1][1],list) and isinstance(data[1][2],list):
                        S = W.simple_reflections()
                        u = W.one()
                        v = W.one()
                        for i in data[1][1]:
                            u = u*S[i]
                        for j in data[1][2]:
                            v = v*S[j]
                        
                        # If the expressions for u and v were reduced, use them as is
                        if u.length() == len(data[1][1]) and v.length() == len(data[1][2]):
                            word = data[1][1]
                            for n in data[1][2]:
                                word.append(-n)
                        

                    elif isinstance(data[1][1], WeylGroupElement) and isinstance(data[1][2], WeylGroupElement):
                        u = data[1][1]
                        v = data[1][2]
                    else:
                        try:
                            u = WeylGroupElement(WeylGroup(CartanType),data[1][1])
                            v = WeylGroupElement(WeylGroup(CartanType),data[1][2])
                        except Exception:
                            print("Invalid input for construction type 'DB'")

                        
            
                self._construction_digraph, self._frozen, self._strings = DoubleBruhatDigraph(CartanType,u,v,word)
                
                # Modifies the original construction so that the digraph has skew-symmetric edges
                newEdges = []
                for edge in self._construction_digraph.edges():
                    if edge[0] in self._frozen or edge[1] in self._frozen:
                        newEdges.append((edge[0],edge[1],(1,-1)))
                    else:
                        newEdges.append(edge)
                self._digraph = DiGraph(newEdges)
                # Ensures isolated frozen vertices are in the digraph
                self._digraph.add_vertices(self._frozen)
                u = u.reduced_word()
                v = v.reduced_word()
                

                self._cartan_type = CartanType        
                self._details = " using the reduced words " + str(u) + " and " +  str(v)
                
                
        if not self._description:
            self._description = str(type(data))
        
    def _repr_(self):
        """
        Return the string representation of ``self``.
        """
        if self._description == "Seed":
            return "a cluster seed"
        elif self._description == "Quiver":
            return "a quiver"
        elif self._description == "Digraph":
            return "a digraph"
        elif self._description == "Matrix":
            return "a skew-symmetrizable matrix"
        elif self._description == "DB":
            return "a double Bruhat cell in a group of type " + str(self._cartan_type)     


    
    def full_description(self):
        """
        Prints a string with detailed information about the construction_type.

        EXAMPLES::

        sage: ClusterQuiver(['DB',[['D',4],[1,2,3,4],[4,3,2,1]]])._construction_type.full_description()
        Quiver constructed from a double Bruhat cell in a group of type ['D', 4] using the reduced words [1, 2, 4, 3] and [4, 3, 2, 1].
        """
        print "Quiver constructed from " + str(self) + self._details  + "." 
    
def _construction_type_error(data):
    return ValueError(data)
    
    
