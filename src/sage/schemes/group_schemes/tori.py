from __future__ import print_function, absolute_import



from sage.schemes.generic.scheme import Scheme
from sage.categories.morphism import SetMorphism
"""
from sage.modules.free_module import FreeModule_generic_pid
"""
from sage.matrix.action import MatrixMulAction
from sage.categories.action import PrecomposedAction
from six import integer_types
from sage.rings.integer_ring import ZZ
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.categories.map import Map
from sage.modules.free_module import FreeModule_ambient_pid
import sage.groups.libgap_group
import sage.matrix.matrix0
from sage.categories.homset import Hom
from sage.groups.matrix_gps.catalog import GL
from sage.libs.gap.element import GapElement
import sage.groups.libgap_wrapper 
import sage.groups.libgap_group
import sage.libs.gap.util
from sage.libs.gap.libgap import libgap
from sage.interfaces.gap import gap
from sage.matrix.constructor import matrix
from sage.combinat.permutation import Permutation
from sage.matrix.matrix_space import MatrixSpace


def is_GroupScheme(X):
    return isinstance(X,GroupScheme)



class GroupScheme(Scheme):
    """Group_scheme"""
    def __init__(self, X=None):
        Scheme.__init__(self,X)
        

















class GAPMap_toGLn(Map):
    """Get a Sage map from a GAP group homomoprhism with values in matrices  """
    def __init__(self, galois,rank, hom):
        Map.__init__(self,Hom(galois,GL(rank,ZZ)))
        self._morphism=hom
        self._domain=galois
        self._codomain=GL(rank,ZZ)
    def _call_(self,element):
        Elt=Permutation(element)
        return matrix((gap.Image(self._morphism,Elt)).sage())
        



class CharacterLattice(FreeModule_ambient_pid):
    """
    Generates the Character lattice  with Galois action
    If the parameter action is an integer, it will be the rank of the ZZ-module, seen with trivial Galois action. Otherwise, action should be an r-tuple of matrice where r is the number of generators or the group.



    """




    def __init__(self, galois, action, check=True):
        self._galois_group=galois
        self._generators=galois.gens()
        if action  in ZZ: 
            FreeModule_ambient_pid.__init__(self,ZZ,action)
            self._rank=action
            self._action_matrices=[matrix.identity(action) for i in range(len(self._generators))]
        elif not len(action) :
            raise ValueError('The module is missing the action of the Galois group.')
        else:
            self._action_matrices=action    
            for x in action: 
                if not x.nrows()==x.ncols():
                    raise ValueError('The matrices for the action need to be squares.')
                elif not x.nrows()==action[0].nrows(): 
                    raise ValueError('The matrices fot the action need to have same dimension.')
                elif not x in GL(action[0].nrows(),ZZ):
                    raise ValueError('The matrices need to be invertible in ZZ.')
                else :
                    self._rank=action[0].nrows()
                    FreeModule_ambient_pid.__init__(self,ZZ,self._rank)
                    self._action=action
            if not len(self._generators)==len(action):
                raise ValueError('The number of action matrices needs to match the number of generators of the Galois group.')
        G=libgap(galois)
        GenG=gap(self._generators)
        Mats=libgap(self._action_matrices)
        self._action_morphism=gap.GroupHomomorphismByImages(G,libgap.GL(self._rank,ZZ),GenG,Mats)

        if check and self._action_morphism==gap.fail:
            raise ValueError('The action is not well defined') 

        self._GAPMap=GAPMap_toGLn(self._galois_group,self._rank,self._action_morphism)
        MZ=MatrixSpace(ZZ,self._rank)
        A=MZ.get_action(self)
        self._action=PrecomposedAction(A,self._GAPMap,None)
    def _act(self,g,e):
        m=self._GAPMap
        mat=m._value(g)
        return mat * e










		
		
class AlgebraicTorus(GroupScheme):
    """
    Creates an algebraic torus through its equivalence of categories with the action of a Galois Group on an integral lattice.




    """

    def __init__(self, lattice):
        """
        Here G is the Galois group, which can be either an abstract group, or actual galois group, with corresponding base field and/or extension.

        The argument A is the action of this group, which has to be specified for each generator of G.
        """

        if lattice is None:
        	raise ValueError('You have to specify a lattice.')
        else:
            self._lattice = lattice


