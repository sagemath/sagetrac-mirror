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
from sage.modules.free_module import FreeModule_generic
from sage.modules.free_module import FreeModule_submodule_pid
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
from sage.groups.matrix_gps.finitely_generated import MatrixGroup
from sage.groups.perm_gps.permgroup import load_hap
from sage.matrix.special import block_matrix



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
        if element is tuple:
            Elt=Permutation(element)
        else:
            Elt=element
        return matrix((gap.Image(self._morphism,Elt)).sage())
        

class Lattice_generic(FreeModule_generic):




    def __init__(self, galois, action=1, check=True):
        if type(galois)==list :
            self._group=MatrixGroup(galois)
            self._generators=galois
            self._rank = galois[0].nrows()
            self._action_matrices=galois
            FreeModule_generic.__init__(self,ZZ,self._rank,self._rank)
            A=libgap(galois)
            G=gap.GroupByGenerators(A)
            self._action_morphism=gap.GroupHomomorphismByImages(G,G,A,A)
        else:
            self._group=galois
            self._generators=galois.gens()
            if action  in ZZ: 
                FreeModule_generic.__init__(self,ZZ,action,action)
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
                        FreeModule_generic.__init__(self,ZZ,self._rank,self._rank)
                if not len(self._generators)==len(action):
                    raise ValueError('The number of action matrices needs to match the number of generators of the Galois group.')
            G=gap(self._group)
            GenG=gap(self._generators)
            Mats=gap(self._action_matrices)
            if self._rank==1 : 
                gl=gap.Group([ [ [ -1 ] ] ])
            else :
                gl=gap.GL(self._rank,ZZ)
            self._action_morphism=gap.GroupHomomorphismByImages(G,gl,GenG,Mats)
            if check and gap.IsBool(self._action_morphism) :
                raise ValueError('The action is not well defined') 

        self._GAPMap=GAPMap_toGLn(self._group,self._rank,self._action_morphism)
        MZ=MatrixSpace(ZZ,self._rank)
        A=MZ.get_action(self)
        self._action=PrecomposedAction(A,self._GAPMap,None)
    def _act(self,g,e):
        return  self._action.act(g,e)






    ##### METHODS FOR LATTICES WITH GROUP ACTIONS#######
    def group(self):
        return self._group

    def isomorphic_ambient_lattice(self):
        return self.ambientify()

    def GAPMatrixGroup(self):
        f=self._GAPMap
        return gap.Image(f._morphism)

    def TateCohomology(self,n):
        return self.TCohomology(n)

    def induced_lattice(self,group,build=True):
        return self.induction(group,build)

    def zero_sum_sublattice(self):
        return self.zero_sublattice()
    def norm_one_restriction(self,group,build=True):
        return self.norm_one_restriction_of_scalars(self,group,build)









class Lattice_ambient(FreeModule_ambient_pid,Lattice_generic):
    """
    Generates the Character lattice  with Galois action
    If the parameter action is an integer, it will be the rank of the ZZ-module, seen with trivial Galois action. Otherwise, action should be an r-tuple of matrice where r is the number of generators or the group.


    """

    def __init__(self, galois, action=1, check=True):
        Lattice_generic.__init__(self,galois,action,check)
        FreeModule_ambient_pid.__init__(self,ZZ,self._rank)


    def ambientify(self):
        return self




    def colattice(self):
        newacts = [ i.inverse().transpose() for i in self._action_matrices]
        return Lattice_ambient(self._group,newacts)



    def TCohomology(self,n):
        
        MG=self.GAPMatrixGroup()
        G=libgap(self._group) if type(self._group[0]) is PermutationGroupElement else gap.Group(self._action_matrices)
        #This is a slight adaptation of the code from Hoshi and Yamasaki done on GAP
        if n==0 :
            M=matrix.zero(self._rank)
            if type(self._group[0]) is PermutationGroupElement :
                Lst=[gap(i) for i in self._group]
            else:
                Lst=[i.gap() for i in self._group]
            for i in Lst:
                M+=matrix(gap.Image(self._action_morphism,i).sage())
            #M=matrix((libgap.Sum(MG)).sage())
            S=M.smith_form(False,True)
            R=S.rank()
            RR=[S[i][i] for i in range(R)]
            return [i for i in RR if i>1]
        elif n==-1 :
            m=gap([])
            for i in gap.GeneratorsOfGroup(MG) : 
                m=gap.Concatenation(m,i-gap.Identity(MG))
            ms=matrix(m.sage())
            s=ms.smith_form(False,True)
            r=s.rank()
            rr=gap([s[i][i] for i in range(r)])
            return [i for i in rr if i>1]
        else  :
            load_hap()
            if n>0 :
                #This computes the standard resolution of G in HAP
                R = gap.ResolutionFiniteGroup(G,n+1)
                #Then applies the map to the action to the resolution
                TR = gap.HomToIntegralModule(R,self._action_morphism)
                #Might have a problem because gap does only right actions ?
                return (gap.Cohomology(TR,n))
            else :
                R = gap.ResolutionFiniteGroup(G,-n)
                TR = gap.TensorWithIntegralModule(R,self._action_morphism)
                return (gap.Homology(TR,-n-1)).sage()



    #.self._induction(containing_group,build) computes the lattice obtained by induction from 
    #the group defining the lattice to a bigger group "containing_group" 


    #The build argument will ensure that we build the lattice if it is true, 
    #otherwise it just gives us the matrices needed to build the lattice with
    #the group we input, this will be useful to get other constructions like 
    #the character lattice of the norm 1 restriction of scalars torus



    def induction(self,group,build=True):
        LCos=group.cosets(self._group,'left')
        LCosReps=[i[0] for i in LCos]
        LCosnum=len(LCos)
    #get_coset_index will take an element g of G and return the index of the coset it belongs to
        def get_coset_index(g):
            for i in range(LCosnum):
                if g in LCos[i]:
                    return i
            raise ValueError('Found an element which does not belong to any coset')

    #decomp takes an element g of the big group G, and returns a pair
        def decomp(g):
            i=get_coset_index(g)
            return [i,(LCosReps[i].inverse())*g]
        r"""
        decomp2(element,index) will do the main computation for the group action on the 
        induced lattice (we deal with finite groups so induced is the same as coinduced)
        -\ZZ[G]\otimes_{\ZZ[H]} L -

         To do so, we look at the multiplication g*gi where gi is a representative of the
         ith left coset
        """

        def decomp2(g,i):
            return decomp(g*LCosReps[i])


        #Now we create block matrices corresponding to actions of each generator of G, it 
        #is the composition of a block unimodular matrix for permutation of blocks
        #and the block diagonal matrix applying the correct element of H to each block

        def make_matrix(g):
            Lst=[]
            Bigmatlist = []
            for i in range(LCosnum) : 
                Lst.append(decomp2(g,i))
            Matlist=[matrix.zero(self._rank) for j in range(LCosnum)]
            Matlist=[ matrix((gap.Image(self._action_morphism,gap(j[1]))).sage())   for j in Lst]
            for i in range(LCosnum):
                Bigmatlist+=[Matlist[j] if Lst[j][0]==i else matrix.zero(self._rank)   for j in range(LCosnum)]
            return block_matrix(LCosnum,LCosnum,Bigmatlist)

  

        GGen=self._generators




    



        gens_induced_act = [make_matrix(g) for g in group.gens()]
        if build :
            return Lattice_ambient(group,gens_induced_act)   
        else:
            return gens_induced_act
    #now we build the character lattice for the norm 1 restrictions of scalars        
    def zero_sum_sublattice(self,ambient=False):
        if ambient:
            A=self
            acts=self._action_matrices
            newacts=[A[[i for i in range(A.ncols()-1)],[i for i in range(A.ncols()-1)]]-A[[i for i in range(A.ncols()-1)],[A.ncols()-1]]*matrix(1,[1 for i in range(A.ncols()-1)]) for A in acts]
            return Lattice_ambient(self._group,newacts)
        else :
            oldBasis=self.basis()
            newBasis = [ oldBasis[i]-oldBasis[len(oldBasis)-1]     for i in range(len(oldBasis)-1)  ]
            return SubLattice(self,newBasis)
        r"""
        This does the following, for each nxn matrix defining the action of the induced representation, 
        it extracts two blocks
        (        |   )
        (    A   | B )     A which is of size n-1 x n-1
        (        |   ) and B which is of size n-1 x 1 
        (------------)
        (    C   | D )

        Then it computes the action for the norm 1 restriction of scalars, which is
        A - B * (1,1,...,1)


        Indeed, since the norm 1 restriction lattice will be the sublattice 
        of elements summing up to 1, so we only consider the first n-1 basis elements, 
        doing the matrix computation, the element l=(l_1,...l_(n-1)) is sent to 
        A *l + B * (-sum_i l_i), which is the transformation A-B*(1,...,1)
        """

    def norm_one_restriction_of_scalars(self,group,build=True):
        return self.zero_sum_sublattice(self.induction(group))


 


class SubLattice(Lattice_generic,FreeModule_submodule_pid):
    def __init__(self, lattice,basis,check=True):
        Lattice_generic.__init__(self,lattice._group,lattice._action_matrices)
        FreeModule_submodule_pid.__init__(self,lattice,basis)
        self._parent_lattice=lattice
        if check:
            for i in lattice._group.gens():
                for j in basis:
                    if not lattice._act(i,j) in self:
                        raise ValueError("The basis is not stable under the action of the group")

    def ambientify(self):
        oldBasis=self.basis()
        act_builder=[]
        for g in self._group.gens() :
            mat_builder=[]
            for i in oldBasis:
                mat_builder.append(self.coordinate_vector(self._act(g,i)))
            act_builder.append(matrix(mat_builder))
        return Lattice_ambient(self._group,act_builder)

    def TCohomology(self,n):
        return self.ambientify().TCohomology(n)

    def induction(self,group):
        return self.ambientify().induction(group)

    def norm_one_restriction_of_scalars(self,group):
        return self.ambientify().induction(group)
    


        """
        TO FINISH
    def zero_sum_sublattice(self):
        oldBasis=self.basis()
        newBasis=[]
        diagonal=0
        for i in self._parent_lattice.basis():
            diagonal+=i
        for i in self.basis():
            newBasis.append()
        return 1
        """



		
		
class AlgebraicTorus(GroupScheme):
    """
    Creates an algebraic torus through its equivalence of categories with the action of a Galois Group on an integral lattice.




    """

    def __init__(self, lattice):
        """
        TO DO

        """
        GroupScheme.__init__(self)
        if lattice is None:
        	raise ValueError('You have to specify a lattice.')
        else:
            self._lattice = lattice
    def _repr_(self):
        return "Algebraic Torus of rank %s defined by the following lattice :\n"%(self.rank())+self._lattice._repr_()+"\nand an action by the galois group of the form :\n"+self._lattice._group._repr_() 
    def Tate_Cohomology(self,n):
        #This gives the isomorphism type of the nth cohomology group using Tate-Nakayama duality 
        return self._lattice.TateCohomology(2-n)
    def galois_group(self):
        return self._lattice._group
    def character_lattice(self):
        return self._lattice
    def rank(self):
        return self._lattice._rank
    def cocharacter_lattice(self):
        return self._lattice.colattice()

        #gives the torus representing the Restriction of scalars. 
        #Right now, for a torus defined over K, splitting over L, 
        #to compute the restriction of scalars to k inside K, 
        #the user has to enter the galois group of the extension L/k
        #In the future, when we will have a better notion for Galois group
        #perhaps we can deal with fields directly.


    def restriction_of_scalars(self,group):
        return AlgebraicTorus(self._lattice.induced_lattice(group))
    def norm_one_restriction(self,group):
        return AlgebraicTorus(self._lattice.norm_one_restriction(group))