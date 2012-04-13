from sage.modular.btquotients.btquotient import *
from sage.modular.btquotients.ocmodule import *
#########################################################################
#       Copyright (C) 2011 Cameron Franc and Marc Masdeu
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################
from collections import namedtuple
from sage.structure.element import Element
from sage.structure.element import ModuleElement
from sage.modules.module import Module
from sage.rings.all import Integer
from sage.structure.element import Element
from sage.matrix.constructor import Matrix, zero_matrix
from sage.matrix.matrix_space import MatrixSpace_generic
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.all import Qp
from sage.rings.all import RationalField
from sage.rings.number_field.all import NumberField
from copy import copy
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.modular.hecke.all import (AmbientHeckeModule, HeckeSubmodule, HeckeModuleElement)
import sage.rings.arith as arith
import sage.modular.hecke.hecke_operator


class HarmonicCocycleElement(HeckeModuleElement):
    r"""
    Objects of this type are Gamma-invariant harmonic cocycles on the 
    Bruhat-Tits tree. Gamma-invariance is necessary so that the cocycle
    can be stored in terms of a finite amount of data.

    More precisely, given a BTQuotient T, we store harmonic cocycles as 
    a list of values in some coefficient module (e.g. for weight 2 forms
    can take Cp) indexed by edges of a fundamental domain for T in the 
    Bruhat-Tits tree. Evaluate the cocycle at other edges using Gamma
    invariance (although the values may not be equal over an orbit of
    edges as the coefficient module action may be nontrivial).

    INPUT:
    fixme: describe these???
     - ``vec`` - (default: None) 

     - ``from_values`` -  (default: False) 

    EXAMPLES:

    ::

    AUTHORS:


    - Cameron Franc (2012-02-20)
    - Marc Masdeu
    """
    def __init__(self,_parent,vec=None,from_values=False):
        HeckeModuleElement.__init__(self,_parent,None)
        self._parent=_parent
        if from_values:
            self._R=_parent._U._R
            self._wt=_parent._k
            self._nE=len(_parent._E)
            self._F=copy(vec)
            return
        if isinstance(vec, self.__class__):
            #The first argument is just a modular form that we copy as is
            self._R=vec._R
            self._nE=vec._nE
            self._wt=vec._wt
            self._F=[_parent._U.element_class(_parent._U,vec._F[ii],quick=True) for ii in range(self._nE)]
            return
        # When vec contains coordinates for the basis
        self._R=self._parent._R
        try:
            v=[self._R(x) for x in vec.list()]
        except AttributeError:
            v=[self._R(vec) for ii in range(self._parent.dimension())]
        self._wt=_parent._k
        self._nE=len(_parent._E)
        vmat=Matrix(self._R,1,self._parent.dimension(),v)
        tmp=(vmat*self._parent.ambient_module().basis_matrix()).row(0)
        self._F=[_parent._U.element_class(_parent._U,Matrix(self._R,self._wt-1,1,tmp[e*(self._wt-1):(e+1)*(self._wt-1)]),quick=True) for e in range(self._nE)]
        return

    def _add_(self,g):
        r"""
        This function adds two cocycles componentwise.
        """
        #Should ensure that self and g are modular forms of the same weight and on the same curve
        C=self.__class__
        return C(self._parent,self.element()+g.element())
    
    def _sub_(self,g):
        r"""
        This function computes the difference of two cocycles.
        """
        #Should ensure that self and g are modular forms of the same weight and on the same curve
        C=self.__class__
        return C(self._parent,self.element()-g.element())

    def _rmul_(self,a):
        r"""
        This function multiplies a cocycle by a scalar.
        """
        #Should ensure that 'a' is a scalar
        C=self.__class__
        return C(self._parent,a*self.element())


    def _repr_(self):
        r"""
        Retuns a string representing the values of the cocycle on the edges.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        tmp='Element of '+str(self._parent)
        return tmp

    def __eq__(self,other):
        r"""
        Test for equality with another cocycle. Two cocycles are equal if they
        take the same values on all the edges in a fundamental domain.

        EXAMPLES:

        This example illustrates ...

        ::
        """
        return all([self._F[e].__eq__(other._F[e]) for e in range(self._nE)])

    def __ne__(self,other):
        r"""
        Test for non-equality with another cocycle. Two cocycles are non-equal if they
        take the different values on at least one edge in a fundamental domain.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return any([self._F[e].__ne__(other._F[e]) for e in range(self._nE)])

    def __nonzero__(self):
        r"""
        Test for being non-zero.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return any([self._F[e].__nonzero__() for e in range(self._nE)])

    def valuation(self):
        r"""
        Returns the valuation of the cocycle, defined as the minimum of the values
        it takes on a set of representatives.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        if not self.__nonzero__():
            return oo
        else:
            return min([self._F[e].valuation() for e in range(self._nE)])

    def _compute_element(self):
        r"""
    
        """
        R=self._R
        A=self._parent.basis_matrix().transpose()
        B=Matrix(R,self._nE*(self._parent._k-1),1,[self._F[e]._val[ii,0]  for e in range(self._nE) for ii in range(self._parent._k-1) ])
        res=(A.solve_right(B)).transpose()
        return self._parent.free_module()(res.row(0))

    #In HarmonicCocycle
    def evaluate(self,e1):
        r"""
        This function evaluates the cocycle on an edge of the Bruhat-Tits tree.

        INPUT:

        - ``e1`` - a matrix corresponding to an edge of the Bruhat-Tits tree

        OUTPUT

          An element of the coefficient module of the cocycle which describes 
          the value of the cocycle on e1

        EXAMPLES:
        """
        X=self._parent._X
        p=X._p
        u=DoubleCosetReduction(X,e1)
        return (u.sign()*self._F[u.label]).l_act_by(u.igamma(self._parent.embed_quaternion)*(p**(-u.power)))

    #In HarmonicCocycle
    def riemann_sum(self,f,center=1,level=0,E=None):
        r"""
        fixme: what is this? where are you evaluating the form??? Is f the VALUE? This is terrible if so...
        This function evaluates the rigid analytic modular form.
        
        EXAMPLES:

        This example illustrates ...

        ::
        """
        R1=LaurentSeriesRing(f.base_ring(),'r1')
        R1.set_default_prec(self._parent._k-1)
        R2=PolynomialRing(f.base_ring(),'r2')

        if E is None:
            E=self._parent._X._BT.get_balls(center,level)
        else:
            E=self._parent._X._BT.subdivide(E,level)
        value=0
        ii=0
        for e in E:
            ii+=1
            exp=R1([e[1,1],e[1,0]])**(self._parent._k-2)*e.determinant()**(-(self._parent._k-2)/2)*f(R1([e[0,1],e[0,0]])/R1([e[1,1],e[1,0]])).truncate(self._parent._k-1)
            new=self.evaluate(e).l_act_by(e.inverse()).evaluate(exp)
            value+=new
        return value

    def modular_form(self,z=None,level=0):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        return self.derivative(z,level,order=0)

    # In HarmonicCocycle
    def derivative(self,z=None,level=0,order=1):
        r"""
        The derivative of the modular form attached to ``self``.

        EXAMPLES:
        ::
        """
        def F(z):
            R=PolynomialRing(z.parent(),'x,y').fraction_field()
            Rx=PolynomialRing(z.parent(),'x1').fraction_field()
            x1=Rx.gen()
            subst=R.hom([x1,z],codomain=Rx)
            x,y=R.gens()
            center=self._parent._X._BT.find_containing_affinoid(z)
            zbar=z.trace()-z
            f=R(1)/(x-y)
            k=self._parent._k
            V=[f]
            for ii in range(order):
                V=[v.derivative(y) for v in V]+[k/(y-zbar)*v for v in V]
                k+=2
            return sum([self.riemann_sum(subst(v),center,level) for v in V])
        if(z is None):
            return F
        else:
            return F(z)


class HarmonicCocycles(AmbientHeckeModule):
    Element=HarmonicCocycleElement
    r"""
    This object represents a space of Gamma invariant harmonic cocycles valued in
    a cofficient module.

    INPUT:

     - ``X`` - A BTQuotient object

     - ``k`` - integer - The weight.

     - ``prec`` - integer (Default: None). If specified, the precision for the coefficient module

     - ``basis_matrix`` - integer (Default: None)

     - ``base_field`` - (Default: None)

    EXAMPLES:

    ::

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu
    """
    def __init__(self,X,k,prec=None,basis_matrix=None,base_field=None):
        self._k=k
        self._X=X
        self._E=self._X.get_edge_list()
        self._V=self._X.get_vertex_list()

        if prec is None:
            self._prec=None
            if base_field is None:
                try:
                    self._R= X.get_splitting_field()
                except AttributeError:
                    raise ValueError, "It looks like you are not using Magma as backend...and still we don't know how to compute splittings in that case!"
            else:
                pol=X.get_splitting_field().defining_polynomial().factor()[0][0]
                self._R=base_field.extension(pol,pol.variable_name()).absolute_field(name='r')
            self._U=OCVn(self._k-2,self._R)
        else:
            self._prec=prec
            if base_field is None:
                self._R=Qp(self._X._p,prec=prec)
            else:
                self._R=base_field
            self._U=OCVn(self._k-2,self._R,self._k-1)
        if basis_matrix is None:
            if(self._k==2):
                rank=self._X.genus()
            else:
                A=self.basis_matrix()
                rank=A.nrows()
            self.__rank=rank
        else:
            self.__matrix=basis_matrix
            self.__matrix.set_immutable()
            self.__rank=self.__matrix.nrows()
        AmbientHeckeModule.__init__(self, self._R, self.__rank, self._X.prime()*self._X.Nplus()*self._X.Nminus(), weight=self._k)
        self._populate_coercion_lists_()

    def base_extend(self,base_ring):
        r"""
        This function extends the base ring of the coefficient module.
        
        INPUT:

        - ``base_ring`` - a ring that has a coerce map from the current base ring

        EXAMPLES:

        ::
        """
        if not base_ring.has_coerce_map_from(self.base_ring()):
            raise ValueError, "No coercion defined"
        else:
            return self.change_ring(base_ring)

    def change_ring(self, new_base_ring):
        r"""
        This function changes the base ring of the coefficient module.

        INPUT:

        - ``new_base_ring'' - a ring that has a coerce map from the current base ring

        EXAMPLES:
        ::

        """
        if not new_base_ring.has_coerce_map_from(self.base_ring()):
            raise ValueError, "No coercion defined"
        else:
            return self.__class__(self._X,self._k,prec=self._prec,basis_matrix=self.basis_matrix().change_ring(base_ring),base_field=new_base_ring)

    def rank(self):
        r"""
        The rank (dimension) of ``self``.
        
        EXAMPLES:
        ::

        """
        return self.__rank

    def submodule(self,v,check=False):
        r"""
        Return the submodule of ``self`` spanned by ``v``.

        EXAMPLES:
        """
        return HarmonicCocyclesSubmodule(self,v,dual=None,check=check)

    def is_simple(self):
        r"""
        Whether ``self`` is irreducible.

        EXAMPLES:
        ::

        """
        return self.rank()==1

    def _repr_(self):
        r"""
        This returns the representation of self as a string.
        """
        return 'Space of harmonic cocycles of weight %s on %s'%(self._k,self._X)

    def _latex_(self):
        r"""
        A LaTeX representation of ``self``.

        EXAMPLES:
        ::
        """
        s='\\text{Space of harmonic cocycles of weight }'+latex(self._k)+'\\text{ on }'+latex(self._X)
        return s

    def _an_element_(self):
        r"""

        """
        return self.basis()[0]


    def _coerce_map_from_(self, S):
        r"""
        Can coerce from other HarmonicCocycles or from pAutomorphicForms
        """
        if isinstance(S,(HarmonicCocycles,pAutomorphicForms)):
            if(S._k!=self._k):
                return False
            if(S._X!=self._X):
                return False
            return True
        return False

    def __cmp__(self,other):
        r"""

        """
        try:
            res=(self.base_ring()==other.base_ring() and self._X==other._X and self._k==other._k)
            return res
        except:
            return False

    def _element_constructor_(self,x):
        r"""

        """
        #Code how to coherce x into the space
        #Admissible values of x?
        if isinstance(x,HarmonicCocycleElement):
            return HarmonicCocycleElement(self,x)
        elif isinstance(x,pAutomorphicForm):
            tmp=[self._U.element_class(_parent._U,x._F[ii]).l_act_by(self._E[ii].rep) for ii in range(self._nE)]
            return HarmonicCocycleElement(self,tmp,from_values=True)
        else:
            return HarmonicCocycleElement(self,x)


    def free_module(self):
        r"""
        This function returns the underlying free module

        EXAMPLES:
        ::
        """
        try: return self.__free_module
        except AttributeError: pass
        V = self.base_ring()**self.dimension()
        self.__free_module = V
        return V

    def character(self):
        r"""
        Only implemented the trivial character so far.

        EXAMPLES:

        """
        return lambda x:x

    def embed_quaternion(self,g):
        r"""
        Embed the quaternion element ``g`` into the matrix algebra.

        EXAMPLES:
        ::
        """
        return self._X.embed_quaternion(g,exact = self._R.is_exact(), prec = self._prec)

    def basis_matrix(self):
        r"""
        If the coefficient module M is of finite rank then the space of Gamma invariant
        M valued harmonic cocycles can be represented as a subspace of the finite rank
        space of all functions from the finitely many edges in the corresponding 
        BTQuotient into M. This function computes this representation of the space of
        cocycles.

        fixme: if the coefficient module is not finite then this should give an
        error message. At the moment it just breaks.

        OUTPUT:

          A basis matrix describing the cocycles in the spaced of all M valued Gamma
          invariant functions on the tree

        EXAMPLES:

        This example illustrates ...

        ::

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        try: return self.__matrix
        except AttributeError: pass
        nV=len(self._V)
        nE=len(self._E)
        stab_conds=[]
        S=self._X.get_edge_stabs()
        p=self._X._p
        d=self._k-1
        for e in self._E:
            try:
                g=filter(lambda g:g[2],S[e.label])[0]
                C=self._U.l_matrix_representation(self.embed_quaternion(g[0]))
                C-=self._U.l_matrix_representation(Matrix(QQ,2,2,p**g[1]))
                stab_conds.append([e.label,C])
            except IndexError: pass

        n_stab_conds=len(stab_conds)
        self._M=Matrix(self._R,(nV+n_stab_conds)*d,nE*d,0,sparse=True)
        for v in self._V:
            for e in filter(lambda e:e.parity==0,v.leaving_edges):
                C=sum([self._U.l_matrix_representation(self.embed_quaternion(x[0])) for x in e.links],Matrix(self._R,d,d,0))
                self._M.set_block(v.label*d,e.label*d,C)
            for e in filter(lambda e:e.parity==0,v.entering_edges):
                C=sum([self._U.l_matrix_representation(self.embed_quaternion(x[0])) for x in e.opposite.links],Matrix(self._R,d,d,0))
                self._M.set_block(v.label*d,e.opposite.label*d,C)

        for kk in range(n_stab_conds):
            v=stab_conds[kk]
            self._M.set_block((nV+kk)*d,v[0]*d,v[1])

        x1=self._M.right_kernel().matrix()
        K=[c for c in x1.rows()]

        if not self._R.is_exact():
            for ii in range(len(K)):
                s=min([t.valuation() for t in K[ii]])
                for jj in range(len(K[ii])):
                    K[ii][jj]=(p**(-s))*K[ii][jj]

        self.__matrix=Matrix(self._R,len(K),nE*d,K)
        self.__matrix.set_immutable()
        return self.__matrix

    def __apply_atkin_lehner(self,q,f):
        r"""
        This function applies an Atkin-Lehner involution to a harmonic cocycle

        INPUT:

          - ``q`` - an integer dividing the full level p*Nminus*Nplus
          fixme: do we check that q divides?

          - ``f`` - a harmonic cocycle

        OUTPUT:

          The harmonic cocycle obtained by hitting f with the Atkin-Lehner at q

        EXAMPLES:
        ::
        """
        R=self._R
        Data=self._X._get_atkin_lehner_data(q)
        p=self._X._p
        tmp=[self._U.element_class(self._U,zero_matrix(self._R,self._k-1,1),quick=True) for jj in range(len(self._E))]
        d1=Data[1]
        mga=self.embed_quaternion(Data[0])
        for jj in range(len(self._E)):
            t=d1[jj]
            tmp[jj]+=(t.sign()*f._F[t.label]).l_act_by(p**(-t.power)*mga*t.igamma(self.embed_quaternion))
        return HarmonicCocycleElement(self,tmp,from_values=True)

    def test_app_hecke(self,l,f):
        return self.__apply_hecke_operator(l,f)

    def __apply_hecke_operator(self,l,f):
        r"""
        This function applies a Hecke operator to a harmonic cocycle.

        INPUT:

          - ``l`` - an integer

          - ``f`` - a harmonic cocycle

        OUTPUT:

          A harmonic cocycle which is the result of applying the lth Hecke operator
          to f

        EXAMPLES:
        ::

        """
        R=self._R
        HeckeData,alpha=self._X._get_hecke_data(l)
        if(self.level()%l==0):
            factor=QQ(l**(Integer((self._k-2)/2))/(l+1))
        else:
            factor=QQ(l**(Integer((self._k-2)/2)))
        p=self._X._p
        alphamat=self.embed_quaternion(alpha)
        tmp=[self._U.element_class(self._U,zero_matrix(self._R,self._k-1,1),quick=True) for jj in range(len(self._E))]
        for ii in range(len(HeckeData)):
            d1=HeckeData[ii][1]
            mga=self.embed_quaternion(HeckeData[ii][0])*alphamat
            for jj in range(len(self._E)):
                t=d1[jj]
                tmp[jj]+=(t.sign()*f._F[t.label]).l_act_by(p**(-t.power)*mga*t.igamma(self.embed_quaternion))
        return HarmonicCocycleElement(self,[factor*x for x in tmp],from_values=True)

    def _compute_atkin_lehner_matrix(self,d):
        r"""
        When the underlying coefficient module is finite, this function computes the 
        matrix of an Atkin-Lehner involution in the basis provided by the function
        basis_matrix

        INPUT:

          - ``d`` - an integer dividing p*Nminus*Nplus

        OUTPUT:

          The matrix of the AL-involution at d in the basis given by self.basis_matrix

        EXAMPLES:
        ::

        """
        res=self.__compute_operator_matrix(lambda f:self.__apply_atkin_lehner(d,f))
        return res

    def _compute_hecke_matrix_prime(self,l):
        r"""
        When the underlying coefficient module is finite, this function computes the 
        matrix of a (prime) Hecke operator in the basis provided by the function
        basis_matrix

        INPUT:

          - ``l`` - an integer prime

        OUTPUT:

          The matrix of T_l acting on the cocycles in the basis given by 
          self.basis_matrix

        EXAMPLES:
        ::

        """
        res=self.__compute_operator_matrix(lambda f:self.__apply_hecke_operator(l,f))
        return res

    def __compute_operator_matrix(self,T):
        r"""
        Compute the matrix of the operator ``T``.

        EXAMPLES:
        ::

        """
        R=self._R
        A=self.basis_matrix().transpose()
        basis=self.basis()
        B=zero_matrix(R,len(self._E)*(self._k-1),self.dimension())
        for rr in range(len(basis)):
            g=T(basis[rr])
            B.set_block(0,rr,Matrix(R,len(self._E)*(self._k-1),1,[g._F[e]._val[ii,0]  for e in range(len(self._E)) for ii in range(self._k-1) ]))

        try:
            res=(A.solve_right(B)).transpose()
            res.set_immutable()
            return res
        except ValueError:
            print A
            print B
            raise ValueError

class HarmonicCocyclesSubmodule(sage.modular.hecke.submodule.HeckeSubmodule,HarmonicCocycles):
    r"""
    fixme:you can explain this
    This function returns the point `(x^5,y)`.

    INPUT:

     - ``x`` - integer (default: 1) the description of the
       argument x goes here.  If it contains multiple lines, all
       the lines after the first need to be indented.

     - ``y`` - integer (default: 2) the ...

    EXAMPLES:

    AUTHOR:

    - Marc Masdeu (2012-02-20)
    """
    def __init__(self, ambient_module, submodule, dual=None, check=False):
        """
            ambient_module -- HarmonicCocycles
            submodule -- a submodule of the ambient space.
            dual_module -- (default: None) ignored
            check -- (default: False) whether to check that the
                     submodule is Hecke equivariant
        """
        A = ambient_module
        sage.modular.hecke.submodule.HeckeSubmodule.__init__(self, A, submodule, check=check)
        self.__rank=submodule.dimension()
        HarmonicCocycles.__init__(self,X=A._X,k=A._k,prec=A._prec,basis_matrix=submodule.basis_matrix()*A.basis_matrix())

    def rank(self):
        r"""
        Returns the rank (dimension) of the submodule.

        OUTPUT:

        integer - The rank of ``self``.

        """
        return self.__rank

    def _repr_(self):
        r"""
        Returns the representation of self as a string.
        """
        return "Subspace of %s of dimension %s"%(self.ambient(),self.dimension())


class pAutomorphicForm(ModuleElement):
    r"""
    This class is a rudimentary implementation of a class for a p-adic automorphic
    form on a definite quaternion algebra over Q. These are required in order to
    compute moments of measures associated to harmonic cocycles on the BT-tree
    using the overconvergent modules of Darmon-Pollack and Matt Greenberg. See
    Greenberg's thesis for more details.

    INPUT:
    fixme: describe these
     - ``vec`` - 

     - ``quick`` - boolean (default: False) 

    EXAMPLES:

    This example illustrates ...

    ::

        
        
        

    ...

    REFERENCES:

    Matthew Greenberg's thesis (available on his webpage as of 02/12).

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu


    """
    def __init__(self,_parent,vec,quick=False):
        ModuleElement.__init__(self,_parent)
        self._parent=_parent
        self._nE=2*len(_parent._E) # We record the values at the opposite edges
        self._cached_values=dict()
        self._R=Qp(_parent._X._p,prec=_parent._prec)
        if(quick):
                self._F=[v for v in vec]
        else:
            if(isinstance(vec,pAutomorphicForm)):
                self._F=[self._parent._U(vec._F[ii]) for ii in range(self._nE)]
                self._make_invariant()

            elif(isinstance(vec,HarmonicCocycleElement)):
                assert(_parent._U.weight()==vec._wt-2)
                self._F=[]
                assert(2*len(vec._F)==self._nE)
                assert(isinstance(_parent._U,OCVn))
                E=self._parent._E


                MMM=vec._parent._U.element_class
                tmp=[]
                for ii in range(len(vec._F)):
                    newtmp=MMM(vec._parent._U,vec._F[ii]).l_act_by(E[ii].rep.inverse())
                    tmp.append(newtmp)
                    self._F.append(_parent._U(newtmp))
                A=Matrix(QQ,2,2,[0,-1/_parent._X._p,-1,0])
                for ii in range(len(vec._F)):
                    self._F.append(_parent._U(-1*tmp[ii].r_act_by(A)))
                self._make_invariant()

            elif(isinstance(vec,list) and len(vec)==self._nE):
                try:
                    self._F=[self._parent._U(v) for v in vec]
                except:
                    try:
                        veczp=_parent._U._R(vec)
                        self._parent=_parent
                        self._F=[self._parent._U(veczp) for ii in range(self._nE)]
                    except:
                        print vec
                        assert(0)
            else:
                try:
                    veczp=_parent._U._R(vec)
                    self._parent=_parent
                    self._F=[self._parent._U(veczp) for ii in range(self._nE)]
                except:
                    raise ValueError,"Cannot initialize a p-adic automorphic form with the given input="+str(vec)

    def _make_invariant(self):
        r"""
        fixme: have to double check what this does
        This function takes
        
        EXAMPLES:

        This example illustrates ...

        ::

            
            
            
        """
        S=self._parent._X.get_edge_stabs()
        E=self._parent._E
        M=[e.rep for e in E]+[e.opposite.rep for e in E]
        lS=len(S)
        assert(2*len(S)==self._nE)
        MMM=self._parent._U.element_class
        for ii in range(self._nE):
            Si=S[ii%lS]
            if(any([v[2] for v in Si])):
                x=MMM(self._F[ii].parent(),self._F[ii]._val,quick=True)
                self._F[ii]=MMM(self._F[ii].parent(),0)
                s=QQ(0)
                m=M[ii]
                for v in Si:
                    s+=1
                    self._F[ii]+=x.r_act_by(m.adjoint()*self._parent.embed_quaternion(v[0])*m)
                self._F[ii]=(1/s)*self._F[ii]

    def precision(self):
        r"""
        The precision of ``self``, which is the minimum among the precision of the values on a fundamental domain.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return min(x.precision() for x in self._F)

    def _add_(self,g):
        r"""
        This function adds two p-adic automorphic forms.
        
        INPUT:

          - ``g`` - a p-adic automorphic form

        OUTPUT:

          the result of adding g to self

        EXAMPLES:

        This example illustrates ...

        ::
        """
        #Should ensure that self and g are of the same weight and on the same curve
        vec=[self._F[e]+g._F[e] for e in range(self._nE)]
        return pAutomorphicForm(self._parent,vec,quick=True)
    
    def _sub_(self,g):
        r"""
        This function subtracts a p-adic automorphic form from another.
        
        INPUT:
      
          - ``g`` - a p-adic automorphic form

        OUTPUT:

          the result of subtracting g from self

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #Should ensure that self and g are of the same weight and on the same curve
        vec=[self._F[e]-g._F[e] for e in range(self._nE)]
        return pAutomorphicForm(self._parent,vec,quick=True)

    def _getitem_(self,e1):
        r"""
        Evaluates a p-adic automorphic form on a matrix in GL2(Qp).

        INPUT:

          - ``e1`` - a matrix in GL2(Qp)

        OUTPUT:

          the value of self evaluated on e1

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return self.evaluate(e1)

    def evaluate(self,e1):
        r"""
        Evaluates a p-adic automorphic form on a matrix in GL2(Qp).
        
        INPUT:
     
          - ``e1`` - a matrix in GL2(Qp)

        OUTPUT:

          the value of self evaluated on e1

        EXAMPLES:

        This example illustrates ...

        ::

        """
        X=self._parent._X
        u=DoubleCosetReduction(X,e1)
        return (self._F[u.label+u.parity*self._nE/2]).r_act_by((u.t())*X._p**(u.power))

    def _rmul_(self,a):
        r"""
        Multiplies the automorphic form by a scalar.
        
        EXAMPLES:

        This example illustrates ...

        ::

        """
        #Should ensure that 'a' is a scalar
        return pAutomorphicForm(self._parent,[a*self._F[e] for e in range(self._nE)],quick=True)

    def left_act_by(self,alpha):
        r"""
        Act by ``alpha`` on the left.
        
        EXAMPLES:

        This example illustrates ...

        ::

        """
        Y=self._parent._X
        E=self._parent._E
        p=Y._p
        Tf=[]
        for e in E:
            u=DoubleCosetReduction(Y,alpha*e.rep)
            r=u.t()*p**(-(u.power))
            if(u.parity==0):
                tmp=self._F[u.label].r_act_by(r)
            else:
                tmp=self._F[u.label+len(E)].r_act_by(r)
            Tf.append(tmp)
        for e in E:
            u=DoubleCosetReduction(Y,alpha*e.opposite.rep)
            r=u.t()*gg*p**(-(u.power))
            if(u.parity==0):
                tmp=self._F[u.label].r_act_by(r)
            else:
                tmp=self._F[u.label+len(E)].r_act_by(r)
            Tf.append(tmp)
        return pAutomorphicForm(self._parent,Tf,quick=True)


    def _repr_(self):
        r"""
        This returns the representation of self as a string.
        
        EXAMPLES:

        This example illustrates ...

        ::

        """
        tmp='p-adic automorphic form on '+str(self._parent)+':\n'
        tmp+='   e   |   c(e)'
        tmp+='\n'
        for e in range(Integer(self._nE/2)):
            tmp+='  '+str(e)+' | '+str(self._F[e])+'\n'
        return tmp

    def valuation(self):
        r"""
        The valuation of ``self``, defined as the minimum of the valuations of the values that it takes on a set of edge representatives.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        if not self.__nonzero__():
            return oo
        else:
            return(min([self._F[e].valuation() for e in range(self._nE)]))

    def __nonzero__(self):
        r"""
        """
        return(any([self._F[e].__nonzero__() for e in range(self._nE)]))

    def improve(self, verbose = True):
        r"""
        Repeatedly applies the `U_p` operator to a p-adic automorphic form. This
        is used to compute moments of a measure associated to a rigid modular form in the
        following way: lift a rigid modular form to an ``overconvergent'' `p`-adic automorphic
        form in any way, and then repeatedly apply `U_p` to project to the ordinary part.
        The resulting form encodes the moments of the measure of the original rigid modular 
        form (assuming it is ordinary). 


        EXAMPLES:


        REFERENCES:

        For details see Matthew Greenberg's thesis (available on his webpage as of 02/12).
        Alternatively check out Darmon-Pollack for the analogous algorithm in the case of
        modular symbols.

        AUTHORS:


        - Cameron Franc (2012-02-20)
        - Marc Masdeu

        """
        MMM=self._parent
        if verbose == True:
            sys.stdout.flush()
        h2=MMM._apply_Up_operator(self,True)
        if verbose == True:
            sys.stdout.write("#")
            sys.stdout.flush()
        ii=0
        current_val=0
        old_val=-Infinity
        init_val=self.valuation()
        while(current_val>old_val):
            old_val=current_val
            ii+=1
            self._F=[self._parent._U(c) for c in h2._F]
            h2=MMM._apply_Up_operator(self,True)
            current_val=(h2-self).valuation()-init_val
            if current_val is Infinity:
                break
            if verbose == True:
                sys.stdout.write("#")
                sys.stdout.flush()
        if verbose == True:
            print ''
        self._F=[self._parent._U(c) for c in h2._F]

    def integrate(self,f,center=1,level=0,method='moments'):
        r"""
        Calculate
        .. MATH::

            \int_{\PP^1(\QQ_p)} f(x)d\mu(x)

        were `\mu` is the measure associated to ``self``.

        INPUT:

         - ``f`` - An analytic function.

         - ``center`` - 2x2 matrix over Qp (default: 1)

         - ``level`` - integer (default: 0)

         - ``method`` - string (default: 'moments'). Which method of integration to use. Either 'moments' or 'riemann_sum'.


        EXAMPLES:

        ::


        AUTHORS:

        - Marc Masdeu (2012-02-20)
        - Cameron Franc (2012-02-20)

        """
        E=self._parent._X._BT.get_balls(center,level)
        R1=LaurentSeriesRing(f.base_ring(),'r1')
        R2=PolynomialRing(f.base_ring(),'r2')
        value=0
        ii=0
        if(method=='riemann_sum'):
            R1.set_default_prec(self._parent._U.weight()+1)
            for e in E:
                ii+=1
                #print ii,"/",len(E)
                exp=((R1([e[1,1],e[1,0]]))**(self._parent._U.weight())*e.determinant()**(-(self._parent._U.weight())/2))*f(R1([e[0,1],e[0,0]])/R1([e[1,1],e[1,0]]))
                #exp=R2([tmp[jj] for jj in range(self._parent._k-1)])
                new=self.evaluate(e).evaluate(exp.truncate(self._parent._U.weight()+1))
                value+=new
        elif(method=='moments'):
            R1.set_default_prec(self._parent._U.dimension())
            for e in E:
                ii+=1
                #print ii,"/",len(E)
                tmp=((R2([e[1,1],e[1,0]]))**(self._parent._U.weight())*e.determinant()**(-(self._parent._U.weight())/2))*f(R2([e[0,1],e[0,0]])/R2([e[1,1],e[1,0]]))
                exp=R1(tmp.numerator())/R1(tmp.denominator())
                new=self.evaluate(e).evaluate(exp)
                value+=new
        else:
            print 'The available methods are either "moments" or "riemann_sum". The latter is only provided for consistency check, and should never be used.'
            return False
        return value

    def modular_form(self,z=None,level=0,method='moments'):
        r"""
        Returns the modular form corresponding to ``self``.
        
        INPUT:

          - ``z`` - (default: None). If specified, returns the value of the form at the point ``zz`` in the `p`-adic upper half plane.

          - ``level`` - integer (default: 0). If ``method`` is 'riemann_sum', will use a covering of `\PP^1(\QQ_p)` with balls of size `p^-\mbox{level]`.
        
          - ``method`` - string (default: ``moments``). It must be either ``moments`` or ``riemann_sum``.

        OUTPUT:

        A function from the `p`-adic upper half plane to `\CC_p`. If an argument ``z`` was passed, returns instead the value at that point.

        """
        return self.derivative(z,level,method,order=0)

    def derivative(self,z=None,level=0,method='moments',order=1):
        r"""
        Returns the derivative of the modular form corresponding to
        ``self``.

        INPUT:

          - ``z`` - (Default: None). If specified, evaluates the derivative
              at the point ``z`` in the `p`-adic upper half plane.

          - ``level`` - integer (default: 0). If ``method`` is 'riemann_sum', will use a covering of `\PP^1(\QQ_p)` with balls of size `p^-\mbox{level]`.
        
          - ``method`` - string (default: ``moments``). It must be either ``moments`` or ``riemann_sum``.

          - ``order`` - integer (Default: 1). The order of the derivative to be computed.

        OUTPUT:

        A function from the `p`-adic upper half plane to `\CC_p`. If an argument ``z`` was passed, returns instead the value of the derivative at that point.

        EXAMPLES:

        """
        def F(z,level=0,method='moments'):
            R=PolynomialRing(z.parent(),'x,y').fraction_field()
            Rx=PolynomialRing(z.parent(),'x1').fraction_field()
            x1=Rx.gen()
            subst=R.hom([x1,z],codomain=Rx)
            x,y=R.gens()
            center=self._parent._X._BT.find_containing_affinoid(z)
            zbar=z.trace()-z
            f=R(1)/(x-y)
            k=self._parent._n+2
            V=[f]
            for ii in range(order):
                V=[v.derivative(y) for v in V]+[k/(y-zbar)*v for v in V]
                k+=2
            return sum([self.integrate(subst(v),center,level,method) for v in V])
        if z is None:
            return F

        try:
            try: tmp=self._cached_values[z,level,method,order,self.precision()]
            except KeyError:
                tmp=F(z,level,method)
                self._cached_values[z,level,method,order,self.precision()]=tmp
            return tmp
        except TypeError:
            return F(z,level,method)


    # So far we can't break it into two integrals because of the pole at infinity.
    def coleman(self,t1,t2,E=None,method='moments',mult=False,delta=-1,level=0):
        r"""
        If ``self`` is a `p`-adic automorphic form that corresponds to a rigid modular form,
        then this computes the coleman integral of this form between two points on
        the boundary `\PP^1(\QQ_p)` of the `p`-adic upper half plane.

        INPUT:

         - ``t1``, ``t2`` - elements of `\PP^1(\QQ_p)` (the endpoints of integration)

         - ``E`` - (Default: None). If specified, will not compute the covering adapted
            to ``t1`` and ``t2`` and instead use the given one. In that case, ``E``
            should be a list of matrices corresponding to edges describing the open
            balls to be considered.

         - ``method`` - string (Default: 'moments'). Tells which algorithm to use 
           (alternative is 'riemann_sum', which is unsuitable for computations 
           requiring high precision)

         - ``mult`` - boolean (Default: False). Whether to use the multiplicative version.

         - ``delta`` - integer (Default: -1)

         - ``level`` - integer (Default: 0)

        OUTPUT:

          The result of the coleman integral

        EXAMPLES:

        This example illustrates ...

        ::

        TESTS:

        ::
            sage: p = 7
            sage: lev = 2
            sage: prec = 20
            sage: X = BTQuotient(p,lev)
            sage: k = 2
            sage: M = HarmonicCocycles(X,k,prec)
            sage: B = M.basis()
            sage: f = 3*B[0]
            sage: MM = pAutomorphicForms(X,k,prec,overconvergent=True)
            sage: D = -11
            sage: X.is_admissible(D)
            True
            sage: K.<a> = QuadraticField(D)
            sage: CM = X.get_CM_points(D,prec=prec)
            sage: Kp = CM[0].parent()
            sage: P=CM[0]
            sage: Q=P.trace()-P # The conjugate
            sage: F=MM.lift(f, verbose = False) # long time
            sage: J0=F.coleman(P,Q,mult=True) # long time
            sage: E=EllipticCurve([1,0,1,4,-6])
            sage: T=E.tate_curve(p)
            sage: xx,yy=getcoords(T,J0,prec) # long time
            sage: P = E.base_extend(K).lift_x(algdep(xx,1).roots(QQ)[0][0]); P # long time
            (7/11 : 58/121*a - 9/11 : 1)

            sage: p = 13
            sage: lev = 2
            sage: prec = 20
            sage: Y = BTQuotient(p,lev)
            sage: k = 2
            sage: M=HarmonicCocycles(Y,k,prec)
            sage: B=M.basis()

            sage: f = B[1]
            sage: g = -4*B[0]+3*B[1]
            sage: MM = pAutomorphicForms(Y,k,prec,overconvergent=True)
            sage: D = -11
            sage: Y.is_admissible(D)
            True
            sage: K.<a> = QuadraticField(D)
            sage: CM = Y.get_CM_points(D,prec=prec)
            sage: Kp = parent(CM[0])
            sage: P = CM[0]
            sage: Q = P.trace()-P
            sage: F = MM.lift(f, verbose = False) # long time
            sage: J11 = F.coleman(P,Q,mult = True) # long time
            sage: E = EllipticCurve('26a2')
            sage: T = E.tate_curve(p)
            sage: xx,yy = getcoords(T,J11,prec) # long time
            sage: HP = E.base_extend(K).lift_x(algdep(xx,1).roots(QQ)[0][0]); HP # long time
            (-137/11 : 2/121*a + 63/11 : 1)

        AUTHORS:

        - Cameron Franc (2012-02-20)
        - Marc Masdeu (2012-02-20)
        """
        if(mult and delta>=0):
            raise NotImplementedError, "Need to figure out how to implement the multiplicative part."
        p=self._parent._X._p
        K=t1.parent()
        R=PolynomialRing(K,'x')
        x=R.gen()
        R1=LaurentSeriesRing(K,'r1')
        r1=R1.gen()
        if(E is None):
            E=self._parent._X._BT.find_covering(t1,t2)
            # print 'Got %s open balls.'%len(E)
        value=0
        ii=0
        value_exp=K(1)
        if(method=='riemann_sum'):
            R1.set_default_prec(self._parent._U.weight()+1)
            for e in E:
                ii+=1
                b=e[0,1]
                d=e[1,1]
                y=(b-d*t1)/(b-d*t2)
                poly=R1(y.log()) #R1(our_log(y))
                c_e=self.evaluate(e)
                new=c_e.evaluate(poly)
                value+=new
                if mult:
                    value_exp *= K.teichmuller(y)**Integer(c_e[0].rational_reconstruction())

        elif(method=='moments'):
            R1.set_default_prec(self._parent._U.dimension())
            for e in E:
                ii+=1
                f=(x-t1)/(x-t2)
                a,b,c,d=e.list()
                y0=f(R1([b,a])/R1([d,c])) #f( (ax+b)/(cx+d) )
                y0=p**(-y0(0).valuation())*y0
                mu=K.teichmuller(y0(0))
                y=y0/mu-1
                poly=R1(0)
                ypow=y
                for jj in range(1,R1.default_prec()+10):
                    poly+=(-1)**(jj+1)*ypow/jj
                    ypow*=y
                if(delta>=0):
                    poly*=((r1-t1)**delta*(r1-t2)**(self._parent._n-delta))
                c_e=self.evaluate(e)
                new=c_e.evaluate(poly)
                value+=new
                if mult:
                    value_exp *= K.teichmuller((b-d*t1)/(b-d*t2))**Integer(c_e[0].rational_reconstruction())

        else:
            print 'The available methods are either "moments" or "riemann_sum". The latter is only provided for consistency check, and should not be used in practice.'
            return False
        if mult:
            return value_exp * value.exp()
        return value


class pAutomorphicForms(Module):
    Element=pAutomorphicForm
    r"""
    The module of (quaternionic) `p`-adic automorphic forms.

    EXAMPLES:

    ::


    AUTHORS:


    - Cameron Franc (2012-02-20)
    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,X,U,prec=None,t=None,R=None,overconvergent=False):
        if(R is None):
            if(not isinstance(U,Integer)):
                self._R=U.base_ring()
            else:
                if(prec is None):
                    prec=100
                self._R=Qp(X._p,prec)
        else:
            self._R=R
        #U is a CoefficientModuleSpace
        if(isinstance(U,Integer)):
            if(t is None):
                if(overconvergent):
                    t=prec-U+1
                else:
                    t=0
            self._U=OCVn(U-2,self._R,U-1+t)
        else:
            self._U=U
        self._X=X
        self._V=self._X.get_vertex_list()
        self._E=self._X.get_edge_list()
        self._prec=self._R.precision_cap()
        self._n=self._U.weight()
        Module.__init__(self,base=self._R)
        self._populate_coercion_lists_()

    def _repr_(self):
        r"""
        This returns the representation of self as a string.
        
        EXAMPLES:

        This example illustrates ...

        ::
        """
        s='Space of automorphic forms on '+str(self._X)+' with values in '+str(self._U)
        return s

    def _coerce_map_from_(self, S):
        r"""
        Can coerce from other HarmonicCocycles or from pAutomorphicForms
        """
        if(isinstance(S,HarmonicCocycles)):
            if(S._k-2!=self._n):
                return False
            if(S._X!=self._X):
                return False
            return True
        if(isinstance(S,pAutomorphicForms)):
            if(S._n!=self._n):
                return False
            if(S._X!=self._X):
                return False
            return True
        return False

    def _element_constructor_(self,x):
        r"""
        """
        #Code how to coherce x into the space
        #Admissible values of x?
        if isinstance(x,(HarmonicCocycleElement,pAutomorphicForm)):
            return pAutomorphicForm(self,x)

    def _an_element_(self):
        r"""
        Returns an element of the module.
        """
        return pAutomorphicForm(self,1)

    def embed_quaternion(self,g):
        r"""
        Returns the image of ``g`` under the embedding
        of the quaternion algebra into 2x2 matrices.
        """
        return self._X.embed_quaternion(g,prec = self._prec)

    def lift(self,f, verbose = True):
        r"""
        Lifts the harmonic cocycle ``f`` to an
        overconvegent automorphic form, thus computing
        all the moments.
        """
        F=self.element_class(self,f)
        F.improve(verbose = verbose)
        return F

    def _apply_Up_operator(self,f,scale=False):
        r"""
        Apply the Up operator to ``f``.
        """
        HeckeData=self._X._get_Up_data()
        if(scale==False):
            factor=(self._X._p)**(self._U.weight()/2)
        else:
            factor=1
        Tf=[]
        for jj in range(2*len(self._E)):
            tmp=self._U(0)
            for d in HeckeData:
                gg=d[0] # acter
                u=d[1][jj] # edge_list[jj]
                r=self._X._p**(-(u.power))*(u.t()*gg)
                tmp+=f._F[u.label+u.parity*len(self._E)].r_act_by(r)
            Tf.append(factor*tmp)
        return pAutomorphicForm(self,Tf,quick=True)

