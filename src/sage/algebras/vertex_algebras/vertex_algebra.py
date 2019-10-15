r"""
Vertex algebras
AUTHORS

- Reimundo Heluani (08-09-2019): Initial implementation
"""


#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.combinat.partition_tuple import PartitionTuples_level
from sage.categories.vertex_algebras import VertexAlgebras
from sage.categories.lie_conformal_algebras import LieConformalAlgebras, LiftMorphism
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partition
from sage.combinat.partition_tuple import PartitionTuples
from sage.functions.other import binomial, factorial
from sage.misc.misc_c import prod
from sage.rings.infinity import Infinity
from sage.sets.family import Family
from sage.structure.element import parent
from sage.rings.integer import Integer
from sage.misc.cachefunc import cached_method


class VertexAlgebra(Parent, UniqueRepresentation):
    @staticmethod
    def __classcall_private__(cls, R=None, arg0 = None, **kwds):
        if R is None: 
            raise ValueError("First argument must be a ring!")
       
        try: 
            if R.has_coerce_map_from(arg0.base_ring()) and \
                arg0 in LieConformalAlgebras(arg0.base_ring()):
                return UniversalEnvelopingVertexAlgebra(R, arg0, **kwds)
        except AttributeError: 
            pass

        print "Nothing to construct"

    def __init__(self, R, **kwds):
        category = kwds.get('category', None)
        category = VertexAlgebras(R).or_subcategory(category)
        super(VertexAlgebra, self).__init__(base=R, names = 
            kwds.get('names', None), category = category)

    def base_ring(self):
        return self.category().base_ring()


class UniversalEnvelopingVertexAlgebra(VertexAlgebra):

    def __init__(self, R, L, **kwds):

        if L not in LieConformalAlgebras(R).WithBasis().FinitelyGenerated():
            raise ValueError ( "L needs to be a finitely generated " \
                "Lie conformal algebra with basis, got {}".format(L) )

        category=kwds.get('category', None)
        
        category = VertexAlgebras(R).FinitelyGenerated().WithBasis().\
           or_subcategory(category)

        if L in LieConformalAlgebras(R).HGraded():
            category = category.HGraded()

        kwds['category'] = category

        super(UniversalEnvelopingVertexAlgebra, self).__init__(R, **kwds)

        self._lca = L
        cp = kwds.get('central_parameters', {})
        if cp == {}:
            cp = { i:0  for i in L.central_elements() }
        if set(cp.keys()) != set(L.central_elements()):
            raise ValueError ("central_parameters must be parametrized by "\
                              "central elements")

        self._central_parameters = Family(cp)
        #need to call directly this because of 1 generator. 
        #Also:self._module is needed for self.ngens()
        _basis = PartitionTuples_level(L.ngens()-len(L.central_elements()))
        self._module = CombinatorialFreeModule(self.base_ring(), _basis)
        self.register_lift()
    
    def _repr_(self):
        return "The universal enveloping vertex algebra of {}".format(
            self._lca)

    def register_lift(self):
        from sage.categories.homset import Hom
        self._lca.lift = LiftMorphism(Hom(self._lca, self, category = 
                        LieConformalAlgebras(self._lca.base_ring())))
        try: 
            self._lca.lift.register_as_coercion()
        except AssertionError:
            #we already constructed this morphisms and its fine
            pass



    def basis(self,n=None):
        if n == None:
            return Family( PartitionTuples_level(level = self.ngens()), 
                    lambda i : self(i), lazy = True)
        else:
            return Family ( (self(self.module()._from_dict(
                b.monomial_coefficients())) for b in 
                self.get_graded_part(n).basis() ) )

    def gens(self):
        return tuple(self.gen(i) for i in range(self.ngens()))

    def ngens(self):
        """
        The number of generators of this vertex algebra

        We can't call directly gens cause we need this to construc the lift
        """
        return self._lca.ngens() - len(self._lca.central_elements())

    def gen(self,i):
        l = [[]]*self.ngens()
        l[i] = [1]
        return self(l)

    def central_elements(self):
        return tuple ( self(i) for i in self._lca.central_elements())

    def central_parameters(self):
        return self._central_parameters

    def get_degree(self,n):
        """
        return the degree n filtered part of `self`. This is a 
        `CombinatorialFreeModule` with the same dictionary keys as 
        `self.module()`. 
        """
        #TODO: deal with the vacuum in a better way
        basis = [PartitionTuples_level(self.ngens())([[],]*self.ngens()),]
        basis += [PartitionTuples_level(self.ngens())(p) for m in range(1 ,n+1) 
            for p in PartitionTuples(self.ngens(),m) if 
            self(PartitionTuples_level(self.ngens())(p)).weight() <= n ]
        return CombinatorialFreeModule(self.base_ring(), basis)

    def get_graded_part(self,n):
        """
        return the degree n filtered part of `self`. This is a 
        `CombinatorialFreeModule` with the same dictionary keys as 
        `self.module()`. 
        """
        #TODO: deal with the vacuum in a better way
        if n == 0:
            return CombinatorialFreeModule( self.base_ring(),
                [PartitionTuples_level(self.ngens())([[]]*self.ngens()),] )
        else:
            basis = [PartitionTuples_level(self.ngens())(p) for m in range(1 ,n+1) 
            for p in PartitionTuples(self.ngens(),m) if 
            self(PartitionTuples_level(self.ngens())(p)).weight() == n ]
        return CombinatorialFreeModule(self.base_ring(), basis)

    def dimension(self,n):
        return self.get_graded_part(n).dimension()

    def module(self):
        return self._module

    @cached_method
    def vacuum(self):
        vac = [Partition([]),]*self.ngens()
        return self.element_class(self, self.module()(vac))

    @cached_method
    def zero(self):
        return self.element_class(self, self.module().zero())

    def _element_constructor_(self,x):
        if x == self.base_ring().zero():
            return self.zero()
        try:
            v = self._module(x)
        except TypeError:
            raise TypeError("do not know how to convert {0} into an element "\
                            "of {1}".format(x,self))
        return self.element_class(self, v)

    def li_filtration(self,n,k=None):
        A = self.get_graded_part(n)
        ret = {}
        for m in range(n+2) if k==None else range(k,k+1):
            basis = [b for b in A.basis() if self._from_dict(
                b.monomial_coefficients()).li_filtration_degree()>=m]
            ret[m] = A.submodule(basis)
        return ret if k==None else ret[k]
         
    from sage.structure.element_wrapper import ElementWrapper
    class Element(ElementWrapper):
        def _repr_(self):
            p = self.parent()
            if self == p.zero():
                return "0";
            coeff = self.value.monomial_coefficients().items()
            ret = ""
            for i in range(len(coeff)):
                if i > 0  and coeff[i][1] > 0:
                    ret += "+"
                if coeff[i][1] < 0:
                    ret += "-"

                if abs(coeff[i][1]) != 1:
                    ret += "{}*".format(abs(coeff[i][1]))
                for idx,j in enumerate(coeff[i][0]):
                    for k in j.to_list():
                        ret += "{}_".format(p._lca._repr_generator(p._lca.gen(idx)))
                        if p.is_graded():
                            ret+= "{}".format(1 -k-p._lca.gen(idx).degree())
                        else:
                            ret+= "({})".format(-k)
                ret+="|0>"
            return ret

        def _add_(self,right):
            return type(self)(self.parent(), self.value + right.value)

        def __nonzero__(self):
            return bool(self.value)

        def _sub_(self, right):
            return type(self)(self.parent(), self.value - right.value)

        def _bracket_(self, other):
            p = self.parent()
            ret = {}
            svmc = other.value.monomial_coefficients()
            for k in svmc.keys():
                c = svmc[k]
                i = next((j for j,x in enumerate(k) if x), None)
                if i is None:
                    continue
                l = [Partition([]),]*p.ngens()
                l[i] = Partition([k[i].get_part(0)])
                if k == l:
                    svmc2 = self.value.monomial_coefficients()
                    for k2 in svmc2.keys():
                        c2 = svmc2[k2]
                        i2 = next((j for j,x in enumerate(k2) if x), None)
                        if i2 is None:
                            continue
                        l2 = [Partition([]),]*p.ngens()
                        l2[i2] = Partition([k2[i2].get_part(0)])
                        if k2 == l2:
                            br = p._lca.gen(i2).T(l2[i2].size()-1) \
                                .bracket(p._lca.gen(i).T(l[i].size()-1))
                            for j in br.keys():
                                rec = ret.get(j,p.zero())
                                ret[j] = rec + factorial(l[i].size()-1)**(-1)*\
                                factorial(l2[i2].size()-1)**(-1)*c*c2*p(br[j])
                        else:
                            #skew-symmetry
                            br = p(k)._bracket_(p(k2))
                            for l in range(max(br.keys())+1):
                                rec = ret.get(l, p.zero())
                                ret[l] = rec + sum( c*c2*(-1)**(j+1)*br[j].T(j-l) 
                                            for j in br.keys() if j >= l )
                else:
                    #non-commutative Wick Formula
                    sf = k.to_list()
                    sf[i] = Partition(Partition(sf[i]).to_list()[1:])
                    br = self._bracket_(p(l))
                    for j in br.keys():
                        rec = ret.get(j, p.zero())
                        ret[j] = rec + c*br[j]._mul_(p(sf))
                        #integral term
                        br2 = br[j]._bracket_(p(sf))
                        for m in br2.keys():
                            rec = ret.get(j+m+1 ,p.zero())
                            ret[j+m+1] = rec + c*binomial(j+m+1 ,j)*br2[m]
                    br = self._bracket_(p(sf))
                    for j in br.keys():
                        rec = ret.get(j, p.zero())
                        ret[j] = rec + c*p(l)._mul_(br[j])
            return { k:ret[k] for k in ret.keys() if ret[k] } 

        def _mul_(self,right):
            """
            returns the normally ordered product :`self``right`:
            This is where the magic happens
            """
            p = self.parent()
            ret = p.zero()
            svmc = self.value.monomial_coefficients()
            for k in svmc.keys():
                c = svmc[k]
                i = next((j for j,x in enumerate(k) if x), None)
                if i is None:
                    ret += c*right
                else:
                    l = [Partition([]),]*p.ngens()
                    l[i] = Partition([k[i].get_part(0)])
                    if k == l:
                        ni = k[i].get_part(0)
                        rvmc = right.value.monomial_coefficients()
                        for k2 in rvmc.keys():
                            c2 = rvmc[k2]
                            i2 = next((j for j,x in enumerate(k2) if x), None)
                            if i2 is None:
                                ret += c*c2*p(l)
                            elif  i2 < i or (i2==i and k2[i2].get_part(0)>ni):
                                ni2 = k2[i2].get_part(0)
                                l2 = [Partition([]),]*p.ngens()
                                l2[i2] = Partition([ni2])
                                sf = k2.to_list()
                                sf[i2] = Partition(Partition(sf[i2]).to_list()[1:])
                                rec = p(l)._mul_(p(sf))
                                ret += c*c2*p(l2)._mul_(rec)
                                #Here's the only place where we need the
                                #structure constants. We are using that
                                #T^(n)a_{(-1)} = a_{(-1-n)} and the commutator 
                                #formula. Given the implementation it's faster
                                #to compute the bracket in the other direction
                                br = p._lca.gen(i2).bracket(p._lca.gen(i))
                                if br:
                                    ret -= c*c2*sum( binomial(-ni2,j)*\
                                    factorial(ni+ni2+j-1)**(-1)*p(br[j].\
                                    T(ni+ni2+j-1))._mul_(p(sf)) for j in br.keys() )
                            else:
                                sf = k2.to_list()
                                l2 = Partition(k2[i]).to_list()
                                l2.insert(0,ni)
                                sf[i] = Partition(l2)
                                ret += c*c2*p(sf)
                    else:
                        #quasi-associativity
                        l = p(l)
                        sf = k.to_list()
                        sf[i] = Partition(Partition(sf[i]).to_list()[1:])
                        sf = p(sf)
                        rec = sf._mul_(right)
                        ret+= c*l._mul_(rec)
                        sfbr = sf.bracket(right)
                        lbr = l.bracket(right)
                        if sfbr:
                            ret+=c*sum(prod(Integer(1)/i for i in range(1 ,j+2 ))*
                                l.T(j+1)._mul_(sfbr[j]) for j in sfbr.keys())
                        if lbr:
                            ret+=c*sum(prod(Integer(1)/i for i in range(1 ,j+2 ))*
                                sf.T(j+1)._mul_(lbr[j]) for j in lbr.keys())
            return ret

        def T(self,n=1):
            if n==0:
                return self
            if n > 1:
                return self.T().T(n-1)
            coef = self.value.monomial_coefficients()
            p = self.parent()
            ret = p.zero()
            for k in coef.keys():
                c = coef[k]
                for i in range(len(k)):
                    jdic = partderiv(k[i])
                    for j in jdic.keys():
                        pt = k.components()
                        pt[i] = j
                        ret += c*jdic[j]*p(pt)
            return ret

        def __getitem__(self, i):
            return self.value.__getitem__(i)

        def _acted_upon_(self, scalar, self_on_left=False):
            scalar_parent = parent(scalar)
            if scalar_parent != self.parent().base_ring():
                if self.parent().base_ring() \
                        .has_coerce_map_from(scalar_parent):
                    scalar = self.parent().base_ring()( scalar )
                else:
                    return None
            if self_on_left:
                return type(self)(self.parent(), self.value * scalar)
            return type(self)(self.parent(), scalar * self.value)

        def monomial_coefficients(self):
            """
            Return the monomial coefficients of ``self`` as a dictionary.
           """
            c = self.value.monomial_coefficients()
            p = self.parent()
            return { p(k) : c[k] for k in c.keys() } 

        def is_monomial(self):
            return ( len(self.monomial_coefficients()) == 1  or self.is_zero() )

        def index(self):
            if self.is_zero():
                return None
            if not self.is_monomial():
                raise ValueError ("index can only be computed for"
                         " monomials")
            return self.value.monomial_coefficients().keys()[0]

        def weight(self):
            p = self.parent()
            if not p.is_graded():
                raise ValueError("Weight is only defined in H-graded"
                    " vertex algebras")
            if self.is_zero():
                return Infinity
            ls = []
            for idx in self.value.monomial_coefficients().keys():
                ret = sum(len(idx[i])*(p._lca.gen(i).degree()-1) + idx[i].size() 
                    for i in range(p.ngens()))
                ls.append(ret)
            if ls[1:] == ls[:-1]:
                return ls[0]
            raise ValueError("{} is not homogeneous!".format(self))

        def li_filtration_degree(self):
            p = self.parent()
            if self.is_zero():
                return Infinity
            ls = []
            for p in self.value.monomial_coefficients().keys():
                ret = p.size() - sum(j.length() for j in p.components())
                ls.append(ret)
            return max(ls)

           

        def PBW_filtration_degree(self):
            p = self.parent()
            if self.is_zero():
                return -1 
            ls = []
            for p in self.value.monomial_coefficients().keys():
                ret = sum(j.length() for j in p.components())
                ls.append(ret)
            return max(ls)

        def _pbw_one_less(self):
            if self.is_zero():
                return (None, None, None)
            svmc = self.value.monomial_coefficients()
            if len(svmc) != 1:
                raise ValueError("{} is not a monomial.".format(self))
            k = svmc.keys()[0]
            c = svmc[k] 
            p = self.parent()
            i = next((j for j,x in enumerate(k) if x), None)
            if i is None:
                return (p.vacuum(), p.vacuum(), c)
            l = [Partition([]), ]*p.ngens()
            l[i] = Partition([k[i].get_part(0)])
            sf = k.to_list()
            sf[i] = Partition(Partition(sf[i]).to_list()[1:])
            return (p(l),p(sf),c) 

        
class VirasoroVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self, R, c, arg0 = None, **kwds):
        from lie_conformal_algebra import VirasoroLieConformalAlgebra
        ML = VirasoroLieConformalAlgebra(R)
        if arg0 is not None:
            c = 1  - 6*(c-arg0)**2/(c*arg0)
        cp = Family({ML.gen(1):c})
        super(VirasoroVertexAlgebra,self).__init__(R, ML,
                 central_parameters=cp, names = "L")
        self._c = c

    def _repr_(self):
        return "The Virasoro vertex algebra at central charge {}".format(self.central_charge())

    def central_charge(self):
        return self._c

class AffineVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self, R, ct, k, **kwds):
        from lie_conformal_algebra import AffineLieConformalAlgebra
        ML = AffineLieConformalAlgebra(R, ct, prefix='E', bracket='(')
        cp = Family({ML.central_elements()[0]: k})
        super(AffineVertexAlgebra,self).__init__(R, ML, 
            central_parameters = cp, **kwds)

        self._level = k
        if type(ct) is str:
            from sage.combinat.root_system.cartan_type import CartanType
            ct = CartanType(ct)
        self._ct = ct
    
    def level(self):
        return self._level

    def cartan_type(self):
        return self._ct

    def is_critical(self):
        return self.level() == -self.cartan_type().dual_coxeter_number()

    def _repr_(self):
        if self.is_critical():
            return "The universal affine vertex algebra of CartanType {} at critical level".format(self.cartan_type())
        else:
            return "The universal affine vertex algebra of CartanType {0} at level {1}".format(self.cartan_type(), self.level())
    
def partderiv(p):
    l = p.to_exp()
    ret = {}
    for i in range(len(l)):
        if l[i] > 0:
            l2 = list(l)
            l2[i] -= 1 
            try:
                l2[i+1]+=1 
            except IndexError:
                l2.append(1)
            ret[Partition(exp=l2)] = l[i]*(i+1)
    return ret
        

