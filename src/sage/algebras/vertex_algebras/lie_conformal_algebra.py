r"""
Lie Conformal Algebras
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
from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from sage.algebras.vertex_algebras.lie_conformal_algebra_element import \
        LCAStructureCoefficientsElement
from sage.sets.family import Family 
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra
from sage.categories.cartesian_product import cartesian_product
from sage.combinat.free_module import CombinatorialFreeModule
from sage.functions.other import binomial
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer

class LieConformalAlgebra(Parent, UniqueRepresentation):
    r""" 
    A Lie conformal algebra `L`. pingoloeon termina
    """
    # As of now either pass vertex=A vertex algebra 
    # or pass arg0 = a dict with the non-vanishig OPE
    @staticmethod
    def __classcall_private__(cls, R=None, arg0=None, names=None, 
        central_elements=None, index_set=None, category=None, 
        weights=None, **kwds):
        
        if R is None:
            raise ValueError("First argument must be a ring!")

        if isinstance(arg0, dict):
            graded=kwds.get("graded", False)
            if weights is not None or graded:
                #if not arg0:
                    #return AbelianLieConformalAlgebra(R, names, index_set,
                    #    weights, graded=True)
                return GradedLieConformalAlgebra(R, arg0,
                     names=names, index_set=index_set, 
                    central_elements=central_elements,category=category,weights=weights, **kwds)
            else:
                #if not arg0:
                #    return AbelianLieConformalAlgebra(R, names, index_set)
                return LieConformalAlgebraWithStructureCoefficients(R, arg0,
                     names=names, index_set=index_set,central_elements=central_elements,category=category, **kwds)

    def __init__(self, R, names=None, category=None):
        category = LieConformalAlgebras(R).or_subcategory(category)
        super(LieConformalAlgebra,self).__init__(base=R, names=names, 
            category=category)

    def _element_constructor_(self,x):
        if x in self.base_ring():
            if x != 0:
                raise ValueError("can only convert the scalar 0 into a Lie conformal algebra element")
            return self.zero()
        try: 
            x = self.module()(x)
        except TypeError:
            raise ValueError("Don't know how to convert {0} into an element of "\
                             "{1}".format(x,self))
        return self.element_class(self,x)

    @cached_method
    def zero(self):
        return self.element_class(self, {})

    def monomial(self,i):
        return self.element_class(self, {i:self.base_ring().one()})

class LieConformalAlgebraWithBasis(LieConformalAlgebra):
    """
    A Lie conformal algebra with distinguished basis. 
    """
    def __init__(self,R, names=None, index_set=None, 
        category=None, prefix='B', **kwds):
        self._indices = index_set
        category = LieConformalAlgebras(R).WithBasis().or_subcategory(category)
        super(LieConformalAlgebraWithBasis,self).__init__(R, names, category)

    @cached_method
    def lie_conformal_algebra_basis(self):
        return Family(self._indices, self.monomial, name="basis map")

    def indices(self):
        return self._indices

    def module(self):
        return CombinatorialFreeModule(self.base_ring(), basis_keys=self._indices,
            category=self.category())

    def monomial(self,i):
        B = self.module().basis()
        return self.element_class(self, self.module()(B[i]))
    
    def zero(self): 
        return self.element_class(self, self.module().zero())

class LieConformalAlgebraWithGenerators(LieConformalAlgebraWithBasis):
    """
    A Lie conformal algebra with distinguished generators as `R[T]` module. 
    We now only accept direct sums of free modules plus some generators `C` 
    such that `TC = 0`.
    """
    def __init__(self,R, names=None, index_set=None, central_elements=None,      
            category=None, prefix='B', **kwds):
        self._generators = tuple(index_set)

        E = cartesian_product([index_set, NonNegativeIntegers()])
        if central_elements is not None:
            self._generators = self._generators + tuple(central_elements)
            E = DisjointUnionEnumeratedSets((E, cartesian_product([
                tuple(central_elements), {0}])))
    
        super(LieConformalAlgebraWithGenerators,self).__init__(R,names=names, 
            index_set=E, category=category, prefix=prefix, **kwds)

        self._central_elements = tuple(central_elements)

    @cached_method
    def lie_conformal_algebra_generators(self):
        return Family(self._generators, lambda i: self.monomial((i,Integer(0))), 
                name = "generator map")

    @cached_method
    def gens(self):
        G = self.lie_conformal_algebra_generators()
        try:
            return tuple(G[i] for i in self._generators)
        except (KeyError, ValueError):
            return tuple(G)

    def gen(self,i):
        return self.gens()[i]

    def central_elements(self):
        G = self.lie_conformal_algebra_generators()
        return tuple( G[i] for i in self._central_elements)


from sage.structure.indexed_generators import (IndexedGenerators, 
        standardize_names_index_set)

class FinitelyGeneratedLieConformalAlgebra(
        LieConformalAlgebraWithGenerators):
    """
    A finitely generated Lie conformal algebra.
    """
    def __init__(self,R, names=None, index_set=None, central_elements=None,      
            category=None, prefix='B', **kwds):
        super(FinitelyGeneratedLieConformalAlgebra,self).__init__(R,names=names,
            index_set=index_set,central_elements=central_elements, 
            category=category, prefix=prefix, **kwds)
        self.__ngens = len(self._generators)

    def _repr_(self):
        if self.__ngens == 1:
            return "Lie conformal algebra on the generator {0} over {1}".format(
                self.gen(0), self.base_ring())
        return "Lie conformal algebra on {0} generators {1} over {2}.".format(
            self.__ngens, self.gens(), self.base_ring()) 

    @lazy_attribute
    def _ordered_generators(self):
        return tuple(self.lie_conformal_algebra_generators().keys())

    def _an_element_(self):
        return self.sum(self.lie_conformal_algebra_generators())

    def ngens(self):
        return self.__ngens
 
class LieConformalAlgebraWithStructureCoefficients(
    FinitelyGeneratedLieConformalAlgebra, IndexedGenerators):
    r"""
    A Lie conformal algebra with a set of specified structure coefficients.

    The structure coefficients are specified as a dictionary `d` whose
    keys are pairs of generator indices, and whose values are
    dictionaries of dictionaries. A typical entry is of the form `{n:val}` 
    where `n` is a positive number (the `n`-th product) and where `val` is 
    itself a dictionary, whose keys are in turn indexed by pairs `(w,m)` of a
    non-negative number `m` and a generator index `w`. The value of `val` for the
    entry `(w,m)` (for `w` a generator index) is the coefficient of
    `T^{(m)}b_w` in the OPE  `b_{u(n)} b_v` (where `b_x` means the generator
    element with index `x`).
    """

    @staticmethod
    def __classcall__(cls,R,s_coeff, names=None,
            index_set=None, central_elements=None, **kwds):

        names, index_set = standardize_names_index_set(names,index_set)

        if names is not None and names != tuple(index_set):
            names2 = names + tuple(central_elements)
            index_set2 = DisjointUnionEnumeratedSets((index_set,
                Family(tuple(central_elements))))
            d = {x: index_set2[i] for i,x in enumerate(names2)}
            try:
                s_coeff = {(d[k[0]],d[k[1]]): {a : {(d[x],y): 
                    s_coeff[k][a][(d[x],y)] for (x,y) in 
                    s_coeff[k][a] } for a in s_coeff[k] } for k in s_coeff}

            except KeyError:
                # At this point we assume they are given by the index set
                pass

        if names is None:
            names = 'l' 
            ngens=index_set.cardinality()
        else:
            ngens= len(names)

        names,index_set = standardize_names_index_set(names, index_set,ngens)

        s_coeff = LieConformalAlgebraWithStructureCoefficients._standardize_s_coeff(s_coeff, index_set, central_elements)
   
 
        if central_elements is not None:
            names += tuple(central_elements)
        #if s_coeff.cardinality() == 0:
         #   return AbelianLieConformalAlgebra(R, names, index_set, **kwds)

        return super(LieConformalAlgebraWithStructureCoefficients, cls).__classcall__(
            cls, R, s_coeff, names=names, index_set=index_set, 
                central_elements=central_elements,  **kwds)

    @staticmethod
    def _standardize_s_coeff(s_coeff, index_set, ce):
        index_to_pos = {k: i for i,k in enumerate(index_set)}

        sc = {}
        #k has a pair of generators
        for mypair in s_coeff.keys():
            #e.g.  v = { 0: { (L,2):3, (G,3):1}, 1:{(L,1),2} }
            v = s_coeff[mypair]
            if index_to_pos[mypair[0]] > index_to_pos[mypair[1]]:
                key=(mypair[1],mypair[0])
                maxpole = max(v.keys())
                vals={} 
                for k in range(maxpole+1):
                    kth_product = {}
                    for j in range(maxpole+1-k):
                        if k+j in v.keys():
                            for i in v[k+j].keys():
                                if ( i[0] not in ce ) or (i[0] in
                                    ce and i[1] + j == 0):
                                    kth_product[(i[0],i[1]+j)] = kth_product.get((i[0], i[1]+j), 0)
                                    kth_product[(i[0],i[1]+j)] += v[k+j][i]*(-1)**(k+j+1)*binomial(i[1]+j,j)
                    kth_product = { k:v for k,v in kth_product.items() if v }
                    if kth_product is not {}:
                        vals[k]=kth_product
            else:
                if not index_to_pos[mypair[0]] < index_to_pos[mypair[1]]:
                    if mypair[0] == mypair[1]:
                        #Here I should check skew-symmetry!
                        pass
                key = tuple(mypair)
                vals={}
                for l in v.keys():
                    lth_product = { k:y for k,y in v[l].items() if y }
                    if lth_product is not {}:
                        vals[l]=lth_product

            myvals = tuple( (i , tuple( (x,v) for x,v in vals[i].items()) ) 
                        for i in vals.keys() if vals[i]  )
            
            if key in sc.keys():
                if sorted(sc[key]) != sorted(myvals): 
                        raise ValueError("two distinct values given for one and the same bracket")

            if myvals:
                sc[key] = myvals

        return Family(sc)

    def __init__(self, R, s_coeff, names=None, index_set=None,
            category=None, central_elements=None,  bracket=False,
            latex_bracket=False, string_quotes=False, prefix='', **kwds):
        
        self._index_to_pos = {k: i for i,k in enumerate(index_set)}
        if "sorting_key" not in kwds:
            kwds["sorting_key"] = self._index_to_pos.__getitem__

        category = LieConformalAlgebras(R).WithBasis().FinitelyGenerated().or_subcategory(category)
        FinitelyGeneratedLieConformalAlgebra.__init__(self, R, names=names, 
            index_set=index_set, category=category,
            central_elements=central_elements, prefix=prefix, **kwds)
        IndexedGenerators.__init__(self, self._indices, prefix=prefix,
                                   bracket=bracket, latex_bracket=latex_bracket,
                                   string_quotes=string_quotes, **kwds)

        #at this stage we have access to self.module(). We convert the values
        #of s_coeff to be elements in self.
   
        self._s_coeff = Family({ k: tuple( (p[0], sum( c[1]*self.monomial(c[0]) for
                c in p[1] )) for p in s_coeff[k] ) for k in s_coeff.keys() })

    _repr_term = IndexedGenerators._repr_generator
    _latex_term = IndexedGenerators._latex_generator

    def structure_coefficients(self):
        return self._s_coeff

    class Element(LCAStructureCoefficientsElement):
        """
        """ 
   
class GradedLieConformalAlgebra(LieConformalAlgebraWithStructureCoefficients):
    def __init__(self, R, s_coeff, names=None, index_set=None,
            category=None, central_elements=None,  bracket=False,
            latex_bracket=False, string_quotes=False, prefix='', 
            weights=None, **kwds):
        category = LieConformalAlgebras(R).HGraded().WithBasis().FinitelyGenerated().or_subcategory(category)
        LieConformalAlgebraWithStructureCoefficients.__init__(self,R, s_coeff,
            names=names, category=category, central_elements=central_elements, 
            bracket=bracket, latex_bracket=latex_bracket, index_set=index_set,
            string_quotes=string_quotes, prefix=prefix, **kwds)
        if weights is None:
            weights = (1,)* (len(self._generators) - 
                    len(self.central_elements()))
        if len (weights) != (len(self._generators) - 
                    len(self.central_elements())):
            raise ValueError("weights and (non-central) generator lists must be of same length")
        self._weights = weights

    class Element(LCAStructureCoefficientsElement):
        def degree(self):
            if self.is_zero():
                return Infinity
            p = self.parent()
            ls = []
            for idx in self.value.monomial_coefficients().keys():
                ret = idx[1]
                if idx[0] not in p._central_elements:
                    ret+= p._weights[p._index_to_pos[idx[0]]] 
                ls.append(ret)
            if ls[1:] == ls[:-1]:
                return ls[0]
            raise ValueError("{} is not homogeneous!".format(self))
 
def VirasoroLieConformalAlgebra(R):
    
    virdict =  {('L','L'): {0: {('L',1):1}, 1: {('L',0): 2}, 
                3:{('C', 0): R(2).inverse_of_unit()}}}
    return LieConformalAlgebra(R, virdict, names = ('L',), 
        central_elements = ('C'), weights = (2,))

def AffineLieConformalAlgebra(R, ct, **kwds):
    if type(ct) is str: 
        from sage.combinat.root_system.cartan_type import CartanType
        ct = CartanType(ct)
    if not ( ct.is_finite() or ct.is_irreducible ):
        raise ValueError("Only affine algebras of simple finite dimensional"
            "Lie algebras are implemented")
    hv = Integer(ct.dual_coxeter_number())
    g = LieAlgebra(R, cartan_type=ct)
    B = g.basis()
    gdict = {}
    for k1 in B.keys():
        for k2 in B.keys():
            myb = B[k1].bracket(B[k2]).monomial_coefficients()
            myf = R(2).inverse_of_unit()*R(hv).inverse_of_unit()*g.killing_form(B[k1],B[k2])
            if myb or myf:
                gdict[(k1,k2)] = {}
                if myb:
                    gdict[(k1,k2)][0] = { (nk,0): myb[nk] for nk in
                        myb.keys() } 
                if myf:
                    gdict[(k1,k2)][1] = { ('K', 0): myf }

    weights = (1,)*B.cardinality()
    return LieConformalAlgebra(
        R, gdict, index_set=B.keys(),  
        central_elements=('K', ), weights = weights, 
        **kwds)

    


