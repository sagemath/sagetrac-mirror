r"""
Quotients of vertex algebras 
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

from sage.structure.element_wrapper import ElementWrapper
from sage.structure.element import parent
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.modules import Modules
from sage.misc.cachefunc import cached_method 

class VertexAlgebraQuotient_space_element(ElementWrapper):
    def _repr_(self):
        V = self.parent().ambient()
        M = V.module()
        return repr(V(M._from_dict(self.value[1].monomial_coefficients())))

    def _add_(self,right):
        p = self.parent()
        m = max(self.value[0], right.value[0]) 
        S = p._submodule.get_degree(m)
        A = S.ambient()
        Q = A.quotient_module(S)
        #coerce self and right to A
        qx = self.value[1].parent()
        x = A._from_dict(qx.lift(self.value[1]).monomial_coefficients())
        qy = right.value[1].parent()
        y = A._from_dict(qy.lift(right.value[1]).monomial_coefficients())
        return p((m,Q.retract(x+y)))

    def _sub_(self,right):
        p = self.parent()
        m = max(self.value[0], right.value[0]) 
        S = p._submodule.get_degree(m)
        A = S.ambient()
        Q = A.quotient_module(S)
        #coerce self and right to A
        qx = self.value[1].parent()
        x = A._from_dict(qx.lift(self.value[1]).monomial_coefficients())
        qy = right.value[1].parent()
        y = A._from_dict(qy.lift(right.value[1]).monomial_coefficients())
        return p((m,Q.retract(x-y)))

    def _acted_upon_(self,scalar, self_on_left=False):
        scalar_parent = parent(scalar)
        if scalar_parent != self.parent().base_ring():
            if self.parent().base_ring() \
                    .has_coerce_map_from(scalar_parent):
                scalar = self.parent().base_ring()( scalar )
            else:
                return None
        if self_on_left:
            return type(self)(self.parent(), (self.value[0], 
                                            self.value[1]*scalar))
        return type(self)(self.parent(),
                            (self.value[0],scalar*self.value[1]))

    def monomial_coefficients(self):
        p = self.value[1].parent()
        return {type(self)(self.parent(),(self.value[0],p(k))):v for k,v in 
                            self.value[1].monomial_coefficients().items()}

    def length(self):
        return len(self.monomial_coefficients())

    def is_monomial(self):
        return self.length() == 1 

    def __getitem__(self,i):
        return self.value[1].__getitem__(i)

    def filtered_degree(self):
        return self.value[0]

    def __nonzero__(self):
        return bool(self.value[1])

    def lift(self):
        return self.parent().lift(self)
    

class VertexAlgebraQuotient_space(Parent, UniqueRepresentation):

    def __init__(self, submodule, **kwds):
        category = kwds.get('category', Modules(submodule.base_ring()))
        names = kwds.get('names', None)
        super(VertexAlgebraQuotient_space,self).__init__(base = submodule.base_ring(), 
                    category=category, names=names  )
        self._submodule = submodule
        self._ambient = submodule.ambient()

    Element = VertexAlgebraQuotient_space_element

    def ambient(self):
        return self._ambient

    def lift(self, x):
        assert x in self
        V = self._ambient
        return V(V.module()._from_dict(x.value[1].lift().\
                    monomial_coefficients()))
    
    def retract(self, x):
        if x == 0:
            n = 0 
        else:
            n = x.filtered_degree()
        S = self._submodule.get_degree(n)
        A = S.ambient()
        y = A._from_dict(x.value.monomial_coefficients())
        Q = A.quotient_module(S)

        return self.element_class(self, (n,Q.retract(S.reduce(y))))

    def _element_constructor_(self,v):
        if v == self.base_ring().zero():
            return self.zero()
        if v in self._ambient:
            return self.retract(self._ambient(v))
        return self.element_class(self, v)

    def _coerce_map_from_(self,other):
        if self._ambient.has_coerce_map_from(other):
            return True
    
    @cached_method
    def li_filtration(self,n,m=None):
        S = self._submodule.get_graded_part(n)
        A = S.ambient()
        Q = A.quotient_module(S)
        F = self._ambient.li_filtration(n,m)
        if m is not None:
            F = {m:F}
        G = {m: Q.submodule([Q.retract(b.lift()) for b in F[m].basis()]) for
                m in F.keys() }
        if m is not None:
            G = G[m]
        return G
    
    @cached_method
    def get_graded_part(self,n,m=None):
        S = self._submodule.get_graded_part(n)
        A = S.ambient()
        Q = A.quotient_module(S)
        if m == None:
            return Q
        G = self.li_filtration(n,m)
        if m >= n+1:
            return G
        H = self.li_filtration(n,m+1)
        return G.quotient_module([G.retract(b.lift()) for b in H.basis()])
    
    def dimension(self,n):
        return self.get_graded_part(n).dimension()

    @cached_method
    def basis(self, n):
        return [ self((n,b)) for b in self.get_graded_part(n).basis()]
    
    @cached_method
    def zero(self):
        return self.retract(self._ambient.zero())

        
class VertexAlgebraQuotient(VertexAlgebraQuotient_space):
    def __init__(self, ideal):
        category = ideal.ambient().category().Quotients()
        try:
            names = ideal.ambient().variable_names()
        except ValueError:
            names = None

        super(VertexAlgebraQuotient,self).__init__(ideal,
            names = names, category=category)
        
    def _repr_(self):
        return "Quotient of {0} by the ideal generated by {1}".format(self._ambient,
                self._submodule._gens)

    def module(self):
        return self

    def vacuum(self):
        return self.retract(self._ambient.vacuum())

    def gens(self):
        return self._ambient.gens()

    def gen(self,i):
        return self.gens()[i]

    class Element(VertexAlgebraQuotient_space_element):
        def _mul_(self,right):
            p = self.parent()
            return p.retract(self.lift()._mul_(right.lift()))
        def _bracket_(self, other):
            br = self.lift()._bracket_(other.lift())
            return { k:self.parent().retract(v) for k,v in br.items() }

        def T(self,n=1):
            return self.parent().retract(self.lift().T(n))

        def weight(self):
            return self.lift().weight()

               



        

