"""
Quotients of Vertex Algebras.

AUTHORS

- Reimundo Heluani (08-09-2019): Initial implementation.
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

from sage.misc.cachefunc import cached_method
from sage.rings.all import Integer
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.vertex_algebras import VertexAlgebras
from sage.modules.with_basis.indexed_element import IndexedFreeModuleElement
from .energy_partition_tuples import EnergyPartitionTuples_all
from sage.sets.family import Family

class VertexAlgebraQuotientBasis(EnergyPartitionTuples_all):
    def __init__(self,I):
        vgens = I.ambient().gens()
        weights = [g.degree() for g in vgens]
        regular = [2*g.is_even_odd() for g in vgens]
        super(VertexAlgebraQuotientBasis,self).__init__(weights, len(weights),
                                                        regular=regular)
        self._ideal = I

    def _repr_(self):
        return "Basis of the quotient of {} by {}".format(self._ideal.ambient(),
                                                          self._ideal)

    def _element_constructor_(self,x):
        ambient_indices = self._ideal.ambient().indices()
        if not x in ambient_indices:
            raise ValueError("Do not know how to convert {} into {}".format(x,
                             self))

        x = ambient_indices(x)
        if self._ideal._inverse_on_support(x) is None:
                return self.element_class(self,x)
        raise ValueError("Do not know how to convert {} into {}".format(x,self))

    def __contains__(self,x):
        try:
            self(x)
        except ValueError:
            return False
        return True

    def __iter__(self):
        for i in EnergyPartitionTuples_all.__iter__(self):
            if self._ideal._inverse_on_support(i) is None:
                yield(self.element_class(self,i))

    def __getitem__(self,r):
        if isinstance(r,self.element_class):
            return r
        if isinstance(r,(int,Integer)):
            return self.unrank(r)
        elif isinstance(r,slice):
            start=0 if r.start is None else r.start
            stop=r.stop
            if stop is None:
                raise ValueError('infinite set')
            count=0
            parts=[]
            for t in self:
                if count==stop:
                    break
                if count>=start:
                    parts.append(t)
                count+=1
            if count==stop or stop is None:
                return parts
            raise IndexError('value out of range')
        raise NotImplementedError('Do not know how to look for {}'.format(r))

    def subset(self, energy=None):
        if energy==None:
            return self
        return Family([self(i) for i in EnergyPartitionTuples_all.subset(self,
                                            energy=energy) if i in self ])

class VertexAlgebraQuotientElement(IndexedFreeModuleElement):

    def _repr_(self):
        return repr(self.lift())

    def lift(self):
        return self.parent().lift(self)

    def is_even_odd(self):
        #__mro__ picks Modules before VertexAlgebras.Quotients
        return self.lift().is_even_odd()

class VertexAlgebraQuotient(CombinatorialFreeModule):

    def __init__(self, ideal, category=None):
        default_category = VertexAlgebras(ideal.category().base_ring())\
                            .Quotients()
        category = default_category.or_subcategory(category, join=True)
        self._ideal = ideal
        self._ambient = ideal.ambient()
        indices = VertexAlgebraQuotientBasis(ideal)

        try:
            names = self._ambient.variable_names()
        except AttributeError:
            names = None

        CombinatorialFreeModule.__init__(self,ideal.base_ring(),
                                    basis_keys=indices,
                                    element_class=VertexAlgebraQuotientElement,
                                    category=category, names=names)

    def _element_constructor_(self,x):
        if x in self._ambient:
            return self.retract(x)
        return super(VertexAlgebraQuotient,self)._element_constructor_(x)
    def _repr_(self):
        return "Quotient of {} by {}".format(self._ambient,self._ideal)

    def cover_algebra(self):
        """
        The covering vertex algebra of this quotient.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ, 1/2); Q=V.quotient(V.ideal(V.find_singular(6)[0]))
            sage: Q.module()
            Quotient of The Virasoro vertex algebra of central charge 1/2 by the ideal generated by (L_-2L_-2L_-2|0>-33/8*L_-4L_-2|0>+93/64*L_-3L_-3|0>-27/16*L_-6|0>,)
            sage: Q.module().ambient()
            The Virasoro vertex algebra of central charge 1/2
            sage: Q.module().ambient() is V
            True

        """
        return self._ambient

    def defining_ideal(self):
        return self._ideal

    def lift(self, x):
        r"""Return an element on the ambient vertex algebra in the preimage of
        ``x`` by the quotient map

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ, 1/2); Q=V.quotient(V.ideal(V.find_singular(6)[0]))
            sage: L = Q(Q.gen(0))
            sage: v = L*(L*L)
            sage: v
            33/8*L_-4L_-2|0>-93/64*L_-3L_-3|0>+27/16*L_-6|0>
            sage: v.parent()
            Quotient of The Virasoro vertex algebra of central charge 1/2 by the ideal generated by (L_-2L_-2L_-2|0>-33/8*L_-4L_-2|0>+93/64*L_-3L_-3|0>-27/16*L_-6|0>,)
            sage: v.lift().parent()
            The Virasoro vertex algebra of central charge 1/2

        """
        assert x in self
        return self._ambient._from_dict(x._monomial_coefficients)

    def retract(self, x):
        r"""Let this vertex algebra `Q` be the quotient of `V` by the ideal `I`
        and ``x`` be an element of `V`. This method returns the reduction of `x`
        modulo `I`. That is the image of ``x`` under the quotient map.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ, 1/2); Q=V.quotient(V.ideal(V.find_singular(6)[0]))
            sage: L = V.gen(0)
            sage: v = L*(L*L)
            sage: Q.retract(v)
            33/8*L_-4L_-2|0>-93/64*L_-3L_-3|0>+27/16*L_-6|0>
            sage: Q.retract(v).parent()
            Quotient of The Virasoro vertex algebra of central charge 1/2 by the ideal generated by (L_-2L_-2L_-2|0>-33/8*L_-4L_-2|0>+93/64*L_-3L_-3|0>-27/16*L_-6|0>,)
            sage: v.parent()
            The Virasoro vertex algebra of central charge 1/2

        """
        return self._from_dict(self._ideal.reduce(x)._monomial_coefficients)

    def singular_support(self):
        P = self._ambient.singular_support()
        gens = [g.li_filtration_lt() for g in self._ideal.gens()]
        gens = tuple([P(g) for g in gens if g.li_filtration_degree() == 0])
        I = P.ideal(gens)
        return P.quotient(I)

