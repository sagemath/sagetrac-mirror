r"""
Lie conformal algebra element class
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
from sage.functions.other import binomial, factorial
from sage.misc.misc_c import prod
from sage.structure.element import parent
from sage.rings.integer import Integer

class LieConformalAlgebraElementWrapper(ElementWrapper):

    def _repr_(self):
        """
        Return a string representation of ``self``.
       """
        return repr(self.value)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.
       """
        from sage.misc.latex import latex
        return latex(self.value)

    def _ascii_art_(self):
        """
        Return an ascii art representation of ``self``.
       """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art(self.value)

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.
       """
        from sage.typeset.unicode_art import unicode_art
        return unicode_art(self.value)

    def _add_(self, right):
        """
        Add ``self`` and ``rhs``.
       """
        return type(self)(self.parent(), self.value + right.value)

    def _sub_(self, right):
        """
        Subtract ``self`` and ``rhs``.
       """
        return type(self)(self.parent(), self.value - right.value)

    # Need to bypass the coercion model
    def _mul_(left, right):
        return left.lift()*right

    def __neg__(self):
        """
        Return the negation of ``self``.
       """
        return type(self)(self.parent(), -self.value)

    def __getitem__(self, i):
        """
        Redirect the ``__getitem__()`` to the wrapped element.
       """
        return self.value.__getitem__(i)

    def _acted_upon_(self, scalar, self_on_left=False):
        """
        Return the product ``scalar`` times this element of the Lie conformal
        algebra
        """
        scalar_parent = parent(scalar)
        if scalar_parent != self.parent().base_ring():
            if self.parent().base_ring().has_coerce_map_from(scalar_parent):
                scalar = self.parent().base_ring()( scalar )
            else:
                return None
        if self_on_left:
            return type(self)(self.parent(), self.value * scalar)
        return type(self)(self.parent(), scalar * self.value)

    def monomial_coefficients(self):
        """
        Return the monomial coefficients of this element as a dictionary.

        The keys are element of the Lie conformal algebra

        EXAMPLES::

            sage: V = VirasoroLieConformalAlgebra(QQ); V.inject_variables()
            Defining L, C
            sage: v = L + 2*L.T() + 3/2*C; v.monomial_coefficients()
            {TL: 2, L: 1, C: 3/2}

        """
        p = self.parent()
        return { p.monomial(k):v for k,v in 
                self.value.monomial_coefficients().items() }
             

class LCAStructureCoefficientsElement(LieConformalAlgebraElementWrapper):

    """
    An element of a Lie conformal algebra given by structure coefficients.
    """

    def _bracket_(self, right):
        """
        Return the lambda bracket of these two elements.

        The result is a dictionary with non-negative integer keys. The value
        corresponding to the entry `j` is ``self_{(j)}right``

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ)
            sage: L.bracket(L)
            {0: TL, 1: 2*L, 3: 1/2*C}
            sage: L.T().bracket(L)
            {1: -TL, 2: -4*L, 4: -2*C}

        When the Lie conformal algebra is actually a vertex algebra, elements
        are denoted in a different fashion::

            sage: V = VirasoroVertexAlgebra(QQ,1/2)
            sage: V.inject_variables()
            Defining L
            sage: L.bracket(L*L)
            {0: L_-5|0>+2*L_-3L_-2|0>,
             1: 4*L_-2L_-2|0>,
             2: 3*L_-3|0>,
             3: 17/2*L_-2|0>,
             5: 3/2*|0>}

        """
        p = self.parent()
        l1 = self.value.monomial_coefficients()
        l2 = right.value.monomial_coefficients()
        s_coeff = p.structure_coefficients()
        ret = {}
        for i in l1.keys():
            if i[0] in p._central_elements:
                continue
            for j in l2.keys():
                if j[0] in p._central_elements:
                    continue
                if p._index_to_pos[i[0]] > p._index_to_pos[j[0]]:
                    try:
                        v = { k[0]: k[1] for k in s_coeff[j[0],i[0]] }
                    except KeyError:
                        v = {}
                    maxpole = max(v.keys())
                    vals={} 
                    for k in range(maxpole+1):
                        #Do we need characteristic zero or is it the fact that
                        #we chose divided powers as basis?
                        kth_product = sum ( (-1)**(k+j+1)*prod(
                            Integer(1)/i for i in range(1 ,j+1))*
                            v[j+k].T(j) for j in range(maxpole+1 -k) if
                            (j+k) in v.keys() )
                        if kth_product is not p.zero():
                            vals[k]=kth_product
                    myvals = tuple( vals.items() ) 
                else:
                    try:
                        myvals = s_coeff[(i[0],j[0])]
                    except KeyError:
                        myvals = tuple()
                for l in range(j[1]+1):
                    for b in myvals:
                        prev=ret.get(b[0]+j[1]+i[1]-l, p.zero())
                        ret[b[0]+j[1]+i[1]-l]=prev+ \
                            prod(range(b[0]+1 ,b[0]+j[1]+i[1]-l+1))*\
                            (-1)**i[1]*binomial(j[1],l)*\
                            factorial(i[1])**(-1)*factorial(j[1])**(-1)*\
                            l1[i]*l2[j]*b[1].T(l) 
        return { k:ret[k] for k in ret.keys() if ret[k] }
        

    def __getitem__(self, i):
        """
        Return the coefficient of the basis element indexed by ``i``.
       """
        return self.value[self.parent()._indices.index(i)]

    def _repr_(self):
        if self == self.parent().zero():
            return "0";
        coeff = self.value.monomial_coefficients().items()
        ret = ""
        for i in range(len (coeff)):    
            #TODO: deal with this without using an ordering of the ring. 
            if i > 0  and coeff[i][1] > 0 :
                ret += "+"
            if coeff[i][1] < 0 :
                ret += "-"
    
            if abs(coeff[i][1]) != 1 :
                ret += "{}*".format(abs(coeff[i][1]))
            if coeff[i][0][1] == 1 :
                ret += "T"
            elif coeff[i][0][1] > 1 :
                ret += "T^({})".format(coeff[i][0][1])
            ret += "{}".format(coeff[i][0][0])
        return ret
    
    def _latex_(self):
        if self == self.parent().zero():
            return "0";
        coeff = self.value.monomial_coefficients().items()
        ret = ""
        for i in range(len (coeff)):    
            #TODO: deal with this without using an ordering of the ring. 
            if i > 0  and coeff[i][1] > 0 :
                ret += "+"
            if coeff[i][1] < 0 :
                ret += "-"
    
            if abs(coeff[i][1]) != 1 :
                ret += "{}*".format(abs(coeff[i][1]))
            if coeff[i][0][1] == 1 :
                ret += "T"
            elif coeff[i][0][1] > 1 :
                ret += "T^{("+"{}".format(coeff[i][0][1])+")}"
            ret += "{}".format(coeff[i][0][0])
        return ret

 

