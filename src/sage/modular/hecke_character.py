# -*- coding: utf-8 -*-
r"""
Hecke characters

AUTHORS:

- Robert Harron (2016-08-15): initial version

EXAMPLES::

#*****************************************************************************
#       Copyright (C) 2016 Robert Harron <robert.harron@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""

from sage.structure.factory import UniqueFactory
from sage.arith.srange import srange
from sage.arith.misc import lcm
from sage.rings.integer_ring import ZZ
#from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
from sage.groups.abelian_gps.dual_abelian_group_element import DualAbelianGroupElement
from sage.groups.abelian_gps.dual_abelian_group import DualAbelianGroup_class
from sage.lfunctions.dokchitser import Dokchitser

def _mask_to_list(L):
    ret = []
    for i in srange(len(L)):
        if L[i] != 0:
            ret.append(i)
    return ret

class HeckeCharacter(DualAbelianGroupElement):
    def __call__(self, g):
        R = self.parent().group()
        if g.parent() is R:
            DualAbelianGroupElement.__call__(self, g)
        K = R.number_field()
        g = K.ideal(g)
        if self.parent().modulus().finite_part().is_coprime(g):
            return DualAbelianGroupElement.__call__(self, R(g))
        return self.parent().base_ring().zero()
    
    def _log_values_on_gens(self):
        r"""
        Returns a list of integers `(a_j)` such that the value of ``self`` on the jth
        generator of the ray class group is `\exp(2\pi ia_j/d_j)`, where `d_j` is the
        order of the jth generator.
        """
        #R = self.parent().group()
        #orders = R.gens_orders()
        #val_orders = [self(a).multiplicative_order() for a in R.gens_values()]
        #return [orders[j].divide_knowing_divisible_by(val_orders[j]) for j in range(len(orders))]
        return self.exponents()
    
    def modulus(self):
        return self.parent().modulus()
    
    def level(self):
        return self.modulus()
    
    def conductor(self):
        R = self.parent().group()
        K = R.number_field()
        bnr = R.pari_bnr()
        modulus = bnr.bnrconductorofchar(self._log_values_on_gens())
        #print modulus
        return K.modulus(K.ideal(modulus[0]), _mask_to_list(modulus[1])) #Are the infinite places ordered the same in pari?
    
    def is_primitive(self):
        return self.conductor() == self.modulus()
    
    def primitive_character(self):
        #cond = self.conductor()
        #if cond == self.modulus():
        #    return self
        #R = nf_ray_class_group(cond.number_field(), cond)
        raise NotImplementedError
    
    def extend(self, m):
        if isinstance(m, HeckeCharacterGroup_class):
            pass
        raise NotImplementedError
    
    def analytic_conductor(self):
        N = self.conductor().finite_part().norm()
        return N * self.parent().group().number_field().disc().abs()
    
    def root_number(self):
        return self.parent().group().pari_bnr().bnrrootnumber(self._log_values_on_gens()).sage()
    
    def dirichlet_series_coefficients(self, max_n):
        Idict = self.parent().number_field().ideals_of_bdd_norm(max_n)
        ans = [ZZ(1)] + [ZZ(0)] * (max_n - 1)
        for n in range(2, max_n + 1):
            Is = Idict[n]
            if len(Is) == 0:
                continue
            ans[n-1] = sum(self(I) for I in Is)
        return ans
    
    def Lfunction(self, prec=53):
        r"""
        Return 
        
        EXAMPLES:
        
            A totally odd character of a real quadratic field.
            
            sage: F.<a> = NumberField(x^2 - 5)
            sage: mf = F.modulus(F.ideal(4), [0, 1])
            sage: H = HeckeCharacterGroup(mf)
            sage: chi = H.gens()[1]
            sage: L = chi.Lfunction()
            sage: [L(-n) for n in range(3)]
            [1.00000000000000, 0.000000000000000, 15.0000000000000]
        """
        #Figure out Gamma factors for more general characters
        gamma_factors = [0] * self.parent().number_field().degree()
        for i in self.conductor().infinite_part():
            gamma_factors[i] = 1
        L = Dokchitser(self.analytic_conductor(), gamma_factors, 1, self.root_number())
        L.init_coeffs(self.dirichlet_series_coefficients(L.num_coeffs()))
        return L

class HeckeCharacterGroup_class(DualAbelianGroup_class):
    r"""
    EXAMPLES::
    
        sage: F.<a> = NumberField(x^2 - 5)
        sage: mf = F.modulus(F.prime_above(5) * F.prime_above(29), [0,1])
        sage: H = HeckeCharacterGroup(mf); H
        Group of finite order Hecke characters modulo (Fractional ideal (-11/2*a - 5/2)) * infinity_0 * infinity_1
        sage: [[chi(F.ideal(31)), chi(F.ideal(-12672))] for chi in H.gens()]
        [[zeta4, 1], [1, -1]]
    """
    Element = HeckeCharacter
    
    def __init__(self, ray_class_group, base_ring=None, names=None):
        if base_ring is None:
            from sage.rings.number_field.number_field import CyclotomicField
            base_ring = CyclotomicField(lcm(ray_class_group.gens_orders()))
        if names is None:
            names = 'chi'
        DualAbelianGroup_class.__init__(self, ray_class_group, names, base_ring)
    
    def _repr_(self):
        return 'Group of finite order Hecke characters modulo ' + str(self.modulus())
    
    def modulus(self):
        return self.group().modulus()
    
    def level(self):
        return self.modulus()
    
    def number_field(self):
        return self.group().number_field()
    
    def element_from_values_on_gens(self, vals):
        gens_orders = self.gens_orders()
        if len(vals) != len(gens_orders):
            raise ValueError("Incorrect number of values specified. %s specified, but needed %s"%(len(vals), len(gens_orders)))
        exponents = [gens_orders[i].divide_knowing_divisible_by(vals[i].multiplicative_order()) if vals[i] != 1 else ZZ.zero() for i in range(len(gens_orders))]
        #print "in elem..."
        #print vals
        #print gens_orders
        #print exponents
        return self.element_class(self, exponents)
    #def _element_constructor_(self, *args, **kwds):
    #    if isinstance(args[0], basestring):
    #        raise TypeError("Wrong type to coerce into HeckeCharacterGroup.")
    #    try:
    #        n = len(args[0])
    #    except TypeError:
    #        return DualAbelianGroup_class._element_constructor_(self, *args, **kwds)
    #    gens_orders = self.gens_orders()
    #    if n != len(gens_orders):
    #        return DualAbelianGroup_class._element_constructor_(self, *args, **kwds)
    #    exponents = [gens_orders[i].divide_knowing_divisible_by(args[0][i].multiplicative_order()) for i in range(n)]
    #    return self.element_class(self, exponents)

class HeckeCharacterGroupFactory(UniqueFactory):
    def create_key(self, modulus, base_ring=None, names=None):#(m, base_ring=None, names=None):
        #from sage.structure.category_object import normalize_names
        #names = normalize_names()
        if names is None:
            names = 'chi'
        return (base_ring, modulus, names)
    
    def create_object(self, version, key, **extra_args):
        base_ring, modulus, names = key
        return HeckeCharacterGroup_class(modulus.number_field().ray_class_group(modulus), base_ring, names)

HeckeCharacterGroup = HeckeCharacterGroupFactory("HeckeCharacterGroup")
