r"""
Quotient rings for orders inside number fields.

Implements functionality for quotients of orders inside number fields,
that do not hold for general quotients. The classes in this modulo build
upon the classes in sage.rings.quotient_Ring and sage.rings.quotient_ring_element.

AUTHORS:

- Joey van Langen (2017-06-19): initial version

EXAMPLES::
    
        sage: K.<a> = NumberField(x^2 - x + 3)
        sage: R = K.ring_of_integers()
        sage: S = R.quotient(12); S
        Quotient of Maximal Order in Number Field in a with defining polynomial x^2 - x + 3 by the ideal (12)
        
    Another example with a non-principal ideal::
    
        sage: K.<a> = NumberField(x^2 - x + 6)
        sage: R = K.ring_of_integers()
        sage: I = K.prime_above(3) * K.prime_above(7)
        sage: S = R.quotient(I); S
        Quotient of Maximal Order in Number Field in a with defining polynomial x^2 - x + 6 by the ideal (21, 7*a)
        
    One can use this quotient for iteration::
    
        sage: K.<a> = NumberField(x^4 + x^3 + x^2 + x + 1)
        sage: R = K.ring_of_integers()
        sage: S = R.quotient(4)
        sage: [b for b in S if b^2 == 0]
        [0, 2*a^3, 2*a^2, 2*a^3 + 2*a^2, 2*a, 2*a^3 + 2*a, 2*a^2 + 2*a, 2*a^3 + 2*a^2 + 2*a, 2, 2*a^3 + 2, 2*a^2 + 2, 2*a^3 + 2*a^2 + 2, 2*a + 2, 2*a^3 + 2*a + 2, 2*a^2 + 2*a + 2, 2*a^3 + 2*a^2 + 2*a + 2]
        
    It is also possible to request some elementary properties::
    
        sage: K.<a> = NumberField(x^3 - x^2 + 5*x + 1)
        sage: R = K.ring_of_integers()
        sage: I = K.factor(24)[0][0]^4 * K.factor(24)[2][0]
        sage: S = R.quotient(I)
        sage: S.cardinality()
        48
        sage: S.characteristic()
        12
        sage: S.is_finite()
        True
        sage: S.cover()
        Ring morphism:
          From: Maximal Order in Number Field in a with defining polynomial x^3 - x^2 + 5*x + 1
          To:   Quotient of Maximal Order in Number Field in a with defining polynomial x^3 - x^2 + 5*x + 1 by the ideal (-2*a + 2)
          Defn: Natural quotient map

    There is also a naturally defined lifting map, taking elements
    to their representative::
    
        sage: K.<a> = NumberField(x^2 - x + 22)
        sage: I = prod([f[0] for f in K.factor(21)])
        sage: S = K.ring_of_integers().quotient(I)
        sage: S.lifting_map()
        Set-theoretic ring morphism:
          From: Quotient of Maximal Order in Number Field in a with defining polynomial x^2 - x + 22 by the ideal (21, 7*a + 7)
          To:   Maximal Order in Number Field in a with defining polynomial x^2 - x + 22
          Defn: Choice of lifting map
        sage: S.lift(S(a))
        a

"""
#*****************************************************************************
#       Copyright (C) 2013 Joey van Langen <j.m.van.langen@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage
from sage.rings.quotient_ring import QuotientRing_generic
from sage.rings.quotient_ring_element import QuotientRingElement
from sage.rings.number_field.order import Order

from sage.rings.integer import Integer
from sage.rings.infinity import PlusInfinity

from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.fields import Fields

class OrderQuotientElement(QuotientRingElement):
    r"""
    An element of a quotient ring of an order.
    
    INPUT:

    - ``parent`` - the quotient ring of which this is an element.

    - ``rep`` - a representative of the element in the order; this is used
      as the internal representation of the element

    - ``reduce`` - bool (optional, default: True) - if True, then the
      internal representation of the element is ``rep`` reduced modulo
      the ideal `I`
    
    NOTE:
    
    This is simply an extension of the original QuotientRingElement class
    found in sage.rings.quotient_ring_element. This implementation changes
    the method :func:`lift` such that the method :func:`lifting_map` works
    for the parent, which must be an order.
    
    EXAMPLES::
    
        sage: K.<a> = NumberField(x^2 - 15)
        sage: I = K.prime_above(2) * K.prime_above(3) * K.prime_above(17)
        sage: S = K.maximal_order().quotient(I)
        sage: S(a^3 - 4)
        -4*a - 1
    
    An example testing the lift::
        
        sage: K.<a> = NumberField(x^3 - 4*x + 17)
        sage: R = K.ring_of_integers()
        sage: D = K.factor(30)
        sage: I = D[0][0] * D[3][0] * D[4][0]
        sage: S = R.quotient(I)
        sage: S(a^3 - 4)
        4*a - 1
        sage: a^3 - 4 - S(a^3 - 4).lift() in I
        True

    """
    
    def __init__(self, parent, rep, **kwds):
        r"""
        Initializes an element of a quotient ring of an order.
        
        INPUT:
        
        - ``parent`` - The quotient ring of which this is an element.
        
        - ``x`` - A representative of this element in the ambient ring
          of the quotient ring.
          
        - ``reduce`` - boolean (default: True), if ``True`` this element
          will be represented by an element of the ambient ring reduced
          modulo the defining ideal of the quotient ring.
          
        EXAMPLES:
        
            sage: K.<a> = QuadraticField(-13)
            sage: R = K.maximal_order()
            sage: S = R.quotient(a+13); S
            Quotient of Maximal Order in Number Field in a with defining polynomial x^2 + 13 by the ideal (a + 13)
            sage: S(a - 4)
            a - 4
        """
        
        if rep not in parent.ambient():
            raise ValueError("%s can not be viewed as an element of %s"%(
                                                          rep, parent.ambient()))
        QuotientRingElement.__init__(self, parent, rep, **kwds)
        
    def lift(self):
        r"""
        Returns a representative of this element in the ambient ring
        of the quotient ring of which this is an element.
        
        OUTPUT:
        
        An element of :func:`parent().ambient`.
        
        EXAMPLES::
        
            sage: K.<a> = CyclotomicField(7)
            sage: R = K.maximal_order()
            sage: I = R.ideal(a^3 - 6*a + 13)
            sage: S = R.quotient(I); S
            Quotient of Maximal Order in Cyclotomic Field of order 7 and degree 6 by the ideal (a^3 - 6*a + 13)
            sage: S(a).lift()
            1888566*a^5
            sage: S(a).lift() - a in I
            True
        """
        return self.parent().ambient()(self._QuotientRingElement__rep)

class OrderQuotientRing(QuotientRing_generic, UniqueRepresentation):
    r"""
    A quotient of an order by an ideal.
    
    EXAMPLES::
        
        sage: K.<a> = NumberField(x^2 - x - 1)
        sage: R = K.maximal_order()
        sage: I = K.prime_above(5)^3
        sage: S = R.quotient(I); S
        Quotient of Maximal Order in Number Field in a with defining polynomial x^2 - x - 1 by the ideal (-10*a + 5)
        sage: S.cardinality()
        125
        
    An example for an order with non-principal ideals::
        
        sage: K.<a> = NumberField(x^2 - x + 113)
        sage: I = K.prime_above(7) * K.prime_above(11)
        sage: S = K.maximal_order().quotient(I); S
        Quotient of Maximal Order in Number Field in a with defining polynomial x^2 - x + 113 by the ideal (77, a + 16)
        sage: S.characteristic().factor()
        7 * 11
        sage: S(a+3)^17
        -25*a
    """
    
    Element = OrderQuotientElement
    def __init__(self, order, ideal, names=None, category=None):
        r"""
        Initialization of the quotient ring.
        
        INPUT:

        -  ``order`` -- An order, i.e. an instance of
           sage.rings.number_field.Order

        -  ``ideal`` -- An ideal of the order in the first argument.

        - ``names`` -- (default None) a list of generator names.
        
        - ``category`` -- (default None) a specific category this
          object must belong to.
        
        NOTE:
            Most of the structure is left to the QuotientRing_generic
            structure.
            
        SEE_ALSO:
            sage.rings.quotient_ring.QuotientRing_generic
            
        EXAMPLES::
            
            sage: K.<a> = NumberField(x^3 - x^2 + 3*x + 6)
            sage: R = K.maximal_order()
            sage: I = R.ideal(17)
            sage: R.quotient(I)
            Quotient of Maximal Order in Number Field in a with defining polynomial x^3 - x^2 + 3*x + 6 by the ideal (17)
        """
        if not isinstance(order, Order):
            raise TypeError("%s is not an order"%(order))
        if ideal not in order.ideal_monoid():
            raise TypeError("%s is not an ideal of %s"%(ideal, order))
        if not ideal.is_zero() and ideal.is_maximal():
            category = category or Fields()
        QuotientRing_generic.__init__(self,
                                      order,
                                      ideal,
                                      names=names,
                                      category=category)
    
    def cardinality(self):
        r"""
        Return the cardinality of the underlying set.
        
        In this case that would be the norm of the underlying
        ideal or "+Infinity" if it is the zero ideal.
        
        OUTPUT:
            
        Either a positive integer of "+Infinity".
        
        EXAMPLES::
        
            sage: K = QuadraticField(-7)
            sage: R = K.maximal_order()
            sage: I = K.factor(12)[0][0] * K.factor(12)[2][0]
            sage: S = R.quotient(I)
            sage: S.cardinality()
            18
        
        A case in which we can see the degrees of the corresponding
        residue fields::
            
            sage: K.<a> = CyclotomicField(7)
            sage: I = K.prime_above(3) * K.prime_above(2)
            sage: S = K.maximal_order().quotient(I)
            sage: S.cardinality().factor()
            2^3 * 3^6

        """
        n = self.defining_ideal().norm()
        if n == 0:
            return PlusInfinity
        else:
            return Integer(n)
            
    def characteristic(self):
        r"""
        Return the characteristic of the quotient ring.
        
        This is a generator for the kernel of the
        natural map of ZZ into this ring.
        
        OUTPUT:
        
        The smallest positive integer in the ideal underlying
        the quotient ring or "0" if none.
        
        EXAMPLES::
        
            sage: K.<a> = NumberField(x^2 + 5*x + 10)
            sage: R = K.maximal_order()
            sage: S = R.quotient(6)
            sage: S.characteristic()
            6

        An example where the result is not directly obvious::
        
            sage: K.<a> = NumberField(x^3 - x^2 - 3*x + 1)
            sage: E = EllipticCurve(K, [a^2 - a - 2, a^2 - 2*a - 3, a^2 - 1, -16*a^2 - 19*a + 8, 49*a^2 + 57*a - 23])
            sage: I = E.conductor()
            sage: S = K.maximal_order().quotient(I)
            sage: S.characteristic()
            4
        """
        return self.defining_ideal().smallest_integer()
        
    def krull_dimension(self):
        r"""
        Return the krull dimension of this commutative ring.
        
        The krull dimension is the length of the longest ascending
        chain of prime ideals.
        
        OUTPUT:
        
        An integer.
        
        NOTE:
        
        As an order has krull dimension '1' so the krull dimension of
        a quotient is always '0' except if we mod out by the zero
        ideal, in which case it again is '1'.
        
        EXAMPLES::
            
            sage: K.<a> = NumberField(x^2 - x + 4)
            sage: R = K.maximal_order().quotient(3*a - 5)
            sage: R.krull_dimension()
            0
            
        It returns 1 if the ideal is the zero ideal. Note that in
        that case the quotient method of the order will just return
        the order itself, hence we have to be specific to test it.::
        
            sage: K = NumberField(x^2 - x + 1,"a")
            sage: R = K.maximal_order()
            sage: I = R.zero_ideal()
            sage: S = sage.rings.number_field.order_quotients.OrderQuotientRing(R, I)
            sage: S.krull_dimension()
            1
        """
        if self.defining_ideal().is_zero():
            return 1
        else:
            return 0
            
    def order(self):
        """
        The number of elements of ``self``.
        
        This returns the same result as :func:`cardinality`.
        
        OUTPUT:
        
        An integer or `+Infinity`.
        
        .. SEEALSO::

            :func:`cardinality`

        EXAMPLES::
        
            sage: K.<a> = NumberField(x^2 + 5*x + 10)
            sage: R = K.maximal_order()
            sage: S = R.quotient(6)
            sage: S.order()
            36

        An example where the result is not directly obvious::
        
            sage: K.<a> = NumberField(x^3 - x^2 - 3*x + 1)
            sage: E = EllipticCurve(K, [a^2 - a - 2, a^2 - 2*a - 3, a^2 - 1, -16*a^2 - 19*a + 8, 49*a^2 + 57*a - 23])
            sage: I = E.conductor()
            sage: S = K.maximal_order().quotient(I)
            sage: S.order()
            16
        """
        return self.cardinality()
        
    def is_finite(self):
        r"""
        Return ``True`` if this ring is finite
        and ``False`` otherwise.
        
        NOTE:
        
        This will always be the case unless the ideal is
        the zero ideal, which must be explicitly
        constructed.
        
        EXAMPLES::
        
            sage: K = NumberField(x^2 - 15,"a")
            sage: I = K.prime_above(2) * K.prime_above(5)
            sage: S = K.maximal_order().quotient(I)
            sage: S.is_finite()
            True
            
        An example where the result is ``False`` must be
        explicitly constructed::
        
            sage: K = QuadraticField(17)
            sage: R = K.maximal_order()
            sage: I = R.zero_ideal()
            sage: from sage.rings.number_field.order_quotients import OrderQuotientRing
            sage: S = OrderQuotientRing(R, I)
            sage: S.is_finite()
            False
        """
        return not self.defining_ideal().is_zero()
        
    def __iter__(self):
        for a in self.defining_ideal().residues():
            yield self(a)
