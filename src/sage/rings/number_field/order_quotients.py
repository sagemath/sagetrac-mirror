r"""
Quotient rings for orders inside number fields.

Implements functionality for quotients of orders inside number fields,
that does not hold for general quotients. The classes in this modulo
build upopn the classes in sage.rings.quotient_Ring and sage.rings.quotient_ring_element.

AUTHORS:

- Joey van Langen (2017-06-19): initial version

EXAMPLES:



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
from sage.structure.unique_representation import UniqueRepresentation

class OrderQuotientElement(QuotientRingElement):

    def __init__(self, parent, x):
        if x not in parent.ambient():
            raise ValueError("%s can not be viewed as an element of %s"%(
                                                          x, parent.ambient()))
        QuotientRingElement.__init__(self, parent, x)
        
    def lift(self):
        return self.parent().ambient()(self._QuotientRingElement__rep)

class OrderQuotientRing(QuotientRing_generic, UniqueRepresentation):
    
    Element = OrderQuotientElement
    def __init__(self, order, ideal, names=None, category=None):
        if not isinstance(order, Order):
            raise TypeError("%s is not an order"%(order))
        if ideal not in order.ideal_monoid():
            raise TypeError("%s is not an ideal of %s"%(ideal, order))
        if ideal.is_maximal():
            category = category or Fields()
        QuotientRing_generic.__init__(self,
                                      order,
                                      ideal,
                                      names=names,
                                      category=category)
    
    def cardinality(self):
        n = self.defining_ideal().norm()
        if n == 0:
            return +Infinity
        else:
            return Integer(n)
            
    def characteristic(self):
        return self.defining_ideal.smallest_integer()
        
    def krull_dimension(self):
        if self.is_zero:
            return 0
        else:
            return 1
            
    def order(self):
        return self.cardinality()
        
    def is_finite(self):
        return self.cardinality() != +Infinity
        
    def __iter__(self):
        for a in self.defining_ideal().residues():
            yield self(a)
