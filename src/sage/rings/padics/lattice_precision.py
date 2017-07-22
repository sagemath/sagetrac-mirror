"""
Toy implementation of lattice precision for p-adic numbers.

Here is a small demo::

    sage: R = ZpLP(3)
    sage: x = R(1,10)

Of course, when we multiply by 3, we gain one digit of absolute 
precision::

    sage: 3*x
    3 + O(3^11)

The lattice precision machinery sees this even if we decompose
the computation into several steps::

    sage: y = x+x
    sage: y
    2 + O(3^10)
    sage: x+y
    3 + O(3^11)

The same works for the multiplication::

    sage: z = x^2
    sage: z
    1 + O(3^10)
    sage: x*z
    1 + O(3^11)

This comes more funny when we are working with elements given
at different precisions::

    sage: R = ZpLP(2)
    sage: x = R(1,10)
    sage: y = R(1,5)
    sage: z = x+y; z
    2 + O(2^5)
    sage: t = x-y; t
    0 + O(2^5)
    sage: z+t  # observe that z+t = 2*x
    2 + O(2^11)
    sage: z-t  # observe that z-t = 2*y
    2 + O(2^6)

    sage: x = R(28888,15)
    sage: y = R(204,10)
    sage: z = x/y; z
    242 + O(2^9)
    sage: z*y  # which is x
    28888 + O(2^15)

The SOMOS sequence is the sequence defined by the recurrence::

..MATH::

    u_n = \frac {u_{n-1} u_{n-3} + u_{n-2}^2} {u_{n-4}}

It is known for its numerical instability.
On the one hand, one can show that if the initial values are
invertible in `\mathbb{Z}_p` and known at precision `O(p^N)`
then all the next terms of the SOMOS sequence will be known 
at the same precision as well.
On the other hand, because of the division, when we unroll
the recurrence, we loose a lot of precision. Observe::

    sage: R = Zp(2, print_mode='terse')
    sage: a,b,c,d = R(1,15), R(1,15), R(1,15), R(3,15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    4 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    13 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    55 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    21975 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    6639 + O(2^13)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    7186 + O(2^13)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    569 + O(2^13)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    253 + O(2^13)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    4149 + O(2^13)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    2899 + O(2^12)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    3072 + O(2^12)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    349 + O(2^12)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    619 + O(2^12)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    243 + O(2^12)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    3 + O(2^2)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    2 + O(2^2)

If instead, we use the lattice precision, everything goes well::

    sage: R = ZpLP(2)
    sage: a,b,c,d = R(1,15), R(1,15), R(1,15), R(3,15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    4 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    13 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    55 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    21975 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    23023 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    31762 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    16953 + O(2^15)
    sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print d
    16637 + O(2^15)

    sage: for _ in range(100):
    ....:     a,b,c,d = b,c,d,(b*d+c*c)/a
    sage: a
    15519 + O(2^15)
    sage: b
    32042 + O(2^15)
    sage: c
    17769 + O(2^15)
    sage: d
    20949 + O(2^15)
"""

import _weakref as weakref
from sage.misc.prandom import randint
from sage.rings.integer_ring import ZZ
from sage.rings.padics.generic_nodes import pAdicRingBaseGeneric
from sage.rings.padics.padic_generic_element import pAdicGenericElement

class pAdicRingLattice(pAdicRingBaseGeneric):
    # Internal variables:
    #  . self._prec_cap
    #    a cap for the (working) precision
    #    meaning that the precision lattice always contains p^(self._prec_cap)
    #  . self._elements
    #    list of weak references of elements in this parent
    #  . self._precision
    #    an upper triangular matrix over ZZ representing the 
    #    lattice of precision
    #    (its columns are indexed by self._elements)

    def __init__(self, p, prec=50, print_mode={'mode':'terse'}, label=None):
        if label is None:
            self._label = None
        else:
            self._label = str(label)
        self._elements = [ ]
        self._precision = { }
        self._prec_cap = prec
        self._modulus = p**prec
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, '', pAdicLatticeElement)

    def _prec_type(self):
        return 'lattice'

    def label(self):
        return self._label

    def precision_cap(self):
        return self._prec_cap

    def _repr_(self):
        s = "%s-adic Field with lattice precision" % self.prime()
        if self._label is not None:
            s += ", called '%s'" % self._label
        return s

    def _new_element(self, x, prec, dx):
        n = len(self._elements)
        x_ref = weakref.ref(x, self._del_element)
        self._elements.append(x_ref)
        if prec is None or prec > self._prec_cap:
            prec = self._prec_cap
        modulus = self.prime()**prec
        col = n * [ZZ(0)]
        for ref,scalar in dx:
            c = self._precision[ref]
            for i in range(len(c)):
                col[i] += scalar * c[i]
        for i in range(n):
            col[i] = ZZ(col[i]) % modulus
            # We should handle the case where col[i] is not an integer but I'm lazy...
        col.append(modulus)
        self._precision[x_ref] = col
        
    def _del_element(self, ref):
        try:
            index = len(self._precision[ref]) - 1
        except IndexError:
            return
        del self._elements[index]
        del self._precision[ref]

        # Now, we echelonize...
        p = self.prime()
        n = len(self._elements)
        modulus = self._modulus
        for i in range(index,n):
            col = self._precision[self._elements[i]]
            d, u, v = col[i].xgcd(col[i+1])
            up, vp = col[i+1]/d, -col[i]/d
            col[i] = d
            del col[i+1]
            for j in range(i+1,n):
                col = self._precision[self._elements[j]]
                col[i], col[i+1] = (u*col[i] + v*col[i+1]) % modulus, (up*col[i] + vp*col[i+1]) % modulus
        # Actually, the matrix we obtained this way is not in Hermite normal form
        # because its entries above the diagonal are not necessarily completely reduced
        # (though they are at least modulo p^precision_cap)
        # Should we reduce them further (this affects seriously the complexity)?

    def precision_absolute(self, x):
        ref = weakref.ref(x)
        p = self.prime()
        return min([ c.valuation(p) for c in self._precision[ref] ])

    def working_precision(self, x):
        # Do something better here
        return self._prec_cap

ZpLP = pAdicRingLattice


# The elements
##############

class pAdicLatticeElement(pAdicGenericElement):
    # Internal variable:
    #  . self._approximation
    #    an approximation of this p-adic number
    #    it is defined modulo p^(working_precision)

    def __init__(self, parent, x, prec=None, dx={}):
        self._parent = parent
        self._hash = hash((x, randint(1, 100000)))
        pAdicGenericElement.__init__(self, parent)
        if prec is None:
            prec = parent.precision_cap()
        self._parent._new_element(self, prec, dx)
        self._approximation = x

    def __hash__(self):
        return self._hash

    def working_precision(self):
        return self.parent().working_precision(self)

    def approximation(self, reduce=True):
        if reduce:
            workprec = self.working_precision()
            p = self._parent.prime()
            self._approximation %= p**workprec
        return self._approximation

    def precision_absolute(self):
        return self.parent().precision_absolute(self)

    def valuation(self, secure=False):
        p = self._parent.prime()
        val = self._approximation.valuation(p)
        prec = self.precision_absolute()
        if val < prec: 
            return val
        elif secure:
            raise PrecisionError("Not enough precision")
        else:
            return prec

    def precision_relative(self, secure=False):
        return self.precision_absolute() - self.valuation(secure=secure)

    def _repr_(self):
        prec = self.precision_absolute()
        p = self._parent.prime()
        app = self._approximation % (p**prec)
        return "%s + O(%s^%s)" % (app, p, prec)

    def _add_(self, other):
        x = self.approximation() + other.approximation()
        dx = [  [weakref.ref(self), 1], 
               [weakref.ref(other), 1] ]
        return self.__class__(self._parent, x, dx=dx)

    def _sub_(self, other):
        x = self.approximation() - other.approximation()
        dx = [  [weakref.ref(self), 1 ], 
               [weakref.ref(other), -1] ]
        return self.__class__(self._parent, x, dx=dx)

    def _mul_(self, other):
        x_self = self.approximation()
        x_other = other.approximation()
        x = x_self * x_other
        dx = [  [weakref.ref(self), x_other],
               [weakref.ref(other), x_self ] ]
        return self.__class__(self._parent, x, dx=dx)

    def _div_(self, other):
        p = self._parent.prime()
        x_self = self.approximation()
        x_other = other.approximation()
        wp_other = other.working_precision()
        d, inv, _ = x_other.xgcd(p**wp_other)
        q, r = x_self.quo_rem(d)
        if r != 0:
            raise ValueError("division is not exact")        
        x = q * inv
        # dx = (1/other)*dself - (self/other^2)*dother
        dx = [  [weakref.ref(self), inv/d    ],
               [weakref.ref(other), -x*inv/d ] ]
        return self.__class__(self._parent, x, dx=dx)
