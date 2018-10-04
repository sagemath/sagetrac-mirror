"""
Ideals and fractional ideals of discrete valuation rings.

AUTHORS:

- David Roe (initial version)
"""

#*****************************************************************************
#       Copyright (C) 2011 David Roe <roed.math@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.infinity import infinity
from sage.rings.ideal import Ideal_pid, Ideal_generic, Ideal_fractional
import sage.misc.latex as latex

class Ideal_local_generic(Ideal_generic):
    """
    An ideal or fractional ideal in a discrete valuation ring.

    EXAMPLES::

        sage: Zp(5).ideal(25)
        Principal ideal (5^2) of 5-adic Ring with capped relative precision 20
    """
    def __init__(self, ring, gen, coerce=True, construct=False):
        """
        Initialization.

        INPUT:

        - ring -- a p-adic ring or field
        - gen -- an element or ring, or a list or tuple of elements of ring
        - coerce -- boolean (default True).  Whether to convert the entries of gen into ring.
        - construct -- boolean (default False).  If True, assumes gen is an integer or 
          infinity giving the valuation of this ideal.

        EXAMPLES::

            sage: Zp(3).ideal(36)
            Principal ideal (3^2) of 3-adic Ring with capped relative precision 20
        """
        if construct:
            self._v = gen
            gen = ring.uniformizer_pow(gen)
        else:
            if isinstance(gen, (list, tuple)):
                self._v = infinity
                izero = False
                for g in gen:
                    if coerce:
                        g = ring(g)
                    v = g.valuation()
                    if v < self._v:
                        self._v = v
                        if g._is_inexact_zero():
                            izero = True
                        else:
                            izero = False
                    elif izero and v == self._v and (not g._is_inexact_zero()):
                        izero = False
                if izero:
                    raise ValueError, "Given generators do not determine a well defined ideal"
                gen = ring.uniformizer_pow(self._v)
            else:
                if coerce:
                    gen = ring(gen)
                self._v = gen.valuation()
                if gen._is_inexact_zero():
                    raise ValueError, "Given generator does not determine a well defined ideal"
        Ideal_generic.__init__(self, ring, gen, coerce=False)

    def valuation(self):
        """
        The valuation of this ideal, which is equal to the valuation of any generator 
        (for the case of principal ideals).

        EXAMPLES::
        
            sage: Zp(5).ideal(25).valuation()
            2
        """
        return self._v

    def __cmp__(self, other):
        """
        Comparison.

        EXAMPLES::

            sage: I1 = Zp(5).ideal(5)
            sage: I2 = Zp(5).ideal(10)
            sage: I3 = Zp(5).ideal(25)
            sage: I1 == I2
            True
            sage: I1.gen() == I2.gen()
            True
            sage: I1 < I3
            True
        """
        if isinstance(other, Ideal_local_generic) and self.ring() is other.ring():
            return cmp(self._v, other._v)
        return Ideal_pid.__cmp__(self, other)

    def _repr_short(self):
        """
        Called by __repr__ in Ideal_generic.

        EXAMPLES::
        
            sage: Zp(5).ideal(17) # indirect doctest
            Principal ideal (1) of 5-adic Ring with capped relative precision 20
        """
        if self._v is infinity:
            return "(0)"
        elif self._v == 0:
            return "(1)"
        elif self._v == 1:
            return "(%s)"%(self.ring()._printer._ram_name())
        else:
            return "(%s^%s)"%(self.ring()._printer._ram_name(), self._v)

    def _contains_(self, x):
        """
        Tests containment for an element of the ring containing this ideal.

        EXAMPLES::

            sage: 20 in Zp(5).ideal(15) # indirect doctest
            True
        """
        return x.valuation() >= self._v

    def __nonzero__(self):
        """
        Returns True if this ideal is not the zero ideal.

        EXAMPLES::

            sage: bool(Zp(5).ideal(0))
            False
            sage: bool(Zp(5).ideal(1))  
            True
        """
        return self._v is not infinity

    def _latex_(self):
        """
        Latex representation of this ideal.

        EXAMPLES::

            sage: latex(Zp(5).ideal(50))
            \left(5^{2}\right)\ZZ_{5}
            sage: latex(Zp(5).ideal(17))
            \ZZ_{5}
        """
        if self._v is infinity:
            return "(0)"
        elif self._v == 0:
            return latex.latex(self.ring())
        elif self._v == 1:
            return "\\left(%s\\right)%s"%(self.ring()._printer._ram_name(), latex.latex(self.ring()))
        else:
            return "\\left(%s^{%s}\\right)%s"%(self.ring()._printer._ram_name(), self._v, latex.latex(self.ring()))

    def reduce(self, f):
        """
        Reduce an element of the associated ring to a chosen representative.

        EXAMPLES::

            sage: Zp(5).ideal(75).reduce(44)
            4 + 3*5 + O(5^20)
        """
        f = self.ring()(f)
        if self._v is infinity:
            return f
        else:
            return f % self.ring().uniformizer_pow(self._v) # not the most efficient method

    def is_maximal(self):
        """
        Returns whether this ideal is the maximal ideal.

        EXAMPLES::

            sage: Zp(5).ideal(7).is_maximal()
            False
            sage: Zp(5).ideal(40).is_maximal()
            True
        """
        return self._v == 1

    def is_prime(self):
        """
        Returns whether this ideal is prime.

        EXAMPLES::

            sage: Zp(5).ideal(0).is_prime()
            True
        """
        return self._v == 1 or self._v is infinity

    def is_trivial(self):
        """
        Returns whether this ideal is (1) or (0)

        EXAMPLES::

            sage: Zp(5).ideal(17).is_trivial()
            True
            sage: Zp(5).ideal(0).is_trivial()
            True
            sage: Zp(5).ideal(35).is_trivial()
            False
        """
        return self._v == 0 or self._v is infinity
    
    def __add__(self, other):
        """
        Returns self or other, whichever has smaller valuation.

        EXAMPLES::

            sage: I = Zp(5).ideal(35)
            sage: J = Zp(5).ideal(225)
            sage: I + J # indirect doctest
            Principal ideal (5) of 5-adic Ring with capped relative precision 20
        """
        if not (isinstance(other, Ideal_local_generic) and other.ring() is self.ring()):
            other = self.ring().ideal(other)
        if self._v <= other._v:
            return self
        else:
            return other

    def __radd__(self, other):
        """
        Returns self or other, whichever has smaller valuation.

        EXAMPLES::

            sage: I = Zp(5).ideal(35)
            sage: J = Zp(5).ideal(225)
            sage: J.__radd__(I)
            Principal ideal (5) of 5-adic Ring with capped relative precision 20
        """
        if not (isinstance(other, Ideal_local_generic) and other.ring() is self.ring()):
            other = self.ring().ideal(other)
        if self._v <= other._v:
            return self
        else:
            return other

    def gcd(self, other):
        """
        Returns self or other, whichever has smaller valuation.

        EXAMPLES::

            sage: I = Zp(5).ideal(95)
            sage: J = Zp(5).ideal(75)
            sage: I.gcd(J)
            Principal ideal (5) of 5-adic Ring with capped relative precision 20
        """
        if not (isinstance(other, Ideal_local_generic) and other.ring() is self.ring()):
            other = self.ring().ideal(other)
        if self._v <= other._v:
            return self
        else:
            return other

    def __mul__(self, other):
        """
        Returns the product.

        EXAMPLES::

            sage: I = Zp(5).ideal(105)
            sage: J = Zp(5).ideal(50)
            sage: I * J # indirect doctest
            Principal ideal (5^3) of 5-adic Ring with capped relative precision 20
        """
        if not (isinstance(other, Ideal_local_generic) and other.ring() is self.ring()):
            other = self.ring().ideal(other)
        return self.__class__(self.ring(), self._v + other._v, construct=True)

    def __rmul__(self, other):
        """
        Returns the product.

        EXAMPLES::

            sage: I = Zp(5).ideal(105)
            sage: J = Zp(5).ideal(50)
            sage: J.__rmul__(I)
            Principal ideal (5^3) of 5-adic Ring with capped relative precision 20
        """
        if not (isinstance(other, Ideal_local_generic) and other.ring() is self.ring()):
            other = self.ring().ideal(other)
        return self.__class__(self.ring(), self._v + other._v, construct=True)

    def __invert__(self):
        """
        Returns the inverse of this ideal as a fractional ideal.

        EXAMPLES::

            sage: ~Zp(5).ideal(375)
            Fractional ideal (5^-3) of 5-adic Field with capped relative precision 20
        """
        return Ideal_local_fractional(self.ring().fraction_field(), -self._v, construct=True)

class Ideal_local(Ideal_local_generic, Ideal_pid):
    """
    An actual ideal in a DVR.

    EXAMPLES::

        sage: Zp(5).ideal(10)
        Principal ideal (5) of 5-adic Ring with capped relative precision 20
    """
    pass

class Ideal_local_fractional(Ideal_local_generic, Ideal_fractional):
    """
    A fractional ideal in the fraction field of a DVR.

    EXAMPLES::

        sage: K = Qp(5)
        sage: ~K.fractional_ideal(75)
        Fractional ideal (5^-2) of 5-adic Field with capped relative precision 20
    """
    pass
