#############################################################################
#    Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

include "../ext/stdsage.pxi"

import operator, sage

from infinity import Infinity
from integer import Integer
from rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.structure.element cimport Element, AlgebraElement

from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields

from sage.rings.padics.precision_error import PrecisionError


def is_BoundedSeries(x):
    return isinstance(x, BoundedSeries)

cdef class BoundedSeries(AlgebraElement):
    def __init__(self, parent, f, valuation_final_terms, prec, is_gen=False):  # NB: f is ignored here
        AlgebraElement.__init__(self, parent)
        self.__is_gen = is_gen
        if not (prec is Infinity):
            prec = int(prec)
        self._prec = prec
        self._valuation_final_terms = valuation_final_terms

    #def __hash__(self):
        
    #def __reduce__(self):

    def is_sparse(self):
        return self._parent.is_sparse()

    def is_dense(self):
        return self._parent.is_dense()

    def is_gen(self):
        return self.__is_gen

    def base_extend(self, R):
        S = self._parent.base_extend(R)
        return S(self)

    def change_ring(self, R):
        S = self._parent.change_ring(R)
        return S(self)
    
    cpdef int _cmp_(self, Element right) except -2:
        prec = min(self.prec(), right.prec())
        if prec is Infinity:
            x = self.list()
            y = right.list()
        else:
            x = self.list()[:prec]
            y = right.list()[:prec]
        zero = self.base_ring()(0)
        dl = len(x) - len(y)
        if dl < 0:
            x += (-dl) * [ zero ]
        elif dl > 0:
            y += dl * [ zero ]
        for i in range(len(x)):
            c = cmp(x[i], y[i])
            if c: return c
        return 0
    
    def __call__(self, x): # you *MUST* override this in the derived class
        raise NotImplementedError

    def coefficients(self):
        zero = self.parent().base_ring()(0)
        return [c for c in self.list() if c != zero]

    def exponents(self):
        zero = self.parent().base_ring()(0)
        l = self.list()
        return [i for i in range(len(l)) if l[i] != zero]
    
    def list(self): # you *MUST* override this in the derived class
        raise NotImplementedError

    def polynomial(self): # you *MUST* override this in the derived class
        raise NotImplementedError

    def __copy__(self):
        return self

    def base_ring(self):
        return self._parent.base_ring()

    def variable(self):
        return self._parent.variable_name()

    def prec(self):
        return self._prec

    def log_radius(self):
        if self._prec is Infinity:
            return Infinity
        else:
            return self.parent().log_radius()

    def change_prec(self, new_prec=None, valuation_final_terms=None):
        f = self.polynomial()
        if new_prec is None:
            new_prec = self.parent().default_prec()
        if self._valuation_final_terms is None:
            valuation_final_terms = None
        if new_prec is not Infinity:
            if valuation_final_terms is None:
                if self._valuation_final_terms is Infinity:
                    valuation_final_terms = 0
                else:
                    valuation_final_terms = self._valuation_final_terms
            else:
                valuation_final_terms = min(valuation_final_terms, self._valuation_final_terms)
            log_radius = self._parent._log_radius
            for i in range(new_prec, f.degree()+1):
                valuation_final_terms = min(valuation_final_terms, f[i].valuation() + i*log_radius)
        return self.__class__(self._parent, f, valuation_final_terms, new_prec)

    def add_bigoh(self, new_precs, shift=0):
        if new_precs is Infinity:
            return self
        f = self.polynomial()
        valuation_final_terms = self._valuation_final_terms
        log_radius = self._parent._log_radius
        prec = self._prec
        if isinstance(new_precs, list):
            new_precs = [ Infinity ] * shift + new_precs
            if prec is not Infinity and valuation_final_terms is not None:
                for i in range(prec, len(new_precs)):
                    valuation_final_terms = min(valuation_final_terms, new_precs[i] - i*log_radius)
        else:
            valuation_final_terms = min(valuation_final_terms, new_precs)
            if prec is Infinity:
                prec = max(f.degree() + 1, self.parent().default_prec())
            new_precs = [ Infinity ] * shift + [ new_precs - i*log_radius for i in range(shift, prec) ]
        return self.__class__(self._parent, f, valuation_final_terms, prec, new_precs)

    def weierstrass_degree(self, secure=False):
        log_radius = self._parent.log_radius()
        val = Infinity; deg = -1;
        coeffs = self.list()
        for i in range(len(coeffs)):
            v = coeffs[i].valuation() - i*log_radius
            if v < val:
                val = v
                deg = i
        if secure and val > self._valuation_final_terms:
            raise PrecisionError("Unable to determine for sure the Weierstrass degree")
        return deg

    def degree(self, secure=False):
        return self.weierstrass_degree(secure=secure)

    def valuation_final_terms(self):
        return self._valuation_final_terms

    def gauss_valuation(self, secure=False):
        log_radius = self._parent.log_radius()
        val = Infinity
        coeffs = self.list()
        for i in range(len(coeffs)):
            v = coeffs[i].valuation() - i*log_radius
            if v < val: val = v
        if secure and val > self._valuation_final_terms:
            raise PrecisionError("Unable to determine for sure the Gauss valuation")
        return val

    def valuation(self, secure=False):
        return self.gauss_valuation(secure=secure)

    def _repr_(self):
        from sage.rings.integer_ring import ZZ
        X = self._parent.variable_name()
        if self._valuation_final_terms is None:
            bigoh = "O(unknown)"
        elif self._prec is Infinity or self._valuation_final_terms is Infinity:
            bigoh = ""
        else:
            log_radius = self.parent().log_radius()
            unif = self.base_ring().uniformizer()
            power_unif = (self._valuation_final_terms - log_radius * self._prec) / unif.valuation()
            try:
                unif = self.base_ring().variable_name()
            except AttributeError:
                pass
            bigoh = "O("
            if power_unif == 0:
                pass
            elif power_unif == 1:
                bigoh += "%s*" % unif
            elif power_unif in ZZ and power_unif > 0:
                bigoh += "%s^%s*" % (unif, power_unif)
            else:
                bigoh += "%s^(%s)*" % (unif, power_unif)
            bigoh += "%s^%s)" % (X, self._prec)

        if self.is_zero():
            return bigoh if bigoh != "" else "0"

        try:
            atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        except KeyError:
            atomic_repr = False

        s = " "
        v = self.list()
        m = len(v)
        first = True
        for n in xrange(m):
            x = v[n]
            if not x.is_zero():
                x = repr(x)
                if not first:
                    s += " + "
                if not atomic_repr and n > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(X,n)
                elif n==1:
                    var = "*%s"%X
                else:
                    var = ""
                s += "%s%s"%(x,var)
                first = False

        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if bigoh != "": s += " + %s"%bigoh
        return s[1:]

    #def _latex_(self):

    def __getitem__(self,n):
        if n<0:
            return self.base_ring()(0)
        c = self.list()
        if n >= len(c):
            if self._prec > n:
                return self.base_ring()(0)
            else:
                raise IndexError, "coefficient not known"
        return c[n]

    def is_zero(self, secure=False):
        raise NotImplementedError
    
    def __nonzero__(self):
        return not self.is_zero()

    def is_secure(self):
        return self._valuation_final_terms is not None and self.valuation() <= self._valuation_final_terms

    def secure_log_radius(self):
        if self._valuation_final_terms is None:
            return Infinity
        log_radius = self._parent._log_radius
        coeffs = self.list()
        prec = self._prec
        secure_log_radius = Infinity
        valuation_final_terms = self._valuation_final_terms - log_radius*prec
        for i in range(len(coeffs)):
            c = coeffs[i]
            if c.is_zero(): continue
            val = c.valuation()
            slope = (val - valuation_final_terms) / (prec - i)
            if slope < secure_log_radius:
                secure_log_radius = slope
        return max(log_radius, secure_log_radius)

    def change_log_radius(self, log_radius):
        if log_radius < self.parent().log_radius():
            raise TypeError("log radius can't decrease")
        R = self.parent().change_log_radius(log_radius)
        return R(self)

    def restriction(self, log_radius):
        return self.change_log_radius(log_radius)

    def is_unit(self, secure=False):
        return self.weierstrass_degree(secure=secure) == 0

    def newton_polygon(self,secure=True):
        from sage.geometry.newton_polygon import NewtonPolygon
        if secure and not self.is_secure():
            raise PrecisionError("The Newton polygon is not determined")
        coeffs = self.list()
        lastslope = -self._parent.log_radius()
        vertices = [ ]
        for x in range(len(coeffs)):
            c = coeffs[x]
            if c.is_zero(): continue
            vertices.append((x, c.valuation()))
        polygon = NewtonPolygon(vertices, last_slope=lastslope)
        if secure:
            vertices_prec = [ (x, coeffs[x].precision_absolute()) for x in range(len(coeffs)) ]
            if self._prec is not Infinity:
                vertices_prec.append((self._prec, self._valuation_final_terms - self._prec * self._parent._log_radius))
            polygon_prec = NewtonPolygon(vertices_prec, last_slope=lastslope)
            vertices = polygon.vertices()
            vertices_prec = polygon_prec.vertices()
            if vertices[0][0] < vertices_prec[0][0]:
                raise PrecisionError("The Newton polygon is not determined")
            for (x,y) in vertices:
                if polygon_prec(x) <= y:
                    raise PrecisionError("The Newton polygon is not determined")
            (x, y) = vertices_prec[-1]
            if polygon(x) > y:
                raise PrecisionError("The Newton polygon is not determined")
        return polygon

    def minimal_newton_polygon(self):
        from sage.geometry.newton_polygon import NewtonPolygon
        coeffs = self.list()
        lastslope = -self._parent.log_radius()
        vertices = [ (x, coeffs[x].valuation()) for x in range(len(coeffs)) ]
        prec = self.prec()
        valuation_final_terms = self._valuation_final_terms
        if valuation_final_terms is not None and prec is not Infinity:
            vertices.append((prec, valuation_final_terms - prec*self.log_radius()))
        return NewtonPolygon(vertices, last_slope=lastslope)

    def newton_slopes(self,repetition=True,secure=True):
        polygon = self.newton_polygon(secure=secure)
        return [ -s for s in polygon.slopes(repetition=repetition) ]

    def __invert__(self):
        raise NotImplementedError

    def inverse(self, secure=False):
        raise NotImplementedError

    def inverse_of_unit(self, secure=False):
        return self.inverse(self, secure=secure)

    def weierstrass_preparation(self, monic=False, secure=False):
        raise NotImplementedError

    def quo_rem(self, right, secure=False):
        raise NotImplementedError

    def __floordiv__(self,other):
        q,r = self.quo_rem(other)
        return q
            
    def __mod__(self,other):
        q,r = self.quo_rem(other)
        return r

    def __call__(self, x, **kwds):
        from bounded_series_ring import BoundedSeriesRing_generic
        parent_self = self.parent()
        parent_x = x.parent()
        if isinstance(parent_x, BoundedSeriesRing_generic):
            if self._valuation_final_terms is None:
                raise NotImplementedError("Composition of series involving O(unknown) not implemented")
            # First, we consider the case where self is actually a polynomial
            if self._prec is Infinity:
                return self.polynomial()(x)
            # We compute the radius of convergence
            prec = x.prec()
            unif = parent_self.base_ring().uniformizer()
            e = parent_x.base_ring()(unif).valuation() / unif.valuation()
            v = e * parent_self.log_radius()
            coeffs = x.list()
            log_radius = -coeffs[0].valuation()
            if log_radius >= v:
                raise ValueError("The composite does not converge")
            for i in range(1, len(coeffs)):
                slope = (v - coeffs[i].valuation()) / i
                if log_radius < slope:
                    log_radius = slope
            if prec is not Infinity:
                valuation_final_terms = x.valuation_final_terms()
                slope = (v - valuation_final_terms)/prec - x.log_radius()
                if log_radius < slope:
                    log_radius = slope
            parent_x = parent_x.change_log_radius(log_radius)
            x = x.change_log_radius(log_radius)
            check = False
        else:
            check = True
        from bounded_series_ring_morphism import BoundedSeriesHomomorphism_im_gens
        from sage.categories.homset import Hom
        morphism = BoundedSeriesHomomorphism_im_gens(Hom(parent_self, parent_x), x, check=check)
        return morphism(self, **kwds)

    def derivative(self):
        return self.__class__(self._parent, self.polynomial().derivative(), self.valuation_final_terms()+self.log_radius(), max(0,self._prec-1), check=False)
    
    def __setitem__(self, n, value):
        raise IndexError, "power series are immutable"
