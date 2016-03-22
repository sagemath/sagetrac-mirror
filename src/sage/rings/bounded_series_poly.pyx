#############################################################################
#    Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

include "../ext/stdsage.pxi"

from infinity import Infinity
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.rings.padics.precision_error import PrecisionError

from sage.misc.misc import newton_method_sizes

from bounded_series_ring_element cimport BoundedSeries
from bounded_series_ring_element import is_BoundedSeries


cdef class BoundedSeries_poly(BoundedSeries):
    def __init__(self, parent, f=0, valuation_final_terms=None, prec=None, prec_coeffs=None, check=True, is_gen=0):
        R = parent._polynomial_ring
        base = parent.base_ring()
        self._f = R(f, check=check)

        if prec is None:
            if isinstance(f, BoundedSeries):
                prec = f._prec
            else:
                prec = Infinity
        if prec is not Infinity:
            if check:
                self._f = self._f.truncate(prec)
        elif valuation_final_terms is not None:
            valuation_final_terms = Infinity

        if prec_coeffs is not None:
            coeffs = self._f.list()
            zero = base(0)
            for i in range(len(prec_coeffs)):
                if i < len(coeffs):
                    coeffs.append(zero)
                if prec_coeffs[i] is not Infinity:
                    coeffs[i] = coeffs[i].add_bigoh(prec_coeffs[i].ceil())
            self._f = R(coeffs)

        BoundedSeries.__init__(self, parent, None, valuation_final_terms, prec, is_gen)


    def __hash__(self):
        return hash(self._f)

    #def __reduce__(self):

    cpdef BoundedSeries _new_constant_series(self,RingElement a,Parent P,char check=0):
        cdef BoundedSeries_poly f = self.__new__(self.__class__)
        f._parent = P
        f._f = P._polynomial_ring(a)
        f._prec = Infinity
        f._valuation_final_terms = Infinity
        f._degree = 0
        return f

    def polynomial(self):
        return self._f

    cdef _compute_valuation_degree(self):
        log_radius = self._parent.log_radius()
        val = Infinity; deg = -1;
        f = self._f
        if f == 0:
            self._valuation = Infinity
            self._degree = -1
            self._mu = None
            return
        df = f.degree()
        val = f[df].valuation(); deg = df
        mu = Infinity
        for d in range(df-1, -1, -1):
            v = f[d].valuation()
            slope = (v - val) / (deg - d)
            if slope <= log_radius:
                val = v
                deg = d
                mu = Infinity
            elif slope < mu:
                mu = slope
        self._valuation = val + log_radius*deg
        self._degree = deg
        self._mu = mu

    def is_secure(self):
        if self._valuation is None:
            self._compute_valuation_degree()
        return self._valuation_final_terms is not None and self._valuation <= self._valuation_final_terms

    def secure_log_radius(self):
        if self._valuation_final_terms is None:
            return Infinity
        log_radius = self._parent._log_radius
        if self._valuation is not None and self._valuation <= self._valuation_final_terms:
            return log_radius
        f = self._f
        df = f.degree()
        prec = self._prec
        secure_log_radius = Infinity
        valuation_final_terms = self._valuation_final_terms - log_radius*prec
        for i in range(df+1):
            c = f[i]
            if c.is_zero(): continue
            val = c.valuation()
            slope = (val - valuation_final_terms) / (prec - i)
            if slope < secure_log_radius:
                secure_log_radius = slope
        return max(log_radius, secure_log_radius)

    def gauss_valuation(self, secure=False):
        if self._valuation is None:
            self._compute_valuation_degree()
        if secure and not self.is_secure():
            raise PrecisionError("Unable to determine for sure the Gauss valuation")
        return self._valuation

    def minimal_gauss_valuation(self):
        if self._valuation_final_terms is None:
            raise ValueError("This element is not secure")
        if self._valuation is None:
            self._compute_valuation_degree()
        return min(self._valuation, self._valuation_final_terms)

    def weierstrass_degree(self, secure=False):
        if self._valuation is None or self._degree is None:
            self._compute_valuation_degree()
        if secure and not self.is_secure():
            raise PrecisionError("Unable to determine for sure the Weierstrass degree")
        return self._degree

    def highest_slope(self, secure=False):
        if self._mu is None:
            self._compute_valuation_degree()
        if secure and not self.is_secure():
            raise ValueError("Unable to determine for sure the highest slope")
        return -self._mu

    def is_zero(self, secure=False):
        true_zero = True
        l = self._f.list()
        for c in l:
            if not c.is_zero():
                return False
            elif secure:
                prec = c.precision_absolute()
                if prec is not Infinity: true_zero = False
        if secure and not true_zero:
            raise PrecisionError("impossible to distinguish this element from zero")
        return True

    def __getitem__(self, n):
        valuation_final_terms = self._valuation_final_terms
        if n >= self._prec and valuation_final_terms is not Infinity:
            if valuation_final_terms is None:
                raise IndexError("coefficient not known")
            else:
                return self.base_ring()(0).add_bigoh(valuation_final_terms + n*self._parent._log_radius)
        return self._f[n]

    def __neg__(self):
        return BoundedSeries_poly(self._parent, -self._f, self._prec, check=False)
        
    cpdef ModuleElement _add_(self, ModuleElement right_m):
        cdef BoundedSeries_poly right = <BoundedSeries_poly>right_m
        if self._valuation_final_terms is None or right._valuation_final_terms is None:
            valuation_final_terms = None
        else:
            valuation_final_terms = min(self._valuation_final_terms, right._valuation_final_terms)
        res = self._f + right._f
        prec = min(self._prec, right._prec)
        return BoundedSeries_poly(self._parent, res, valuation_final_terms, prec)
                                         
    cpdef ModuleElement _sub_(self, ModuleElement right_m):
        cdef BoundedSeries_poly right = <BoundedSeries_poly>right_m
        if self._valuation_final_terms is None or right._valuation_final_terms is None:
            valuation_final_terms = None
        else:
            valuation_final_terms = min(self._valuation_final_terms + right.gauss_valuation(), right._valuation_final_terms)
        res = self._f - right._f
        prec = min(self._prec, right._prec)
        return BoundedSeries_poly(self._parent, res, valuation_final_terms, prec)

    cpdef RingElement _mul_(self, RingElement right_m):
        cdef BoundedSeries_poly right = <BoundedSeries_poly>right_m
        if self._valuation_final_terms is None or right._valuation_final_terms is None:
            valuation_final_terms = None
        else:
            valuation_final_terms = min(self._valuation_final_terms + right.minimal_gauss_valuation(), right._valuation_final_terms + self.minimal_gauss_valuation())
        res = self._f * right._f
        prec = min(self._prec, right._prec)
        return BoundedSeries_poly(self._parent, res, valuation_final_terms, prec)
                                         
    cpdef ModuleElement _rmul_(self, RingElement c):
        if self._valuation_final_terms is None:
            valuation_final_terms = None
        else:
            valuation_final_terms = self._valuation_final_terms + c.valuation()
        res = c * self._f
        return BoundedSeries_poly(self._parent, res, valuation_final_terms, self._prec)

    cpdef ModuleElement _lmul_(self, RingElement c):
        return self._rmul_(c)

    def truncate(self, prec=Infinity):
        if prec is Infinity:
            return self._f
        else:
            return self._f.truncate(prec)

    def list(self):
        return self._f.list()

    def inverse(self, secure=False):
        if not self.is_unit(secure=secure):
            raise ZeroDivisionError("This series in not invertible in the ring of bounded convergent series on {val >= %s}" % self.parent().log_radius())
        prec = self._prec
        if prec is Infinity:
            prec = self._parent.default_prec()
        f = self._f
        res = f.parent()(~f[0])
        for next_prec in newton_method_sizes(prec)[1:]:
            z = res.square() * f.truncate(next_prec)
            res = 2*res - z.truncate(next_prec)
        if self.is_secure():
            valuation_final_terms = -self._valuation
        else:
            valuation_final_terms = None
        return BoundedSeries_poly(self._parent, res, valuation_final_terms, prec)

    def weierstrass_preparation(self, monic=False, secure=False, compute_inverse=False):
        if secure and not self.is_secure():
            raise PrecisionError("This series is not secure")
        d = self.weierstrass_degree()

        if self._a is None or self._b is None:

            f = self._f
            prec = self._prec
            log_radius = self._parent._log_radius

            # Compute loss of precision coming from self._valuation_final_terms
            ###################################################################

            mu = self._mu

            if prec is Infinity or mu is Infinity:
                precs_res = [ ]
            else:
                lastprec = self._valuation_final_terms
                if lastprec is None:
                    lastprec = self._valuation
                lastprec -= prec * log_radius
                precs_res = [ lastprec + i*mu for i in range(prec) ]
                precs_res.reverse()
            for i in range(d, len(precs_res)):
                precs_res[i] -= self._valuation

            # Compute the decomposition
            ###########################

            parent_poly = f.parent()
            if d == 0:
                a = parent_poly(1)
                b = f
                v = parent_poly(~f[0])
            else:
                a = f._factor_of_degree(d)
                b = f // a

            parent = self._parent
            if self._valuation_final_terms is None or self._valuation > self._valuation_final_terms:
                self._a = BoundedSeries_poly(parent, a, valuation_final_terms=None, prec=Infinity, prec_coeffs=precs_res[:d])
            else:
                self._a = BoundedSeries_poly(parent, a, valuation_final_terms=Infinity, prec=Infinity, prec_coeffs=precs_res[:d])
            valb = self._valuation + d*log_radius
            if self._valuation_final_terms is None or self._valuation > self._valuation_final_terms:
                valuation_final_terms = None
            else:
                valuation_final_terms = valb
            self._b = BoundedSeries_poly(parent, b, valuation_final_terms=valuation_final_terms, prec=prec-d, prec_coeffs=precs_res[d:])

        # Compute the inverse of b (if asked)
        #####################################

        if compute_inverse and self._binv is None:
            self._binv = self._b.inverse()

        # Return the result
        ###################

        a = self._a
        b = self._b
        binv = self._binv
        if monic:
            lc = a[d]
            a /= lc
            b *= lc
            if compute_inverse:
                binv /= lc
                return a, b, binv
            else:
                return a, b
        else:
            if compute_inverse:
                return a, b, binv
            else:
                return a, b


    def quo_rem(self, right, secure=False):
        if right.parent() != self.parent():
            raise TypeError("Parents do not coincide")
        prec = self._prec
        log_radius = self._parent._log_radius
        valright = right.gauss_valuation(secure=secure)
        d = right.weierstrass_degree()

        # Compute loss of precision coming from self._valuation_final_terms
        ###################################################################

        mu = -right.highest_slope()
        if prec is Infinity or mu is Infinity:
            precs_res = [ ]
        else:
            lastprec = self._valuation_final_terms
            if lastprec is None:
                lastprec = self._valuation
            lastprec -= prec * log_radius
            precs_res = [ lastprec + i*mu for i in range(prec) ]
            precs_res.reverse()
        for i in range(d, len(precs_res)):
            precs_res[i] -= valright

        # Compute the result
        ####################

        parent = self._parent
        a, _, binv = right.weierstrass_preparation(compute_inverse=True)
        q, r = self._f.quo_rem(a.polynomial())
        if self._valuation_final_terms is None or binv.valuation_final_terms() is None:
            q = BoundedSeries_poly(parent, q, valuation_final_terms=None, prec=Infinity)
            r = BoundedSeries_poly(parent, r, valuation_final_terms=None, prec=Infinity, prec_coeffs=precs_res[:d])
        else:
            q = BoundedSeries_poly(parent, q, valuation_final_terms=Infinity, prec=Infinity)
            r = BoundedSeries_poly(parent, r, valuation_final_terms=Infinity, prec=Infinity, prec_coeffs=precs_res[:d])
        q *= binv
        q = q.add_bigoh(precs_res[d:])

        return q, r
