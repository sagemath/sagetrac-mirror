#############################################################################
#    Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

include "../ext/stdsage.pxi"

import sage

from infinity import Infinity
from integer import Integer

from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields
from sage.rings.morphism cimport RingHomomorphism 
from sage.rings.morphism cimport RingHomomorphism_im_gens 
from sage.categories.homset import Hom



cdef class BoundedSeriesBaseringInjection(RingHomomorphism):
    def __init__(self, domain, codomain):
        assert codomain.base_ring() is domain, "domain must be basering"
        RingHomomorphism.__init__(self, Hom(domain,codomain))
        self._an_element = codomain.gen()
        self._repr_type_str = "Bounded Series base injection"
        self._new_constant_series_ = self._an_element._new_constant_series

    cpdef Element _call_(self, x):
        return self._new_constant_series_(x, self._codomain)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        try:
            return self._codomain._element_constructor_(x, *args, **kwds)
        except AttributeError:
            # if there is no element constructor, there is a custom call method.
            return self._codomain(x, *args, **kwds)


cdef class BoundedSeriesRestriction(RingHomomorphism):
    def __init__(self, domain, codomain, check=True):
        self._diff_log_radius = codomain.log_radius() - domain.log_radius()
        if check:
            from bounded_series_ring import BoundedSeriesRing_generic
            if not isinstance(domain, BoundedSeriesRing_generic):
                raise TypeError("The domain must be a ring of Bounded Series")
            if not isinstance(codomain, BoundedSeriesRing_generic):
                raise TypeError("The codomain must be a ring of Bounded Series")
            if domain.base_ring() is not codomain.base_ring():
                raise TypeError("The domain and the codomain must share the same base ring")
            if self._diff_log_radius < 0:
                raise TypeError("Log radius of convergence must increase")
        RingHomomorphism.__init__(self, Hom(domain,codomain))
        self._repr_type_str = "Restriction morphism:\n From: %s\n To: %s" % (domain, codomain)
        self._series_class = codomain._series_class

    def _repr_(self):
        return self._repr_type_str

    cpdef Element _call_(self, x):
        prec = x.prec()
        valuation_final_terms = x.valuation_final_terms()
        if prec is not Infinity and valuation_final_terms is not None:
            valuation_final_terms += self._diff_log_radius * prec
        return self._series_class(self.codomain(), x.polynomial(), valuation_final_terms, prec)


cdef class BoundedSeriesHomomorphism_im_gens(RingHomomorphism_im_gens):
    def __init__(self, parent, im_gens, morphism_on_coefficients=None, check=True):
        RingHomomorphism.__init__(self, parent)
        domain = self.domain(); codomain = self.codomain()
        from bounded_series_ring import BoundedSeriesRing_generic

        if check:
            if not isinstance(domain, BoundedSeriesRing_generic):
                raise TypeError("the domain must be a Bounded Series Ring")

        if isinstance(im_gens, (int, Integer)):
            image = im_gens
            if check and image < 1:
                raise ValueError("When im_gens is an integer, it must be positive")
            self._zeroes = (image - 1) * [ codomain.base_ring()(0) ]
            im_gens = [ codomain.gen() ** im_gens ]
        elif not isinstance(im_gens, sage.structure.sequence.Sequence_generic):
            if not isinstance(im_gens, (tuple, list)):
                im_gens = [im_gens]
            im_gens = sage.structure.all.Sequence(im_gens, codomain)
            if check and len(im_gens) != 1:
                raise ValueError("too many images")
            image = im_gens[0]

        base = domain.base_ring()
        unif = base.uniformizer()
        if codomain in CompleteDiscreteValuationFields():
            e = codomain(unif).valuation() / unif.valuation()
            self._function = self._call_cdvf
            self._gain_precision = image.valuation() - e*domain.log_radius()
            if check:
                if self._gain_precision <= 0:
                    raise ValueError("Image does not define a valid morphism")
                base_codomain = codomain.base_ring()
                if (morphism_on_coefficients is not None and 
                   (not isinstance(morphism_on_coefficients, RingHomomorphism) or morphism_on_coefficients.domain() != base or morphism_on_coefficients.codomain() != base)):
                    raise ValueError("morphism_on_coefficients must be an endomorphism of the base ring of the domain (= %s)" % base)
        elif isinstance(codomain, BoundedSeriesRing_generic):
            base_codomain = codomain.base_ring()
            if morphism_on_coefficients is None:
                e = base_codomain(unif).valuation() / unif.valuation()
            else:
                e = morphism_on_coefficients(unif).valuation() / unif.valuation()
            if self._zeroes is None:
                self._gain_precision = image.valuation(secure=True) - e*domain.log_radius()
                self._function = self._call_composition
                self._vertices = image.minimal_newton_polygon().vertices()
            else:
                self._gain_precision = image*codomain.log_radius() - e*domain.log_radius()
                self._function = self._call_composition_monomial
            if check:
                if self._gain_precision < 0 or (self._zeroes is None and self._gain_precision == 0 and image.weierstrass_degree() == 0):
                    raise ValueError("Image does not define a valid morphism")
                if (morphism_on_coefficients is not None and 
                   (not isinstance(morphism_on_coefficients, RingHomomorphism) or morphism_on_coefficients.domain() != base or morphism_on_coefficients.codomain() != base_codomain)):
                    raise ValueError("morphism_on_coefficients must be a morphism from the base ring of the domain (= %s) to the base ring of the codomain (= %s)" % (base, base_codomain))
        else:
            raise NotImplementedError

        self._e = e
        self._image = image
        self.__im_gens = im_gens
        if morphism_on_coefficients is not None and not morphism_on_coefficients.is_identity():
            self._morphism = morphism_on_coefficients

    def _repr_defn(self):
        s = "%s |--> %s" % (self._domain.variable_name(), self.__im_gens[0])
        if self._morphism is not None:
            s += "\nAction on coefficients: %s" % (self._morphism._repr_short())
        return s

    def morphism_on_coefficients(self):
        return self._morphism

    cpdef Element _call_(self, x):
        return self._function(x)

    def _call_cdvf(self, x):
        image = self._image
        codomain = self._codomain
        pow = codomain(1)
        res = 0
        morphism = self._morphism
        if morphism is None:
            morphism = codomain.base_ring()
        for c in x.list():
            res += morphism(c)*pow
            pow *= image
        valuation_final_terms = x.valuation_final_terms()
        if valuation_final_terms is None:
            return res
        prec = self._e * valuation_final_terms + x.prec() * self._gain_precision 
        if prec is not Infinity:
            from sage.functions.other import ceil
            res = res.add_bigoh(ceil(prec))
        return res

    def _call_composition(self, f):
        from sage.functions.other import ceil
        image = self._image
        codomain = self._codomain
        base_codomain = codomain.base_ring()
        prec = f.prec()
        vertices = self._vertices
        d = vertices[-1][0]
        trunc = d*prec
        pow = codomain(1)
        series = codomain(0)
        morphism = self._morphism
        if morphism is None:
            morphism = base_codomain
        for c in f.list():
            series += morphism(c)*pow
            pow *= image
        valuation_final_terms = f.valuation_final_terms()
        coeffs = series.list()
        if prec is not Infinity:
            if len(coeffs) < trunc:
                coeffs.extend((trunc-len(coeffs))*[base_codomain(0)])
            e = self._e
            (ax,ay) = vertices[0]
            for (x,y) in vertices[1:]:
                slope = e * (ay-y) / (ax-x)
                bigoh = e * prec * ay
                start = prec*ax
                end = prec*x
                for i in range(start, min(end,trunc)):
                    coeffs[i] = coeffs[i].add_bigoh(ceil(bigoh))
                    bigoh += slope
                if end > trunc: break
                (ax,ay) = (x,y)
            if valuation_final_terms is not None:
                valuation_final_terms = e * valuation_final_terms + prec * self._gain_precision
        valuation_final_terms = min(valuation_final_terms, series.valuation_final_terms())
        prec = min(trunc, series.prec())
        return codomain._series_class(codomain, coeffs, valuation_final_terms, prec, check=False)


    def _call_composition_monomial(self, x):
        codomain = self._codomain
        base_codomain = codomain.base_ring()
        zeroes = self._zeroes
        ans = [ ]
        coeffs = x.list()
        morphism = self._morphism
        if morphism is None:
            for c in coeffs:
                ans.append(base_codomain(c))
                ans.extend(zeroes)
        else:
            for c in coeffs:
                ans.append(morphism(c))
                ans.extend(zeroes)
        prec = x.prec()
        valuation_final_terms = x.valuation_final_terms()
        if valuation_final_terms is not None:
            valuation_final_terms = self._e * valuation_final_terms + prec * self._gain_precision
        return codomain._series_class(codomain, ans, valuation_final_terms, self._image * prec, check=False)

    def _composition(self, right_m):
        cdef BoundedSeriesHomomorphism_im_gens right
        if isinstance(right_m, BoundedSeriesHomomorphism_im_gens):
            right = <BoundedSeriesHomomorphism_im_gens>right_m
            if self._morphism is None:
                morphism_on_coefficients = right._morphism
            elif right._morphism is None:
                morphism_on_coefficients = self._morphism
            else:
                morphism_on_coefficients = self._morphism * right._morphism
            if self._zeroes is not None and right._zeroes is not None:
                im_gens = self._image * right._image
            else:
                im_gens = self._function(right._image)
            return BoundedSeriesHomomorphism_im_gens(self.parent(), im_gens, morphism_on_coefficients)
        else:
            return RingHomomorphism._composition(self, right_m)

    cpdef int _cmp_(self, Element other) except -2:
        c = cmp(self.domain(), other.domain())
        if c: return c
        c = cmp(self.codomain(), other.codomain())
        if c: return c
        try:
            c = cmp(self.morphism_on_coefficients(), other.morphism_on_coefficients())
        except AttributeError:
            c = (self._morphism is None)
        if c: return c
        c = cmp(self.im_gens(), other.im_gens())
        return c

    def __hash__(self):
        return hash((self.domain(), self.codomain(), (self._morphism, self.__im_gens[0])))

