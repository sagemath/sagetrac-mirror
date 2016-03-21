#############################################################################
#    Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation 
from sage.structure.parent_gens import normalize_names

import sage.categories.basic as categories

from integer import Integer
from rational_field import QQ 
from sage.structure.element import Element 
from ring import Field, CommutativeAlgebra
from infinity import Infinity 
import sage.misc.latex as latex

from polynomial.polynomial_ring_constructor import PolynomialRing 
from polynomial.polynomial_element import Polynomial 
from sage.structure.parent_gens import ParentWithGens

from bounded_series_ring_element import BoundedSeries 
from bounded_series_ring_morphism import BoundedSeriesBaseringInjection
from bounded_series_ring_morphism import BoundedSeriesRestriction


def BoundedSeriesRing(base_ring, log_radius=0, name=None, names=None,
                      sparse=False, default_prec=None):
    if default_prec is None:
        default_prec = 20

    if name is None:
        name = names
    if name is None:
        raise TypeError, "You must specify the name of the indeterminate of the Bounded Power series ring."
    try:
        name = normalize_names(1, name)[0]
    except IndexError:
        raise NotImplementedError("Multivariate bounded power series rings are not implemented yet.")
    except TypeError:
        raise TypeError, "illegal variable name"
    
    if not (isinstance(base_ring, Field)): # and base_ring.is_cdvf()):  <-- TODO: fix this
        raise TypeError("base_ring must be a complete discrete valuation field")
    return BoundedSeriesRing_generic(base_ring, log_radius, name, default_prec, sparse=sparse)


def is_BoundedSeriesRing(R):
    return isinstance(R, BoundedSeriesRing_generic)


class BoundedSeriesRing_generic(CommutativeAlgebra, UniqueRepresentation):
    @staticmethod
    def __classcall__(cls, base_ring, log_radius=0, name=None, default_prec=20, sparse=False, element_class=None):
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse bounded series are not implemented")
            else:
                from bounded_series_poly import BoundedSeries_poly
                element_class = BoundedSeries_poly
        return super(BoundedSeriesRing_generic,cls).__classcall__(cls, base_ring, log_radius, name, default_prec, sparse, element_class)
        
    def __init__(self, base_ring, log_radius=0, name=None, default_prec=20, sparse=False, element_class=None, category=None):
        self.__is_sparse = sparse
        self._series_class = element_class
        self._by_one = False
        self._log_radius = log_radius
        self._default_prec = default_prec
        self._polynomial_ring = PolynomialRing(base_ring, name=name)

        # Algebra.__init__ also calls __init_extra__ of Algebras(...).parent_class, which tries to provide a conversion from the base ring, if it does not exist. This is for algebras that only do the 
        # generic stuff in their initialisation. But here, we want to use PolynomialBaseringInjection. Hence, we need to wipe the memory and construct the conversion from scratch.
        CommutativeAlgebra.__init__(self, base_ring, names=name, normalize=True, category=category)
        self.__generator = self._series_class(self, self._polynomial_ring.gen(), valuation_final_terms=Infinity, prec=Infinity, is_gen=True)
        self._base_inject = BoundedSeriesBaseringInjection(base_ring,self)
        self._coercions_log_radius = { }

        if log_radius in QQ:
            self._refine_category_(categories.EuclideanDomains())

    #def __reduce__(self):

    def has_coerce_map_from(self, S):
        base = self.base_ring()
        if base.has_coerce_map_from(S):
            return True
        elif isinstance(S, BoundedSeriesRing_generic):
            return self.base_ring() is S.base_ring() and self.variable_name() is S.variable_name() and self._log_radius >= S.log_radius()

    def coerce_map_from(self, S):
        base = self.base_ring()
        if base is S:
            return self._base_inject
        elif base.has_coerce_map_from(S):
            return self._base_inject * base.coerce_map_from(S) 
        elif isinstance(S, BoundedSeriesRing_generic):
            log_radius = S.log_radius()
            if self.base_ring() is S.base_ring() and self.variable_name() is S.variable_name() and self._log_radius >= log_radius:
                key = (self.default_prec(), log_radius)
                try:
                    return self._coercions_log_radius[key]
                except KeyError:
                    map = BoundedSeriesRestriction(S, self, check=False)
                    self._coercions_log_radius[key] = map
                    return map

    def _element_constructor_(self, x=None, valuation_final_terms=Infinity, prec=Infinity, check=True, is_gen=False):
        if is_gen:
            x = [ 0,1 ]
            prec = Infinity
            valuation_final_terms = Infinity
        elif isinstance(x, Element):
            if isinstance(x, BoundedSeries):
                if x.log_radius() > self.log_radius():
                    raise TypeError("Impossible to decrease the logarithmic radius of convergence")
                prec = x.prec()
                valuation_final_terms = x.valuation_final_terms()
                if valuation_final_terms is not None:
                    valuation_final_terms += prec * (self.log_radius() - x.log_radius())
                x = x.polynomial()
            #elif isinstance(x, PowerSeries):
            if isinstance(x, Polynomial):
                pass
            else:
                x = [ x ]
        elif isinstance(x, int):
            x = [ x ]
        return self._series_class(self, x, valuation_final_terms, prec)

    def _repr_(self):
        s = "Bounded Convergent Series Ring in %s on {val > %s} over %s"%(self.variable_name(), self.log_radius(), self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    def is_sparse(self):
        return self.__is_sparse

    def is_dense(self):
        return not self.__is_sparse

    def log_radius(self):
        return self._log_radius

    def default_prec(self):
        return self._default_prec

    def _latex_(self):
        return "%s\left<%s\right>"%(latex.latex(self.base_ring()), self.latex_variable_names()[0])

    #def _is_valid_homomorphism_(self, codomain, im_gens):

    def base_extend(self, R):
        raise NotImplementedError
    
    def change_ring(self, R):
        raise NotImplementedError

    def change_log_radius(self, log_radius):
        return self.__class__(self.base_ring(), log_radius, self.variable_name(), sparse=self.is_sparse(), default_prec=self.default_prec())

    def change_default_prec(self, default_prec):
        return self.__class__(self.base_ring(), self.log_radius(), self.variable_name(), sparse=self.is_sparse(), default_prec=default_prec)

    def is_exact(self):
        return False

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "generator n>0 not defined"
        return self.__generator

    def ngens(self):
        return 1

    def random_element(self, prec=None, degree=None, valuation=0, *args, **kwds):
        from sage.functions.other import ceil
        if prec is None:
            prec = self.default_prec()
        if degree is None or degree >= prec:
            degree = prec - 1
        if degree is Infinity:
            degree = self.default_prec() - 1
        val = valuation
        log_radius = self._log_radius
        coeffs = [ ]
        base = self.base_ring()
        integer_base = base.integer_ring()
        for i in range(degree+1):
            coeffs.append(base(integer_base.random_element()) << ceil(val))
            val -= log_radius
        if prec is Infinity:
            valuation_final_terms = Infinity
        else:
            valuation_final_terms = valuation
        return self._series_class(self, coeffs, valuation_final_terms=valuation_final_terms, prec=prec)

    def __cmp__(self, other):
        return self is other

    def uniformizer(self):
        if self._log_radius not in QQ:
            raise ValueError("This ring does not have a uniformizer")
        base = self.base_ring()
        unif = base.uniformizer()
        image_val = unif.valuation()
        fraction = self._log_radius / image_val
        num = fraction.numerator()
        denom = fraction.denominator()
        _, n, v = num.xgcd(denom)
        if n < 0:
            n += denom
            v -= num
        coeffs = n * [ base(0) ] + [ unif ** v ]
        return self._series_class(self, coeffs)

    def is_atomic_repr(self):
        return False

    def is_commutative(self):
        return True

    def is_field(self, proof = True):
        return False

    def is_finite(self):
        return False

    def characteristic(self):
        return self.base_ring().characteristic()

    #def integers(self):
    #    return BoundedByOneSeriesRing_generic(self.base_ring(), self._log_radius(), self.variable_name(), self._default_prec, self.__is_sparse, self._series_class)

    def hom(self, im_gen, morphism_on_coefficients=None):
        from sage.categories.homset import Hom
        from bounded_series_ring_morphism import BoundedSeriesHomomorphism_im_gens
        if isinstance(im_gen, (int, Integer)):
            parent = self
        elif isinstance(im_gen, list):
            parent = im_gen[0].parent()
        else:
            parent = im_gen.parent()
        homset = Hom(self, parent)
        return BoundedSeriesHomomorphism_im_gens(homset, im_gen, morphism_on_coefficients)

    def frobenius_endomorphism(self, n=1, b=None):
        from sage.categories.homset import Hom
        from bounded_series_ring_morphism import BoundedSeriesHomomorphism_im_gens
        base = self.base_ring()
        if b is None:
            b = base.characteristic()
            if b == 0: b = base.prime()
        morphism_on_coefficients = base.frobenius_endomorphism(n=n)
        homset = Hom(self, self)
        return BoundedSeriesHomomorphism_im_gens(homset, b**n, morphism_on_coefficients)
