"""
Relative extensions of `p`-adic rings

We represent general extensions of p-adic rings as a two-step extension:
first an unramified extension of Qp, followed by an Eisenstein extension
of the result.

This file contains the parent classes for such extensions.
"""

#*****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from .generic_nodes import pAdicFixedModRingGeneric, pAdicCappedAbsoluteRingGeneric, pAdicCappedRelativeRingGeneric, pAdicCappedRelativeFieldGeneric, pAdicFloatingPointRingGeneric, pAdicFloatingPointFieldGeneric
from .eisenstein_extension_generic import EisensteinExtensionGeneric
from .relative_ramified_FM import RelativeRamifiedFixedModElement
from .relative_ramified_CA import RelativeRamifiedCappedAbsoluteElement
from .relative_ramified_CR import RelativeRamifiedCappedRelativeElement
from .relative_ramified_FP import RelativeRamifiedFloatingPointElement
from .pow_computer_relative import PowComputer_relative_maker

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.rings.padics.factory import Qq

class pAdicRelativeBaseringInjection(Morphism):
    """
    The injection of the unramified base into the two-step extension.

    INPUT:

    - ``R`` -- an unramified `p`-adic ring or field
    - ``S`` -- an eisenstein extension of ``R``.

    EXAMPLES::

        sage: K.<a> = Qq(125)
        sage: R.<x> = K[]
        sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
        sage: f = W.coerce_map_from(K); f
        Generic morphism:
          From: 5-adic Unramified Extension Field in a defined by x^3 + 3*x + 3
          To:   5-adic Eisenstein Extension Field in w defined by x^3 + 15*a*x - 5*a^2 - 5 over its base field
    """
    def __init__(self, R, S):
        """
        Initialization.

        EXAMPLES::

            sage: K.<a> = Qq(125)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
            sage: f = W.coerce_map_from(K)
            sage: type(f)
            <class 'sage.rings.padics.relative_extension_leaves.pAdicRelativeBaseringInjection'>
        """
        if not R.is_field() or S.is_field():
            Morphism.__init__(self, Hom(R, S))
        else:
            from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
            Morphism.__init__(self, Hom(R, S, SetsWithPartialMaps()))

    def _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: K.<a> = Qq(125,2)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
            sage: f = W.coerce_map_from(K)
            sage: f(a+5)
            a + (4*a^2 + 4*a + 3)*w^3 + (a + 2)*w^4 + (2*a^2 + 4*a + 2)*w^5 + O(w^6)
        """
        if x.is_zero():
            return self.codomain()(0,x.precision_absolute())
        else:
            return self.codomain()([x])

    def _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        EXAMPLES::

            sage: K.<a> = Qq(125,2)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
            sage: f = W.coerce_map_from(K)
            sage: f(5*a,5)
            (4*a^2 + a + 3)*w^3 + (a^2 + 2*a)*w^4 + O(w^5)
            sage: f(5*a,8,2)
            (4*a^2 + a + 3)*w^3 + (a^2 + 2*a)*w^4 + O(w^5)
        """
        return self.codomain()([x], *args, **kwds)

    def section(self):
        """
        Map back to the base ring.

        EXAMPLES::

            sage: K.<a> = Qq(125,2)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
            sage: f = W.coerce_map_from(K)
            sage: g = f.section()
            sage: g(a + w - w)
            a + O(5^2)
        """
        return pAdicRelativeBaseringSection(self.codomain(), self.domain())

class pAdicRelativeBaseringSection(Morphism):
    """
    The map from a two-step extension back to its maximal unramified subextension.

    EXAMPLES::

        sage: K.<a> = Qq(2^10)
        sage: R.<x> = K[]
        sage: W.<w> = K.extension(x^4 + 2*a*x^2 - 16*x - 6)
        sage: f = K.convert_map_from(W); f
        Generic morphism:
          From: 2-adic Eisenstein Extension Field in w defined by x^4 + 2*a*x^2 - 16*x - 6 over its base field
          To:   2-adic Unramified Extension Field in a defined by x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
    """
    def __init__(self, S, R):
        """
        Initialization.

        EXAMPLES::

            sage: K.<a> = Qq(2^10)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^4 + 2*a*x^2 - 16*x - 6*a)
            sage: f = K.convert_map_from(W); type(f)
            <class 'sage.rings.padics.relative_extension_leaves.pAdicRelativeBaseringSection'>
        """
        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
        Morphism.__init__(self, Hom(S, R, SetsWithPartialMaps()))

    def _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: K.<a> = Qq(2^10)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^4 + 2*a*x^2 - 16*x - 6*a)
            sage: f = K.convert_map_from(W)
            sage: f(a + w - w)
            a + O(2^20)
            sage: f(w)
            Traceback (most recent call last):
            ...
            ValueError: Element not contained in base ring
        """
        f = x.polynomial()
        if f.degree() > 0:
            raise ValueError("Element not contained in base ring")
        return f[0]

    def _call_with_args(self, x, args=(), kwds={}):
        """
        Used when specifying absolute or relative precision.

        EXAMPLES::

            sage: K.<a> = Qq(2^10)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^4 + 2*a*x^2 - 16*x - 6*a)
            sage: f = K.convert_map_from(W)
            sage: f(a, 5)
            a + O(2^5)
        """
        return self.codomain()(self._call_(x), *args, **kwds)

class RelativeRamifiedExtensionRingFixedMod(EisensteinExtensionGeneric, pAdicFixedModRingGeneric):
    """
    Two-step extension ring with fixed-mod precision.

    EXAMPLES::

        sage: A.<a> = ZqFM(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Ring in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base field
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: A.<a> = ZqFM(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16) # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(prec = unram_prec+1)
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'fixed-mod')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedFixedModElement)
        from .relative_ramified_FM import pAdicCoercion_ZZ_FM, pAdicConvert_QQ_FM
        self.register_coercion(pAdicCoercion_ZZ_FM(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_FM(self))

class RelativeRamifiedExtensionRingCappedAbsolute(EisensteinExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    """
    Two-step extension ring with capped absolute precision.

    EXAMPLES::

        sage: A.<a> = ZqCA(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Ring in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base field
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: A.<a> = ZqCA(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16) # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(show_prec=False, type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-abs')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedCappedAbsoluteElement)
        from .relative_ramified_CA import pAdicCoercion_ZZ_CA, pAdicConvert_QQ_CA
        self.register_coercion(pAdicCoercion_ZZ_CA(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_CA(self))

class RelativeRamifiedExtensionRingCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeRingGeneric):
    """
    Two-step extension ring with capped relative precision.

    EXAMPLES::

        sage: A.<a> = ZqCR(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Ring in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base field
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: A.<a> = ZqCR(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16) # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(field=False, show_prec=False, type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-rel')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedCappedRelativeElement)
        from .relative_ramified_CR import pAdicCoercion_ZZ_CR, pAdicConvert_QQ_CR
        self.register_coercion(pAdicCoercion_ZZ_CR(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_CR(self))

class RelativeRamifiedExtensionFieldCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    """
    Two-step extension field with capped relative precision.

    EXAMPLES::

        sage: A.<a> = QqCR(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Field in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base field
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: A.<a> = QqCR(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16) # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        #KFP = approx_modulus.base_ring().change(field=False, show_prec=False, type='floating-point')
        KFP = approx_modulus.base_ring().change(show_prec=False, type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, True, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-rel')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedCappedRelativeElement)
        from .relative_ramified_CR import pAdicCoercion_ZZ_CR, pAdicCoercion_QQ_CR
        self.register_coercion(pAdicCoercion_ZZ_CR(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        # We also want to convert down to the ring of integers: this is used in teichmuller expansion
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring().integer_ring(), self))
        self.register_coercion(pAdicCoercion_QQ_CR(self))

class RelativeRamifiedExtensionRingFloatingPoint(EisensteinExtensionGeneric, pAdicFloatingPointRingGeneric):
    """
    Two-step extension ring with floating point precision.

    EXAMPLES::

        sage: A.<a> = ZqFP(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Ring in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base field
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: A.<a> = ZqFP(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16) # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(field=False, show_prec=False)
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-rel')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedFloatingPointElement)
        from .relative_ramified_FP import pAdicCoercion_ZZ_FP, pAdicConvert_QQ_FP
        self.register_coercion(pAdicCoercion_ZZ_FP(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_FP(self))

class RelativeRamifiedExtensionFieldFloatingPoint(EisensteinExtensionGeneric, pAdicFloatingPointFieldGeneric):
    """
    Two-step extension field with floating point precision.

    EXAMPLES::

        sage: A.<a> = QqFP(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Field in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base field
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: A.<a> = QqFP(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16) # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(field=False, show_prec=False)
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, True, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-rel')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedFloatingPointElement)
        from .relative_ramified_FP import pAdicCoercion_ZZ_FP, pAdicCoercion_QQ_FP
        self.register_coercion(pAdicCoercion_ZZ_FP(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        # We also want to convert down to the ring of integers: this is used in teichmuller expansion
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring().integer_ring(), self))
        self.register_coercion(pAdicCoercion_QQ_FP(self))

class RelativeExtensionRingFixedMod(EisensteinExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        pass

class RelativeExtensionRingCappedAbsolute(EisensteinExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        pass

class RelativeExtensionRingCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeRingGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        pass

class RelativeExtensionRingFloatingPoint(EisensteinExtensionGeneric, pAdicFloatingPointRingGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        pass

class RelativeExtensionFieldFloatingPoint(EisensteinExtensionGeneric, pAdicFloatingPointFieldGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        pass

# Currently only for unramified extensions of general extensions
#    L
#  K1  K
#    K0
#    Qp
# K1/K0/Qp is all unramified
# K/K0 and L/K1 are Eisenstein, with same defining polynomial
# 
class RelativeExtensionFieldCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    """
    Three-step extension field with capped relative precision.

    EXAMPLES::

        sage: A.<a> = QqFP(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Field in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base field
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        self._approx_modulus = approx_modulus
        self._given_ground_ring = approx_modulus.base_ring() # this is K, should be an EisensteinExtension
        
        K0 = self._given_ground_ring.ground_ring()
        K = self._given_ground_ring
        
        section_to_unramified_base = pAdicRelativeBaseringSection(K, K0)

        # iterate through coefficients of approx_modulus and see if they lie in the unramified base
        try:
            unramified_extension_defining_poly = PolynomialRing(K0, name='t0')([section_to_unramified_base(coefficient) for coefficient in approx_modulus.list()])
        except:
            raise ValueError("Defining polynomial must be defined over the unramified extension.") # maybe everything works if I just comment out this line?
        
        # make K1, write down the element of K1 that maps to the generator of K0
        self._construct_maximal_unramified_subextension()

        eisenstein_extension_defining_poly = PolynomialRing(self.K1, name='t1')([coefficient.polynomial()(self._K0_gen) for coefficient in K.defining_polynomial().list()])

        unram_prec = prec
        KFP = self.K1.change(show_prec=False, type='floating-point')

        unif = K._exact_modulus.base_ring()(K0.uniformizer())
        shift_seed = PolynomialRing(self.K1, name='t1')([self._injection_from_K0_to_K1(coefficient) for coefficient in (-K._exact_modulus[:K._exact_modulus.degree()] / unif).change_ring(K0)])

        eisenstein_extension_defining_poly.change_ring(KFP)
        shift_seed.change_ring(KFP)

        self.prime_pow = PowComputer_relative_maker(self.K1.prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, True, eisenstein_extension_defining_poly.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-rel')

        self._implementation = 'Relative'
        EisensteinExtensionGeneric.__init__(self, eisenstein_extension_defining_poly, prec, print_mode, (names[0],names[1],self.K1.variable_names()[0],names[3]), RelativeRamifiedCappedRelativeElement)

        # self._construct_given_gen_v2()
        self._construct_given_gen(unramified_extension_defining_poly.list())
        self._construct_user_representation_of_K1_gen()
        self.register_coercion(pAdicUnramifiedOverGeneralBaseringInjection(self._given_ground_ring, self))
    
    def _construct_maximal_unramified_subextension(self):
        K0 = self._given_ground_ring.ground_ring()
        relative_degree = self._approx_modulus.degree()

        # construct K1 (the maximal unramified extension)
        from sage.rings.padics.factory import Qq
        # self.K1 = Qq(K0.residue_characteristic()**(relative_degree * K0.degree()), prec = K0.precision_cap(), type = 'capped-rel', names='a1', print_mode='terse')

        self.K1 = K0.change(q = K0.residue_characteristic()**(relative_degree * K0.degree()), names='a1')

        # find a way to map a generator of K0 into K1
        self._K0_gen = self._find_root([self.K1(coefficient) for coefficient in K0.defining_polynomial().list()], self.K1)

        # cache powers of the generator of K1 in terms of a K0-basis
        K0_gen_coefficients = self._K0_gen.polynomial()
        self._powers_of_K1_gen = [[K0(1)] + [K0(0)] * (relative_degree - 1)]
        for i in range(1, self.K1.degree()):
            next_term = [K0(0)] + self._powers_of_K1_gen[-1][:(relative_degree - 1)]
            for j in range(0, relative_degree):
                next_term[j] = next_term[j] - self._powers_of_K1_gen[-1][-1]*K0_gen_coefficients[j]
            self._powers_of_K1_gen.append(next_term)

    def _construct_user_representation_of_K1_gen(self):
        K0 = self._given_ground_ring.ground_ring()

        given_K1_gen = self._given_gen.polynomial()[0] # self.K1(self._given_gen.expansion()[0])
        relative_degree = self._approx_modulus.degree()

        K0_degree = self._given_ground_ring.ground_ring().degree()
        target_basis = [0] * (K0_degree * relative_degree)
        target_basis[0] = self.K1(1)

        for i in range(1, K0_degree):
            target_basis[i] = target_basis[i - 1] * self._K0_gen
        for j in range(1, relative_degree):
            for i in range(0, K0_degree):
                target_basis[j * K0_degree + i] = target_basis[(j - 1) * K0_degree + i] * given_K1_gen

        target_basis_vectors = []
        for basis_entry in target_basis:
            polynomial_list = basis_entry.polynomial().list()
            polynomial_list += [0] * (self.K1.degree() - len(polynomial_list))
            target_basis_vectors += polynomial_list

        # now target_basis[j * K0_degree + i] should equal
        # (self._K0_gen)^i * (given_K1_gen)^j 
        
        from sage.matrix.constructor import matrix
        target_basis_matrix = matrix(self._given_ground_ring.ground_ring().ground_ring(), self.K1.degree(), self.K1.degree(), target_basis_vectors)
        target_basis_matrix_inverse = target_basis_matrix.inverse()

        R0 = PolynomialRing(K0, name='t0')
        
        self._powers_of_K1_gen_in_given_basis = []
        for i in range(0, self.K1.degree()):
            target_basis_matrix_inverse_row = target_basis_matrix_inverse.row(i)

            self._powers_of_K1_gen_in_given_basis.append(R0([K0(list(target_basis_matrix_inverse_row[j*K0_degree:(j+1)*K0_degree])) for j in range(0, relative_degree)]))

    def _write_in_K0_basis_relative(self, element): # element is in K1
        if (element.parent() != self.K1):
            raise ValueError("Only elements of K1 can be expressed in the K0-basis.")
        element_in_Qp_basis = element.polynomial().list()
        return sum( element_in_Qp_basis[i] * self._powers_of_K1_gen_in_given_basis[i] for i in range(0, len(element_in_Qp_basis))).list()

    def _write_in_K0_basis(self, element): # element is in K1
        if (element.parent() != self.K1):
            raise ValueError("Only elements of K1 can be expressed in the K0-basis.")
        K0 = self._given_ground_ring.ground_ring() 
        relative_degree = self._approx_modulus.degree()

        element_in_Qp_basis = element.polynomial().list()

        return [ sum( element_in_Qp_basis[i] * self._powers_of_K1_gen[i][j] for i in range(0, len(element_in_Qp_basis))) for j in range(0, relative_degree) ]

    def _express_via_given_gen(self, element):
        element_in_K1_basis = element.polynomial().list()
        eisenstein_degree = len(element_in_K1_basis)
        element_in_K0_basis = [self._write_in_K0_basis_relative(entry) for entry in element_in_K1_basis]
        relative_degree = self._approx_modulus.degree()
        K = self._given_ground_ring

        return [K([element_in_K0_basis[j][i] for j in range(0, eisenstein_degree)]) for i in range(0, relative_degree)]

    def _injection_from_K0_to_K1(self, element):
        return element.polynomial()(self._K0_gen)

    def _section_from_K1_to_K0(self, element):
        element_in_K0_basis = self._write_in_K0_basis(element)        

        for j in range(1, len(element_in_K0_basis)):
            if element_in_K0_basis[j] != 0:
                raise ValueError("Element not contained in base field.")

        return element_in_K0_basis[0]
         
#     # if one day, we want to allow unramified extensions by polynomials whose coefficients are not necessarily in K0
#     def _construct_given_gen_v2(self):
#         self._given_gen = self._find_root([_injection_from_given_basering(coefficient) for coefficient in self._approx_modulus], self)
# 
    def _construct_given_gen(self, unramified_extension_defining_poly_coef_list):
        # elements of unramified_extension_defining_poly_coef_list are in K0, need to send them to K1
        coef_list = [self._injection_from_K0_to_K1(coef) for coef in unramified_extension_defining_poly_coef_list]
        # find root in K1
        root_in_K1 = self._find_root(coef_list, self.K1)
        # lift from K1 to L
        self._given_gen = self([root_in_K1])
    
    def _find_root(self, coeff_list, F):
        prec = F.precision_cap()
        coeff_residue_list = [ coefficient.residue() for coefficient in coeff_list ]
        root_guess = F(PolynomialRing(F.residue_field(), name='t1b')(coeff_residue_list).roots()[0][0]).lift_to_precision(prec)

        poly_to_solve = PolynomialRing(F, name='t1')(coeff_list)
        poly_to_solve_derivative = poly_to_solve.derivative()

        root = root_guess
        for i in range(0, prec):
            root = root - poly_to_solve(root) / poly_to_solve_derivative(root)

        return root

    # def base_ring(self):
    #    return self._given_ground_ring

    def absolute_e(self):
        return self._given_ground_ring.absolute_e()

    def absolute_f(self):
        return self.K1.absolute_f()
 
    def maximal_unramified_subring(self):
        return self.K1

    def defining_polynomial(self, exact=False):
        return self._approx_modulus
    
    def gens(self):
        return [self._given_gen]

class pAdicUnramifiedOverGeneralBaseringInjection(Morphism): 
    def __init__(self, R, S):
        if not R.is_field() or S.is_field():
            Morphism.__init__(self, Hom(R, S))
        else:
            from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
            Morphism.__init__(self, Hom(R, S, SetsWithPartialMaps()))

    def _call_(self, element):
        return self.codomain()([self.codomain()._injection_from_K0_to_K1(coefficient) for coefficient in element.polynomial().list()])

    def _call_with_args(self, x, args=(), kwds={}):
        return self.codomain()(self._call_(x), *args, **kwds)

    def section(self):
        return pAdicUnramifiedOverGeneralBaseringSection(self.codomain(), self.domain())

class pAdicUnramifiedOverGeneralBaseringSection(Morphism):
    def __init__(self, S, R):
        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
        Morphism.__init__(self, Hom(S, R, SetsWithPartialMaps()))

    def _call_(self, element):
        element_in_K1_basis = element.polynomial()
        try:
            element_in_K0_basis = [self.domain()._section_from_K1_to_K0(coefficient) for coefficient in element_in_K1_basis]
        except:
            raise ValueError("Element not contained in base ring.")
        
        return self.domain()._given_ground_ring(element_in_K0_basis)

    def _call_with_args(self, x, args=(), kwds={}):
        return self.domain()(self._call_(x), *args, **kwds)
