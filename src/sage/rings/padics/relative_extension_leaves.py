"""
Relative extensions of `p`-adic rings
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

class pAdicRelativeBaseringInjection(Morphism):
    def __init__(self, R, S):
        if not R.is_field() or S.is_field():
            Morphism.__init__(self, Hom(R, S, R.category()))
        else:
            from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
            Morphism.__init__(self, Hom(R, S, SetsWithPartialMaps()))

    def _call_(self, x):
        if x.is_zero():
            return self.codomain()(0,x.precision_absolute())
        else:
            return self.codomain()([x])

    def _call_with_args(self, x, args=(), kwds={}):
        return self.codomain()([x], *args, **kwds)

    def section(self):
        return pAdicRelativeBaseringSection(self.codomain(), self.domain())

class pAdicRelativeBaseringSection(Morphism):
    def __init__(self, S, R):
        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
        Morphism.__init__(self, Hom(S, R, SetsWithPartialMaps()))

    def _call_(self, x):
        f = x.polynomial()
        if f.degree() > 0:
            raise ValueError("Element not contained in base ring")
        return f[0]

    def _call_with_args(self, x, args=(), kwds={}):
        return self.codomain()(self._call_(x), *args, **kwds)

class RelativeRamifiedExtensionRingFixedMod(EisensteinExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(field=False, show_prec=False, type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'fixed-mod')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedFixedModElement)
        from .relative_ramified_FM import pAdicCoercion_ZZ_FM, pAdicConvert_QQ_FM
        self.register_coercion(pAdicCoercion_ZZ_FM(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_FM(self))

class RelativeRamifiedExtensionRingCappedAbsolute(EisensteinExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(field=False, show_prec=False, type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-abs')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedCappedAbsoluteElement)
        from .relative_ramified_CA import pAdicCoercion_ZZ_CA, pAdicConvert_QQ_CA
        self.register_coercion(pAdicCoercion_ZZ_CA(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_CA(self))

class RelativeRamifiedExtensionRingCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeRingGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
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
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(field=False, show_prec=False, type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, True, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-rel')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedCappedRelativeElement)
        from .relative_ramified_CR import pAdicCoercion_ZZ_CR, pAdicCoercion_QQ_CR
        self.register_coercion(pAdicCoercion_ZZ_CR(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_coercion(pAdicCoercion_QQ_CR(self))

class RelativeRamifiedExtensionRingFloatingPoint(EisensteinExtensionGeneric, pAdicFloatingPointRingGeneric):
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
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
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(field=False, show_prec=False)
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, True, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-rel')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedFloatingPointElement)
        from .relative_ramified_FP import pAdicCoercion_ZZ_FP, pAdicCoercion_QQ_FP
        self.register_coercion(pAdicCoercion_ZZ_FP(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
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

        eisenstein_extension_defining_poly = PolynomialRing(self.K1, name='t1')([coefficient.polynomial()(self._K0_gen) for coefficient in K.defining_polynomial().coefficients()])

        EisensteinExtensionGeneric.__init__(self, eisenstein_extension_defining_poly, prec, print_mode, names, RelativeRamifiedCappedRelativeElement)

        # self._construct_given_gen_v2()
        self._construct_given_gen(unramified_extension_defining_poly_coef_list)

        self._construct_user_representation_of_K1_gen()
    
    def _construct_maximal_unramified_subextension(self):
        K0 = self._given_ground_ring.ground_ring()
        relative_degree = self._approx_modulus.degree()

        # construct K1 (the maximal unramified extension)
        self.K1 = Qq(K0.residue_characteristic()**(relative_degree * K0.degree()), prec = K0.precision_cap(), type = 'capped-rel', names='a1')

        # find a way to map a generator of K0 into K1
        self._K0_gen = _find_root([self.K1(coefficient) for coefficient in K0.defining_polynomial().list()], self.K1)

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

        # this is alpha
        given_K1_gen = self.K1(self._given_gen[0])
        relative_degree = self._approx_modulus.degree()
        # this is beta
        self._K0_gen
        K0_degree = self._given_ground_ring.ground_ring().degree()

        target_basis = [0] * (K0_degree * relative_degree)
        target_basis[0] = self.K1(1)

        for i in range(1, K0_degree):
            target_basis[i] = target_basis[i - 1] * self._K0_gen
        for j in range(1, relative_degree):
            for i in range(0, K0_degree):
                target_basis[j * K0_degree + i] = target_basis[(j - 1) * K0_degree + i] * given_K1_gen

        # now target_basis[j * K0_degree + i] should equal
        # (self._K0_gen)^i * (given_K1_gen)^j 

        target_basis_matrix = matrix(self._given_ground_ring.ground_ring().ground_ring(), self.K1.degree(), self.K1.degree(), [basis_entry.polynomial().list() for basis_entry in target_basis])
        target_basis_matrix_inverse = target_basis_matrix.inverse()

        R0 = PolynomialRing(K0, name='t0')

        self._powers_of_K1_gen_in_given_basis = [R0([K0(target_basis_matrix_inverse.column(i)[j*K0_degree:(j+1)*K0_degree]) for j in range(0, relative_degree)]) for i in range(0, self.K1.degree())]

    def _injection_from_K0_to_K1(self, element):
        return element.polynomial()(self._K0_gen)

    def _write_in_K0_basis_relative(self, element): # element is in K1
        element_in_Qp_basis = element.polynomial().list()
        return sum( element_in_Qp_basis[i] * self._powers_of_K1_gen_in_given_basis[i] for i in range(0, self.K1.degree())).list()

    def _write_in_K0_basis(self, element): # element is in K1
        if (element.parent() != self.K1):
            raise ValueError("Only elements of K1 can be expressed in the K0-basis.")
        K0 = self._given_ground_ring.ground_ring() 
        relative_degree = self._approx_modulus.degree()

        element_in_Qp_basis = element.polynomial()

        return [ sum( element_in_Qp_basis[i] * self._powers_of_K1_gen[i][j] for i in range(0, K1.degree())) for j in range(0, relative_degree) ]

    def _section_from_K1_to_K0(self, element):
        element_in_K0_basis = self._write_in_K0_basis(element)        

        for j in range(1, len(element_in_K0_basis)):
            if element_in_K0_basis[j] != 0:
                raise ValueError("Element not contained in base field.")

        return element_in_K0_basis[0]

    def _injection_from_given_basering(self, element):
        return L([self._injection_from_K0_to_K1(coefficient) for coefficient in element.polynomial().list()])

    def _section_to_given_basering(self, element):
        element_in_K1_basis = element.polynomial()
        try:
            element_in_K0_basis = [self._section_from_K1_to_K0(coefficient) for coefficient in element_K1_basis]
        except:
            raise ValueError("Element not contained in base ring.")
        
        return self._given_ground_ring(element_in_K0_basis)

    def _express_via_given_gen(self, element):
        element_in_K1_basis = element.polynomial().list()
        eisenstein_degree = len(element_in_K1_basis)
        element_in_K0_basis = [self._write_in_K0_basis_relative(entry) for entry in element_in_K1_basis]
        relative_degree = self._approx_modulus.degree()
        K = self._given_ground_ring

        return [K([element_in_K0_basis[j][i] for j in range(0, eisenstein_degree)]) for i in range(0, relative_degree)]
        
         
#     # if one day, we want to allow unramified extensions by polynomials whose coefficients are not necessarily in K0
#     def _construct_given_gen_v2(self):
#         self._given_gen = _find_root([_injection_from_given_basering(coefficient) for coefficient in self._approx_modulus], self)
# 
    def _construct_given_gen(self, unramified_extension_defining_poly_coef_list, K1):
        # find root in K1
        root_in_K1 = _find_root(unramified_extension_defining_poly_coef_list, self.K1)
        # lift from K1 to L
        self._given_gen = self([root_in_K1])
    
    def _find_root(coeff_list, F):
        prec = F.precision_cap()
        coeff_residue_list = [ coefficient.residue() for coefficient in coeff_list ]
        root_guess = F(PolynomialRing(F.residue_field(), name='t1b')(coeff_residue_list).roots()[0][0]).lift_to_precision(prec)

        poly_to_solve = PolynomialRing(F, name='t1')(coeff_list)
        poly_to_solve_derivative = poly_to_solve.derivative()

        root = root_guess
        for i in range(0, prec):
            root = root - poly_to_solve(root) / poly_to_solve_derivative(root)

        return root

    # stuff to do:
    #  printing
    #  basering injection, section

    def ground_ring(self):
        return self._given_ground_ring

    def defining_polynomial(self, exact=False):
        return self._approx_modulus
    
    def gens(self):
        return [self._given_gen]
