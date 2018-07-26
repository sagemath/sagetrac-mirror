from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from .padic_extension_generic import pAdicExtensionGeneric

class RelativeExtensionGeneric(pAdicExtensionGeneric): 
    def __init__(self, poly, prec, print_mode, names, element_class):
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
        self._construct_given_gen(self._unramified_extension_defining_poly.list())
        self._construct_user_representation_of_K1_gen()

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
        
    def _construct_unramified_extension_defining_poly(self):
        coefficients_in_unramified_base = []
        for coefficient in self._approx_modulus.list():
            coefficient_polynomial = coefficient.polynomial()
            if coefficient_polynomial.degree() > 0:
                raise ValueError("Defining polynomial must be defined over the unramified extension.")
            coefficients_in_unramified_base.append(coefficient_polynomial[0])
        self._unramified_extension_defining_poly = PolynomialRing(self._given_ground_ring.ground_ring(), name='t0')(coefficients_in_unramified_base)
    
    def _construct_eisenstein_extension_defining_poly(self):
        self._eisenstein_extension_defining_poly = PolynomialRing(self.K1, name='t1')([coefficient.polynomial()(self._K0_gen) for coefficient in self._given_ground_ring.defining_polynomial().list()])

    def _construct_maximal_unramified_subextension(self):
        K0 = self._given_ground_ring.ground_ring()
        relative_degree = self._approx_modulus.degree()
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

        given_K1_gen = self._given_gen.polynomial()[0] 
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
        element_in_K0_basis = [self._write_in_K0_basis_relative(entry) for entry in element.polynomial().list()]

        return [self._given_ground_ring([element_in_K0_basis[j][i] for j in range(0, len(element_in_K0_basis))]) for i in range(0,  len(element_in_K0_basis[0]) )]

    def _injection_from_K0_to_K1(self, element):
        return element.polynomial()(self._K0_gen)

    def _section_from_K1_to_K0(self, element):
        element_in_K0_basis = self._write_in_K0_basis(element)        

        for j in range(1, len(element_in_K0_basis)):
            if element_in_K0_basis[j] != 0:
                raise ValueError("Element not contained in base field.")

        return element_in_K0_basis[0]

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
         
    def _construct_given_gen(self, unramified_extension_defining_poly_coef_list):
        # elements of unramified_extension_defining_poly_coef_list are in K0, need to send them to K1
        coef_list = [self._injection_from_K0_to_K1(coef) for coef in unramified_extension_defining_poly_coef_list]
        # find root in K1
        root_in_K1 = self._find_root(coef_list, self.K1)
        # lift from K1 to L
        self._given_gen = self([root_in_K1])

    def _extension_type(self):
        return 'Unramified'

#     # if one day, we want to allow unramified extensions by polynomials whose coefficients are not necessarily in K0
#     def _construct_given_gen_v2(self):
#         self._given_gen = self._find_root([_injection_from_given_basering(coefficient) for coefficient in self._approx_modulus], self)
# 
#     # Does not work for now, too confusing to override the base_ring of EisensteinExtensionGeneric.
#     def base_ring(self):
#        return self._given_ground_ring
