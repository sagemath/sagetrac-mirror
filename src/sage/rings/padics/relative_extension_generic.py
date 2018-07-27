"""
Relative Extension Generic

This file implements the shared functionality for three-step extensions of the
form unramified/ramified/unramified.

AUTHORS:

- Vishal Arul
"""
from .padic_extension_generic import pAdicExtensionGeneric
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class RelativeExtensionGeneric(pAdicExtensionGeneric): 
    def __init__(self, poly, prec, print_mode, names, element_class):
        """
        Initializes self.
        The tower structure is named as follows.
          L
        K1  K
          K0
          Qp
        K1/K0/Qp is all unramified
        K/K0 and L/K1 are Eisenstein, with same defining polynomial

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2); L
            5-adic Unramified Extension Field in c defined by s^3 + 2*s^2 + 4*s + 2 over its base field
            sage: c^3 + 2*c^2+ 4*c + 2 == 0
            True
        """
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
        self._construct_given_gen(self._unramified_extension_defining_poly.list())
        self._cache_powers_of_K1_gen_in_given_basis()

    def _pre_init(self, exact_modulus, approx_modulus):
        """
        Computes the data for the extension tower. Needs to be called before
        __init__. Should be called exactly once. For internal use.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: L._unramified_extension_defining_poly # indirect doctest
            (1 + O(5^20))*t0^3 + (2 + O(5^20))*t0^2 + (4 + O(5^20))*t0 + 2 + O(5^20)
            sage: L.maximal_unramified_subextension()
            5-adic Unramified Extension Field in a1 defined by x^6 + x^4 + 4*x^3 + x^2 + 2
            sage: L._eisenstein_extension_defining_poly # indirect doctest
            (1 + O(5^20))*t1^2 - 5 + O(5^21)
        """
        self._exact_modulus = exact_modulus # unnecessary?
        self._approx_modulus = approx_modulus
        self._given_ground_ring = approx_modulus.base_ring()         
        self._construct_unramified_extension_defining_poly()
        self._construct_maximal_unramified_subextension()
        self._construct_eisenstein_extension_defining_poly()
        
    def absolute_e(self):
        """
        Return the absolute ramification index of this ring or field.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: K0.absolute_e()
            1
            sage: K.absolute_e()
            2
            sage: L.absolute_e()
            2
        """
        return self._given_ground_ring.absolute_e()

    def absolute_f(self):
        """
        Return the degree of the residue field of this ring/field
        over its prime subfield.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: K0.absolute_f()
            2
            sage: K.absolute_f()
            2
            sage: L.absolute_f()
            6
        """
        return self.K1.absolute_f()
 
    def maximal_unramified_subextension(self):
        """
        Return the maximal unramified subextension.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: L.maximal_unramified_subextension()
            5-adic Unramified Extension Field in a1 defined by x^6 + x^4 + 4*x^3 + x^2 + 2
        """
        return self.K1

    def defining_polynomial(self, exact=False):
        """
        Return the user-defined polynomial for this extension. 

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: L.defining_polynomial()
            (1 + O(b^40))*s^3 + (2 + O(b^40))*s^2 + (4 + O(b^40))*s + 2 + O(b^40)
        """
        if exact:
            return self._exact_modulus
        return self._approx_modulus
    
    def gens(self, n=0):
        """
        Return a list of generators for self as an extension of the
        user-provided ground ring. There should only be a single generator.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: L.gens()
            [c + O(b^40)]
        """
        if n != 0:
            raise IndexError("only one generator")
        return [self._given_gen]
        
    def _extension_type(self):
        """
        Return the type (``Unramified``, ``Eisenstein``) of this 
        extension as a string, if any.

        Used for printing.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: K0._extension_type()
            'Unramified'
            sage: K._extension_type()
            'Eisenstein'
            sage: L.K1._extension_type()
            'Unramified'
            sage: L._extension_type()
            'Unramified'
        """
        return 'Unramified'

    def _construct_unramified_extension_defining_poly(self):
        """
        Construct the polynomial that defines the extension K1/Qp. For internal
        use. Should be called only once (in the constructor).

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: L._unramified_extension_defining_poly # indirect doctest
            (1 + O(5^20))*t0^3 + (2 + O(5^20))*t0^2 + (4 + O(5^20))*t0 + 2 + O(5^20)
        """
        coefficients_in_unramified_base = []
        for coefficient in self._approx_modulus.list():
            coefficient_polynomial = coefficient.polynomial()
            if coefficient_polynomial.degree() > 0:
                raise ValueError("Defining polynomial must be defined over the unramified extension.")
            coefficients_in_unramified_base.append(coefficient_polynomial[0])
        self._unramified_extension_defining_poly = PolynomialRing(self._given_ground_ring.ground_ring(), name='t0')(coefficients_in_unramified_base)
    
    def _construct_maximal_unramified_subextension(self):
        """
        Construct the extension L/K1 and expresses a generator of K0/Qp in
        terms of a generator of K1/Qp. For internal use. Should be called only
        once (in the constructor).

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: L.K1 # indirect doctest
            5-adic Unramified Extension Field in a1 defined by x^6 + x^4 + 4*x^3 + x^2 + 2
            sage: L._given_ground_ring.base_ring().defining_polynomial()(L._K0_gen) # indirect doctest
            0 + O(5^20)
        """
        K0 = self._given_ground_ring.ground_ring()
        relative_degree = self._approx_modulus.degree()
        self.K1 = K0.change(q = K0.residue_characteristic()**(relative_degree * K0.degree()), names=K0.variable_names()[0] + '1')
        self._K0_gen = self._find_root([self.K1(coefficient) for coefficient in K0.defining_polynomial().list()])
        K0_gen_coefficients = self._K0_gen.polynomial()

    def _construct_eisenstein_extension_defining_poly(self):
        """
        Construct the polynomial that defines the extension L/K1. For internal
        use. Should be called only once (in the constructor).

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: L._eisenstein_extension_defining_poly # indirect doctest
            (1 + O(5^20))*t1^2 - 5 + O(5^21)
        """
        self._eisenstein_extension_defining_poly = PolynomialRing(self.K1, name='t1')([coefficient.polynomial()(self._K0_gen) for coefficient in self._given_ground_ring.defining_polynomial().list()])

    def _construct_given_gen(self, unramified_extension_defining_poly_coeff_list):
        """
        Construct a generator of K1/K0 (that is, a root of the polynomial
        provided by the user) in K1. For internal use. Should be called only
        once (in the constructor).

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: PolynomialRing(L,name='t')([L(coeff) for coeff in L._approx_modulus])(L._given_gen) # indirect doctest
            0 + O(b^40)
        """
        coeff_list = [self._injection_from_K0_to_K1(coeff) for coeff in unramified_extension_defining_poly_coeff_list]
        root_in_K1 = self._find_root(coeff_list)
        self._given_gen = self([root_in_K1])

    def _cache_powers_of_K1_gen_in_given_basis(self):
        """
        Cache powers of a generator of K1/Qp as K0-linear combinations of a
        generator of K1/K0. For internal use. Should be called only once (in
        the constructor).

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: R1 = PolynomialRing(L.K1,names='t1')
            sage: given_K1_gen = L._given_gen.polynomial()(0)
            sage: [R1([L._injection_from_K0_to_K1(coeff) for coeff in polynom])(given_K1_gen) for polynom in L._powers_of_K1_gen_in_given_basis] # indirect doctest
            [1 + O(5^20),
             a1 + O(5^20),
             a1^2 + O(5^20),
             a1^3 + O(5^20),
             a1^4 + O(5^20),
             a1^5 + O(5^20)]
            sage: _[1] == L.K1.gen() # indirect doctest
            True
        """
        K0 = self._given_ground_ring.ground_ring()
        given_K1_gen = self._given_gen.polynomial()[0] 
        relative_degree = self._approx_modulus.degree()
        K0_degree = self._given_ground_ring.ground_ring().degree()
        target_basis = [self.K1(1)] + [0] * (K0_degree * relative_degree - 1)
        for i in range(1, K0_degree):
            target_basis[i] = target_basis[i - 1] * self._K0_gen
        for j in range(1, relative_degree):
            for i in range(0, K0_degree):
                target_basis[j * K0_degree + i] = target_basis[(j - 1) * K0_degree + i] * given_K1_gen

        # now target_basis[j * K0_degree + i] should equal
        # (self._K0_gen)^i * (given_K1_gen)^j 

        target_basis_vectors = []
        for basis_entry in target_basis:
            polynomial_list = basis_entry.polynomial().list()
            polynomial_list += [0] * (self.K1.degree() - len(polynomial_list))
            target_basis_vectors += polynomial_list
        
        from sage.matrix.constructor import matrix
        target_basis_matrix = matrix(self._given_ground_ring.ground_ring().ground_ring(), self.K1.degree(), self.K1.degree(), target_basis_vectors)
        target_basis_matrix_inverse = target_basis_matrix.inverse()

        R0 = PolynomialRing(K0, name='t0')        
        self._powers_of_K1_gen_in_given_basis = []
        for i in range(0, self.K1.degree()):
            target_basis_matrix_inverse_row = target_basis_matrix_inverse.row(i)
            self._powers_of_K1_gen_in_given_basis.append(R0([K0(list(target_basis_matrix_inverse_row[j*K0_degree:(j+1)*K0_degree])) for j in range(0, relative_degree)]))

    def _express_via_given_gen(self, element):
        """
        Construct the representation of element as a K-linear combination of
        powers of the user-given generator of L/K.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: x = a^2*c + a^7*b^6*c^5 - c^12 + 7
            sage: S(L._express_via_given_gen(x))(c) == x
            True
        """
        element_in_K0_basis = [self._write_in_K0_basis(entry) for entry in element.polynomial().list()]
        relative_degree = self._approx_modulus.degree()
        for j in range(0, len(element_in_K0_basis)):
            short_length = len(element_in_K0_basis[j])
            element_in_K0_basis[j] += [0] * (relative_degree - short_length)

        maximal_precision = self._given_gen.precision_absolute() / self.absolute_e()
        return [self._given_ground_ring([element_in_K0_basis[j][i] for j in range(0, len(element_in_K0_basis))]).add_bigoh(maximal_precision) for i in range(0, relative_degree)]

    def _write_in_K0_basis(self, element):
        """
        Construct the representation of an element of K1 as a K0-linear
        combination of powers of a generator of K1/K0 (this generator is a root
        of the user-given polynomial.)

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: a1 = L.K1.gen()
            sage: x = a1^7 + L._injection_from_K0_to_K1(a^3-2)*a1^2 + 1
            sage: R1 = PolynomialRing(L.K1,names='t1')
            sage: given_K1_gen = L._given_gen.polynomial()(0)
            sage: R1([L._injection_from_K0_to_K1(coeff) for coeff in L._write_in_K0_basis(x) ])(given_K1_gen) == x
            True
        """
        if (element.parent() != self.K1):
            raise ValueError("Only elements of K1 can be expressed in the K0-basis.")
        element_in_Qp_basis = element.polynomial().list()
        if (len(element_in_Qp_basis) == 0):
            return [self._given_ground_ring.ground_ring()(0)]
        return sum( element_in_Qp_basis[i] * self._powers_of_K1_gen_in_given_basis[i] for i in range(0, len(element_in_Qp_basis))).list()

    def _injection_from_K0_to_K1(self, element):
        """
        Takes an element of K0 and maps it to K1.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: a+2/a
            -4 + O(5^20)
            sage: L._injection_from_K0_to_K1(a+2/a)
            -4 + O(5^20)
        """
        return element.polynomial()(self._K0_gen)

    def _section_from_K1_to_K0(self, element):
        """
        Takes an element of K1 and maps it to K0 if possible.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: L._section_from_K1_to_K0(L._injection_from_K0_to_K1(a+2/a)) == a+2/a
            True
        """
        element_in_K0_basis = self._write_in_K0_basis(element)        

        if (len(element_in_K0_basis) == 0):
            return self._given_ground_ring.base_ring()(0)

        for j in range(1, len(element_in_K0_basis)):
            if element_in_K0_basis[j] != 0:
                raise ValueError("Element not contained in base field.")

        return element_in_K0_basis[0]

    def _find_root(self, coeff_list):
        """
        Given a polynomial over K1 whose coefficients are in coeff_list,
        attempts to find a simple root of this polynomial in K1 (assuming that
        it has a simple root in K1) by finding a root in the residue field and
        then taking a Hensel lift of that. For internal use.

        EXAMPLES::

            sage: K0.<a> = Qq(25, print_pos=False,print_mode='terse')
            sage: R0.<t> = PolynomialRing(K0)
            sage: K.<b> = K0.extension(t^2 - 5)
            sage: R.<s> = PolynomialRing(K)
            sage: L.<c> = K.extension(s^3 + 2*s^2 + 4*s + 2)
            sage: x = L._find_root([L.K1(-2), L.K1(0), L.K1(1)])
            sage: x^2-2 == 0
            True
        """
        prec = self.K1.precision_cap()
        coeff_residue_list = [ coefficient.residue() for coefficient in coeff_list ]
        root_guess = self.K1(PolynomialRing(self.K1.residue_field(), name='t1b')(coeff_residue_list).roots()[0][0]).lift_to_precision(prec)

        poly_to_solve = PolynomialRing(self.K1, name='t1')(coeff_list)
        poly_to_solve_derivative = poly_to_solve.derivative()

        root = root_guess
        for i in range(0, prec):
            root = root - poly_to_solve(root) / poly_to_solve_derivative(root)

        return self.K1(root)
         
#       # Does not work for now, too confusing to override the base_ring of
#       # EisensteinExtensionGeneric.
#       def base_ring(self):
#          return self._given_ground_ring
#   
#       # Not sure if this is needed. Maybe we need this for some coercion issue?
#       def _pushout_(self, F):
#           if (self._given_ground_ring is F):
#               return self
