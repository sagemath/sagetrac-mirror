r"""
Multivariate Polynomials algebra as Module over Symmetric Functions
"""
#*****************************************************************************
#  Copyright (C) 2012 Nicolas Borie <nicolas.borie at math.u-psud dot fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_function
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.bindable_class import BindableClass
from string import join
import sage.libs.symmetrica.all as symmetrica

from sage.categories.algebras import Algebras
from sage.categories.finite_dimensional_algebras_with_basis import FiniteDimensionalAlgebrasWithBasis
from sage.categories.modules import Modules
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.realizations import Realizations

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutation, from_lehmer_code
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import Matrix

@cached_function
def harmonic_polynomial_for_permutation_group_element(R, perm):
    r"""
    Returns the harmonic polynomial indexed by ``perm`` in the ring
    ``R``.

    ..math:: P_{perm} = \rho~\mathbf{x}^{\epsilon} \Delta_n(\mathbf{x})

    where \epsilon is the lehmer_code of the permutation `w_0 * perm`
    and `w_0` is the longest element of the symmetric group of degree
    `n`.

    EXAMPLES::

        sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
        sage: R = A.polynomial_algebra()
        sage: from sage.combinat.multivariate_polynomials.abstract_sym_module import harmonic_polynomial_for_permutation_group_element
        sage: for p in SymmetricGroup(3): harmonic_polynomial_for_permutation_group_element(R, p)
        2
        2*x1 - 2*x2
        2*x2 - 2*x3
        x1^2 - 2*x1*x2 + 2*x2*x3 - x3^2
        2*x1*x2 - x2^2 - 2*x1*x3 + x3^2
        x1^2*x2 - x1*x2^2 - x1^2*x3 + x2^2*x3 + x1*x3^2 - x2*x3^2
    """
    S = perm.parent()
    exp = Permutation(S.long_element()*perm).to_lehmer_code()
    n = len(R.gens())
    # Let us define Delta as the Vandermonde in n variables
    Delta = R.prod([R.gens()[i] - R.gens()[j] for i in range(n) for j in range(n) if i < j])
    # We dummy apply the rigth number of partial derivatives
    for i in range(n):
        for j in range(exp[i]):
            Delta = Delta.derivative(R.gens()[i])
    return Delta

def divided_difference_polynomial_ring(P, i):
    r"""
    Return the `i^{th}` divided difference of ``P`` which must be a
    polynomial ring element.

    EXAMPLES::

        sage: from sage.combinat.multivariate_polynomials.abstract_sym_module import divided_difference_polynomial_ring as dd
        sage: K.<x1,x2,x3> = PolynomialRing(QQ, 'x1,x2,x3')
        sage: P = x1^2*x2; dd(P,1)                                
        x1*x2
        sage: dd(P,2)
        x1^2
        sage: dd(dd(P,1),2)
        x1
        sage: dd(dd(P,2),1)
        x1 + x2
        sage: dd(dd(dd(P,1),2),1)
        1
        sage: dd(dd(dd(P,2),1),2)
        1
    """
    R = P.parent()
    t = R.gens()[:i-1]+(R.gens()[i],)+(R.gens()[i-1],)+R.gens()[i+1:]
    ans = P - P(t)
    return R(ans/(R.gens()[i-1] - R.gens()[i]))

def maximal_symmetrizer(P):
    r"""
    Returns the image of the polynomial ``P`` in `n` variables by the
    maximal symmetrizer operator. The maximal symmetrizer operateur is
    the composition of divided difference indexed by the letter of a
    reduced word of the longest element in the symmetric group of
    degree `n`.

    EXAMPLES::

        sage: from sage.combinat.multivariate_polynomials.abstract_sym_module import maximal_symmetrizer                     
        sage: K.<x1,x2,x3> = PolynomialRing(QQ, 'x1,x2,x3')
        sage: maximal_symmetrizer(x1^2*x2)
        1
        sage: maximal_symmetrizer( (x1^2*x2) * (1) )
        1
        sage: maximal_symmetrizer( (x1^2*x2) * (x1+x2+x3) )
        x1 + x2 + x3
        sage: maximal_symmetrizer( (x1^2*x2) * (x1^2+x2+x3) )
        x1^2 + x1*x2 + x2^2 + x1*x3 + x2*x3 + x3^2
        sage: maximal_symmetrizer( x1*x2 + x1^2 + x2^2 + x1 + x2 + x3 )
        0
    """
    n = len(P.parent().gens())
    for i in range(n-1,0,-1):
        for j in range(1,i+1):
            P = divided_difference_polynomial_ring(P, j)
    return P;

def coeficient_schur_schubert_from_polynomial(P, perm):
    r"""
    Returns the coeficient of ``P`` among the Schubert polynomial
    indexed by ``perm`` in the Schur / Schubert decomposition of the
    polynomial ``P``.

    EXAMPLES::

        sage: from sage.combinat.multivariate_polynomials.abstract_sym_module import coeficient_schur_schubert_from_polynomial
        sage: K.<x1,x2,x3> = PolynomialRing(QQ, 'x1,x2,x3')
        sage: for p in SymmetricGroup(3): print coeficient_schur_schubert_from_polynomial(x1^2+x2+4*x3*x2, p)
        4*s[1, 1]
        s[]
        -s[] - 4*s[1]
        0
        5*s[]
        0
    """
    S = perm.parent()
    R = P.parent()
    Q = (perm.sign())*R(symmetrica.t_SCHUBERT_POLYNOM(list(Permutation(S.long_element()*perm))))
    t = R.gens()[::-1]
    SF = SymmetricFunctions(R.base_ring()).schur()
    return SF.from_polynomial(maximal_symmetrizer(P(t)*Q))

@cached_function
def schubert_polynomial_to_monomial_under_staircase(K, perm):
    r"""
    Return a dictionnary

    EXAMPLES::

        sage: from sage.combinat.multivariate_polynomials.abstract_sym_module import schubert_polynomial_to_monomial_under_staircase
        sage: S = SymmetricGroup(3)
        sage: for s in S: print schubert_polynomial_to_monomial_under_staircase(QQ, s)
        {(): s[]}
        {(1,2): s[], (2,3): s[]}
        {(1,2): s[]}
        {(1,2,3): s[]}
        {(1,3,2): s[]}
        {(1,3): s[]}
    """
    p = Permutation(perm)
    R = PolynomialRing(ZZ, perm.parent().degree(), 'x')
    pol = R(symmetrica.t_SCHUBERT_POLYNOM(p))
    d = dict()
    Sn = SymmetricGroup(perm.parent().degree())
    SF = SymmetricFunctions(K).schur()
    for coef, mon in pol:
        d[Sn(from_lehmer_code(list(mon.exponents()[0])))] = SF(coef);
    return d


class AbstractPolynomialRingAsSymModule(UniqueRepresentation, Parent):
    r"""
    The ring of multivariate polynomials in `n` formal variables view
    as a free module of rank `n!` over the ring of symmetric
    polynomials.

    EXAMPLES::

        sage: A = AbstractPolynomialRingAsSymModule(QQ, 4); A
        Polynomial ring in x1, x2, x3, x4 over Rational Field view as module over symmetric polynomials over Rational Field
        sage: A.category()
        Join of Category of additive unital additive magmas with realizations and Category of modules over Symmetric Functions over Rational Field
        sage: A.polynomial_algebra()
        Multivariate Polynomial Ring in x1, x2, x3, x4 over Rational Field

    TESTS::

        sage: A = AbstractPolynomialRingAsSymModule(QQ, 5)
        sage: TestSuite(A).run()
    """
    def __init__(self, R, nb_variables, main_repr_var = 'x'):
        r"""
        TESTS::

            sage: A = AbstractPolynomialRingAsSymModule(QQ, 5)
        """
        Parent.__init__(self, base=R, category = (Algebras(SymmetricFunctions(R)).WithRealizations(),))
        self._nb_variables = nb_variables
        self._main_repr_var = main_repr_var

        S = self.Schubert()
        M = self.Staircase()
        H = self.Harmonic()
        D = self.Descent()

    def _repr_(self):
        r"""
        TESTS::

            sage: AbstractPolynomialRingAsSymModule(QQ, 4)
            Polynomial ring in x1, x2, x3, x4 over Rational Field view as module over symmetric polynomials over Rational Field
        """
        vars_string = join([self._main_repr_var+str(i) for i in range(1, self.nb_variables()+1)], ", ")
        return "Polynomial ring in %s over %s view as module over symmetric polynomials over %s"%(vars_string, self.base_ring(), self.base_ring())

    def a_realization(self):
        r"""
        Returns a realization of ``self``.

        EXAMPLES::

            sage: A = AbstractPolynomialRingAsSymModule(QQ, 3)
            sage: A.a_realization()
            Multivariate Polynomial in x1, x2, x3 over Symmetric Functions over Rational Field in the Schur basis in Schur Schubert decomposition
        """
        return self.Schubert()

    @cached_method
    def nb_variables(self):
        r"""
        Returns the number of variables of ``self``.

        EXAMPLES::

            sage: A = AbstractPolynomialRingAsSymModule(QQ, 5);
            sage: A.nb_variables()
            5
        """
        return self._nb_variables

    @cached_method
    def polynomial_algebra(self):
        r"""
        Returns the usual polynomial algebra of multivariate polynomials.

        EXAMPLES::

            sage: A = AbstractPolynomialRingAsSymModule(QQ, 5);
            sage: A.polynomial_algebra()
            Multivariate Polynomial Ring in x1, x2, x3, x4, x5 over Rational Field
        """
        return PolynomialRing(self.base_ring(), [self._main_repr_var+str(i) for i in range(1, 1+self.nb_variables())])

    @cached_method
    def quotient_by_symmetric_polynomials(self):
        r"""
        Returns the abstract quotient of multivariate polynomials by
        the ideal of constant free symmetric polynomials.

        EXAMPLES::

        """
        from sage.combinat.multivariate_polynomials.coinvariants_symmetric_group import Abstract_Quotient_MulPol_as_sym_module
        return Abstract_Quotient_MulPol_as_sym_module(self)

    class Bases(Category_realization_of_parent):
        r"""
        The category of bases as free module over symmetric
        polynomials of the multivariate polynomials. 

        TESTS::


        """
        def super_categories(self):
            r"""
            """
            R = self.base().base_ring()
            return [FiniteDimensionalAlgebrasWithBasis(SymmetricFunctions(R)), Realizations(self.base())]

        class ParentMethods:
            @cached_method
            def one_basis(self):
                r"""
                returns the key indexing the one for the
                multiplication in each realization.

                EXAMPLES::

                    sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                    sage: B = A.Schubert()
                    sage: B.one()                             
                    s[]*Y_()
                    sage: B.one().to_multivariate_polynomial()
                    1
                    sage: M = A.Staircase()
                    sage: M.one()
                    s[]*(1)
                    sage: M.one().to_multivariate_polynomial()
                    1
                """
                return self.basis().keys().one()

            @cached_method
            def nb_variables(self):
                r"""
                Returns the number of variables of ``self``.

                EXAMPLES::

                    sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                    sage: M = A.staircase()
                    sage: M.nb_variables()
                    3
                """
                return self.category().base().nb_variables()

            @cached_method
            def polynomial_algebra(self):
                r"""
                Returns the usual polynomial algebra of multivariate polynomials.
                
                EXAMPLES::

                    sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                    sage: M = A.staircase()
                    sage: M.polynomial_algebra()
                    Multivariate Polynomial Ring in x1, x2, x3 over Rational Field                
                """
                return self.category().base().polynomial_algebra()

            @cached_method
            def transition_matrix_from(self, realization):
                r"""


                EXAMPLES::


                """
                assert(realization.realization_of() is self.realization_of()), "%s and %s must be two realizations of the same abstract module"
                if self.has_coerce_map_from(realization):
                    SF = self.base_ring()
                    L = []
                    keys = self.basis().keys()
                    for k1 in keys:
                        V = []
                        img = self(realization.basis()[k1])
                        for k2 in keys:
                            V.append(img.coefficient(k2))
                        L.append(V)
                    return Matrix(SF, L)
                else:
                    if realization.has_coerce_map_from(self):
                        M = realization.transition_matrix_from(self)
                        return (~M.det())*M.adjoint()
                    else:
                        raise ValueError, "No transition matrix can be computed since there is no coercion"

        class ElementMethods:
            def to_multivariate_polynomial(self):
                r"""
                Returns ``self`` expressed in the polynomial algebra
                modelized by the parent of ``self``.

                EXAMPLES::

                    sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                    sage: B = A.schur_schubert()
                    sage: p = B.sum([B.basis()[k] for k in B.basis().keys()])
                    sage: p.to_multivariate_polynomial()
                    x1^2*x2 + x1^2 + x1*x2 + 2*x1 + x2 + 1
                """
                R = self.parent().category().base().polynomial_algebra()
                return R.sum([R(coef.expand(self.parent().nb_variables()))*self.parent().basis_element_to_multivariate_polynomial(key) for key,coef in self])

            def to_quotient_by_symmetric_polynomials(self):
                r"""

                EXAMPLES::

                """
                d = dict()
                for key, coef in self:
                    d[key] = coef.coefficient([])
                return self.parent().quotient_by_symmetric_polynomials()._from_dict(d)

    class Schubert(CombinatorialFreeModule, BindableClass):
        r"""
        The Schur / Schubert basis of the multivariate polynomials in
        a finite number of variables.

        EXAMPLES::


        TESTS::

        """
        def __init__(self, AbsMod):
            r"""
            TESTS::


            """
            CombinatorialFreeModule.__init__(self, SymmetricFunctions(AbsMod.base_ring()).schur(),
               SymmetricGroup(AbsMod.nb_variables()), prefix='Y_', latex_prefix='Y',
               latex_bracket=False, bracket=False, category=AbsMod.Bases())
            self._main_repr_var = AbsMod._main_repr_var

            M = self.realization_of().staircase()
            # Coercion from Schubert to Staircase
            R = self.base_ring().base_ring()
            f = self.module_morphism(on_basis=lambda x : M._from_dict(schubert_polynomial_to_monomial_under_staircase(R, x)), 
                                     codomain=self.realization_of().staircase())
            f.register_as_coercion()

            # Coercion from Staircase to Schubert
            f = M.module_morphism(on_basis=lambda x : self.from_polynomial(M.basis_element_to_multivariate_polynomial(x)), 
                                     codomain=self)
            f.register_as_coercion()

            D = self.realization_of().Descent()
            # Coercion from Descent to Schubert
            f = D.module_morphism(on_basis=lambda x : self.from_polynomial(D.basis_element_to_multivariate_polynomial(x)), 
                                     codomain=self)
            f.register_as_coercion()

            H = self.realization_of().Harmonic()
            # Coercion from Harmonic to Schubert
            f = H.module_morphism(on_basis=lambda x : self.from_polynomial(H.basis_element_to_multivariate_polynomial(x)), 
                                     codomain=self)
            f.register_as_coercion()

        def _repr_(self):
            r"""
            TESTS::

                sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                sage: A.schur_schubert()
                Multivariate Polynomial in x1, x2, x3 over Symmetric Functions over Rational Field in the Schur basis in Schur Schubert decomposition
            """
            vars_string = join([self._main_repr_var+str(i) for i in range(1, self.nb_variables()+1)], ", ")
            return "Multivariate Polynomial in %s over %s in Schur Schubert decomposition"%(vars_string, self.base_ring())

        def basis_element_to_multivariate_polynomial(self, key):
            r"""
            Returns the element basis of ``self`` indexed by the key
            ``key`` expressed in the polynomial algebra modelized by
            ``self``.

            EXAMPLES::
                sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                sage: B = A.schur_schubert()
                sage: for k in B.basis().keys():
                ....:     B.basis_element_to_multivariate_polynomial(k)
                1
                x1 + x2
                x1
                x1*x2
                x1^2
                x1^2*x2
            """
            return self.polynomial_algebra()(symmetrica.t_SCHUBERT_POLYNOM(list(Permutation(key))))

        def quotient_by_symmetric_polynomials(self):
            r"""


            EXAMPLES::


            """
            return self.realization_of().quotient_by_symmetric_polynomials().Schubert()

        def from_polynomial(self, P):
            r"""

            EXAMPLES::


            """
            assert(len(P.parent().gens()) == self.nb_variables())
            d = dict()
            for k in self.basis().keys():
                d[k] = coeficient_schur_schubert_from_polynomial(P, k)
            return self._from_dict(d)

        def product_on_basis(self, k1, k2):
            r"""

            EXAMPLES::


            """
            P1 = self.basis_element_to_multivariate_polynomial(k1)
            P2 = self.basis_element_to_multivariate_polynomial(k2)
            return self.from_polynomial(P1*P2)

    S = schur_schubert = Schubert

    class Harmonic(CombinatorialFreeModule, BindableClass):
        r"""
        

        EXAMPLES::


        TESTS::


        """
        def __init__(self, AbsMod):
            r"""
            TESTS::

                sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                sage: A.Harmonic()
                Multivariate Polynomial in x1, x2, x3 over Symmetric Functions over Rational Field in the Schur basis expressed with harmonic polynomials
            """
            CombinatorialFreeModule.__init__(self, SymmetricFunctions(AbsMod.base_ring()).schur(),
               SymmetricGroup(AbsMod.nb_variables()), prefix='Y_', latex_prefix='Y',
               latex_bracket=False, bracket=False, category=AbsMod.Bases())
            self._main_repr_var = AbsMod._main_repr_var

        def _repr_(self):
            r"""
            TESTS::


            """
            vars_string = join([self._main_repr_var+str(i) for i in range(1, self.nb_variables()+1)], ", ")
            return "Multivariate Polynomial in %s over %s expressed with harmonic polynomials"%(vars_string, self.base_ring())

        def basis_element_to_multivariate_polynomial(self, key):
            r"""
            Returns the element basis of ``self`` indexed by the key
            ``key`` expressed in the polynomial algebra modelized by
            ``self``.

            EXAMPLES::

                sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                sage: H = A.Harmonic()
                sage: for k in H.basis().keys() : print Permutation(k), factor(H.basis_element_to_multivariate_polynomial(k))
                [1, 2, 3] 2
                [1, 3, 2] (2) * (x1 - x2)
                [2, 1, 3] (2) * (x2 - x3)
                [2, 3, 1] (-1) * (-x1 + 2*x2 - x3) * (x1 - x3)
                [3, 1, 2] (-1) * (x2 - x3) * (-2*x1 + x2 + x3)
                [3, 2, 1] (-1) * (x2 - x3) * (-x1 + x2) * (x1 - x3)
            """
            R = self.polynomial_algebra()
            return harmonic_polynomial_for_permutation_group_element(R, key)

        def _repr_term(self, m):
            r"""
            Returns a string representation of a descent monomial
            indexed by a permutation group element.

            TESTS::

                sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                sage: H = A.Harmonic()
                sage: for k in H.basis().keys(): H._repr_term(k)
                'H_()'
                'H_(2,3)'
                'H_(1,2)'
                'H_(1,2,3)'
                'H_(1,3,2)'
                'H_(1,3)'
            """
            return 'H_' + str(m)

        def quotient_by_symmetric_polynomials(self):
            r"""
            EXAMPLES::

            """
            return self.realization_of().quotient_by_symmetric_polynomials().Harmonic()

    H = harmonic = Harmonic

    class Descent(CombinatorialFreeModule, BindableClass):
        r"""


        EXAMPLES::


        TESTS::


        """
        def __init__(self, AbsMod):
            r"""
            TESTS::

                sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                sage: A.Descent()
                Multivariate Polynomial in x1, x2, x3 over Symmetric Functions over Rational Field in the Schur basis expressed with descent monomials
            """
            CombinatorialFreeModule.__init__(self, SymmetricFunctions(AbsMod.base_ring()).schur(),
               SymmetricGroup(AbsMod.nb_variables()), prefix='Y_', latex_prefix='Y',
               latex_bracket=False, bracket=False, category=AbsMod.Bases())
            self._main_repr_var = AbsMod._main_repr_var

        def _repr_(self):
            r"""
            TESTS::


            """
            vars_string = join([self._main_repr_var+str(i) for i in range(1, self.nb_variables()+1)], ", ")
            return "Multivariate Polynomial in %s over %s expressed with descent monomials"%(vars_string, self.base_ring())

        def basis_element_to_multivariate_polynomial(self, key):
            r"""
            Returns the element basis of ``self`` indexed by the key
            ``key`` expressed in the polynomial algebra modelized by
            ``self``.

            EXAMPLES::

                sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                sage: D = A.Descent()
                sage: for k in D.basis().keys() : print Permutation(k),D.basis_element_to_multivariate_polynomial(k)
                [1, 2, 3] 1
                [1, 3, 2] x1*x3
                [2, 1, 3] x2
                [2, 3, 1] x2*x3
                [3, 1, 2] x3
                [3, 2, 1] x2*x3^2
            """
            R = self.polynomial_algebra()
            mon = Permutation(key).descent_polynomial()
            return R(mon)

        def _repr_term(self, m):
            r"""
            Returns a string representation of a descent monomial
            indexed by a permutation group element.

            TESTS::


            """
            return '(' + str(self.basis_element_to_multivariate_polynomial(m)) + ')'

        def quotient_by_symmetric_polynomials(self):
            r"""
            EXAMPLES::

            """
            return self.realization_of().quotient_by_symmetric_polynomials().Descent()

    D = descent = Descent

    class Staircase(CombinatorialFreeModule, BindableClass):
        r"""

        EXAMPLES::


        TESTS::


        """
        def __init__(self, AbsMod):
            r"""
            TESTS::

            """
            CombinatorialFreeModule.__init__(self, SymmetricFunctions(AbsMod.base_ring()).schur(),
               SymmetricGroup(AbsMod.nb_variables()), prefix='', latex_prefix='',
               latex_bracket=False, bracket=False, category=AbsMod.Bases())
            self._main_repr_var = AbsMod._main_repr_var

        def _repr_(self):
            r"""
            TESTS::

                sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                sage: A.Staircase()
                Multivariate Polynomial in x1, x2, x3 over Symmetric Functions over Rational Field in the Schur basis expressed with monmials under the staircase
            """
            vars_string = join([self._main_repr_var+str(i) for i in range(1, self.nb_variables()+1)], ", ")
            return "Multivariate Polynomial in %s over %s expressed with monmials under the staircase"%(vars_string, self.base_ring())

        def basis_element_to_multivariate_polynomial(self, key):
            r"""
            Returns the element basis of ``self`` indexed by the key
            ``key`` expressed in the polynomial algebra modelized by
            ``self``.

            EXAMPLES::


            """
            R = self.polynomial_algebra()
            exp = Permutation(key).to_lehmer_code()
            return R.prod([R.gens()[i]**exp[i] for i in range(self.nb_variables())])

        def _repr_term(self, m):
            r"""
            Returns a string representation of a monomial under the
            staircase indexed by a permutation group element.

            TESTS::

                sage: A = AbstractPolynomialRingAsSymModule(QQ, 3);
                sage: M = A.staircase()
                sage: for k in M.basis().keys(): print M._repr_term(k)
                (1)
                (x2)
                (x1)
                (x1*x2)
                (x1^2)
                (x1^2*x2)
            """
            return '(' + str(self.basis_element_to_multivariate_polynomial(m)) + ')'

        def quotient_by_symmetric_polynomials(self):
            r"""
            EXAMPLES::

            """
            return self.realization_of().quotient_by_symmetric_polynomials().Staircase()

    M = staircase = Staircase
